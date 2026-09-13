# 10_external_validation.R — see docs/10_external_validation.md

library(GEOquery)
library(limma)
library(dplyr)
library(tibble)

deg_list <- read.csv("results/tables/deg_list.csv")
dmg_list <- read.csv("results/tables/dmg_list.csv")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw/geo", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Validation 1: GSE6791 (Pyeon et al.) - microarray expression, HPV+ vs HPV-
# ---------------------------------------------------------------------------

validate_gse6791 <- function() {
  gset <- getGEO("GSE6791", GSEMatrix = TRUE, destdir = "data/raw/geo/GSE6791")
  eset <- gset[[1]]
  pheno <- pData(eset)

  hpv_col <- grep("hpv", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
  tissue_col <- grep("tissue|site|type", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]

  if (is.na(hpv_col)) {
    cat("GSE6791: could not auto-detect an HPV status column. Inspect pData(eset) manually.\n")
    return(NULL)
  }

  cat("GSE6791: using HPV column '", hpv_col, "'",
      if (!is.na(tissue_col)) paste0(", tissue column '", tissue_col, "'") else "", "\n", sep = "")

  hpv_status <- ifelse(grepl("positive|\\+", pheno[[hpv_col]], ignore.case = TRUE),
                        "HPV+",
                        ifelse(grepl("negative|-", pheno[[hpv_col]], ignore.case = TRUE),
                               "HPV-", NA))

  # Exclude cervical samples if a tissue/site column was found and mentions cervix
  keep_samples <- rep(TRUE, ncol(eset))
  if (!is.na(tissue_col)) {
    keep_samples <- !grepl("cervi", pheno[[tissue_col]], ignore.case = TRUE)
  }
  keep_samples <- keep_samples & !is.na(hpv_status)

  eset <- eset[, keep_samples]
  hpv_status <- hpv_status[keep_samples]

  cat("GSE6791: retained", ncol(eset), "samples (",
      sum(hpv_status == "HPV+"), "HPV+ /", sum(hpv_status == "HPV-"), "HPV- )\n")

  if (ncol(eset) < 6) {
    cat("GSE6791: too few samples after filtering - skipping DE test.\n")
    return(NULL)
  }

  group <- factor(hpv_status, levels = c("HPV-", "HPV+"))
  design <- model.matrix(~ group)
  fit <- eBayes(lmFit(exprs(eset), design))
  results <- topTable(fit, coef = "groupHPV+", number = Inf, adjust.method = "fdr") %>%
    rownames_to_column("probe_id")

  # Map probes to gene symbols via the platform's feature data, if present
  fdata <- fData(eset)
  symbol_col <- grep("symbol|gene.?name", colnames(fdata), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(symbol_col)) {
    results$gene_symbol <- fdata[[symbol_col]][match(results$probe_id, rownames(fdata))]
  } else {
    results$gene_symbol <- NA
    cat("GSE6791: no gene symbol column found in platform annotation - probe IDs only.\n")
  }

  sig_results <- results %>% filter(adj.P.Val <= 0.05, abs(logFC) >= 1)
  cat("GSE6791: significant DE genes (FDR<=0.05, |logFC|>=1):", nrow(sig_results), "\n")

  overlap <- intersect(sig_results$gene_symbol, deg_list$gene_symbol)
  cat("GSE6791: overlap with our TCGA DEG list:", length(overlap), "genes\n")
  print(overlap)

  write.csv(sig_results, "results/tables/gse6791_validation.csv", row.names = FALSE)
  sig_results
}

gse6791_results <- tryCatch(validate_gse6791(), error = function(e) {
  cat("GSE6791 validation failed:", conditionMessage(e), "\n")
  NULL
})

# ---------------------------------------------------------------------------
# Validation 2 & 3: GSE38266 (Lechner et al.) and GSE95036 (Esposti et al.)
# - 450K methylation, compare promoter methylation of key doubly-selected
# genes between HPV+ and HPV-
# ---------------------------------------------------------------------------

validate_methylation_geo <- function(gse_id, genes_of_interest) {
  gset <- getGEO(gse_id, GSEMatrix = TRUE, destdir = file.path("data/raw/geo", gse_id))
  eset <- gset[[1]]
  pheno <- pData(eset)
  beta <- exprs(eset)

  hpv_col <- grep("hpv", colnames(pheno), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(hpv_col)) {
    cat(gse_id, ": could not auto-detect an HPV status column. Inspect pData(eset) manually.\n")
    return(NULL)
  }
  hpv_status <- ifelse(grepl("positive|\\+", pheno[[hpv_col]], ignore.case = TRUE),
                        "HPV+",
                        ifelse(grepl("negative|-", pheno[[hpv_col]], ignore.case = TRUE),
                               "HPV-", NA))

  cat(gse_id, ": ", sum(hpv_status == "HPV+", na.rm = TRUE), "HPV+ / ",
      sum(hpv_status == "HPV-", na.rm = TRUE), "HPV- samples\n")

  fdata <- fData(eset)
  gene_col <- grep("UCSC_RefGene_Name|Gene_Symbol|symbol", colnames(fdata),
                    ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gene_col)) {
    cat(gse_id, ": no gene annotation column found - skipping.\n")
    return(NULL)
  }

  results <- lapply(genes_of_interest, function(gene) {
    probe_ids <- rownames(fdata)[grepl(paste0("(^|;)", gene, "(;|$)"), fdata[[gene_col]])]
    probe_ids <- intersect(probe_ids, rownames(beta))
    if (length(probe_ids) == 0) return(NULL)

    gene_beta <- colMeans(beta[probe_ids, , drop = FALSE], na.rm = TRUE)
    hpv_pos_vals <- gene_beta[hpv_status == "HPV+"]
    hpv_neg_vals <- gene_beta[hpv_status == "HPV-"]

    if (length(hpv_pos_vals) < 2 || length(hpv_neg_vals) < 2) return(NULL)

    test <- tryCatch(t.test(hpv_pos_vals, hpv_neg_vals), error = function(e) NULL)
    if (is.null(test)) return(NULL)

    data.frame(gene_symbol = gene, gse_id = gse_id,
               mean_beta_hpv_pos = mean(hpv_pos_vals),
               mean_beta_hpv_neg = mean(hpv_neg_vals),
               p_value = test$p.value)
  }) %>% bind_rows()

  results
}

genes_of_interest <- unique(c(dmg_list$gene_symbol, "SYCP2", "PITX2", "GJB6",
                               "HSF4", "MYO15B", "SERINC4", "CCNA1"))

gse38266_results <- tryCatch(validate_methylation_geo("GSE38266", genes_of_interest),
                              error = function(e) {
  cat("GSE38266 validation failed:", conditionMessage(e), "\n")
  NULL
})

gse95036_results <- tryCatch(validate_methylation_geo("GSE95036", genes_of_interest),
                              error = function(e) {
  cat("GSE95036 validation failed:", conditionMessage(e), "\n")
  NULL
})

methylation_validation <- bind_rows(gse38266_results, gse95036_results)
if (nrow(methylation_validation) > 0) {
  write.csv(methylation_validation, "results/tables/methylation_geo_validation.csv",
            row.names = FALSE)
  cat("\nSaved methylation validation results to results/tables/methylation_geo_validation.csv\n")
  print(methylation_validation)
} else {
  cat("\nNo methylation validation results produced - see per-dataset messages above.\n")
}
