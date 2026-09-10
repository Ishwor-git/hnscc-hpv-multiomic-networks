# 06_deg_dmg_analysis.R — see docs/06_deg_dmg_analysis.md

library(limma)
library(dplyr)
library(tibble)

expression_processed  <- readRDS("data/processed/expression_normalized.rds")
methylation_processed <- readRDS("data/processed/methylation_normalized.rds")

log_cpm          <- expression_processed$log_cpm
expr_sample_meta <- expression_processed$sample_meta
gene_methylation   <- methylation_processed$gene_methylation
methyl_sample_meta <- methylation_processed$sample_meta

# Ensembl -> symbol mapping for expression
expr_raw <- readRDS("data/raw/tcga/expression/expression_raw.rds")
gene_map <- as.data.frame(SummarizedExperiment::rowData(expr_raw)) %>%
  select(gene_id, gene_name) %>%
  distinct()
rm(expr_raw); gc()

symbol_lookup <- gene_map$gene_name[match(rownames(log_cpm), gene_map$gene_id)]

expr_by_symbol <- as.data.frame(log_cpm) %>%
  mutate(gene_symbol = symbol_lookup) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix()

# Align samples
common_patients <- intersect(expr_sample_meta$patient_barcode,
                              methyl_sample_meta$patient_barcode)

expr_cols   <- expr_sample_meta$sample_barcode[match(common_patients, expr_sample_meta$patient_barcode)]
methyl_cols <- methyl_sample_meta$sample_barcode[match(common_patients, methyl_sample_meta$patient_barcode)]

expr_aligned   <- expr_by_symbol[, expr_cols]
methyl_aligned <- gene_methylation[, methyl_cols]

hpv_group <- factor(
  expr_sample_meta$hpv_status[match(common_patients, expr_sample_meta$patient_barcode)],
  levels = c("HPV-", "HPV+")
)

# DEG
design <- model.matrix(~ hpv_group)
fit_expr <- eBayes(lmFit(expr_aligned, design))
deg_results <- topTable(fit_expr, coef = "hpv_groupHPV+", number = Inf,
                         adjust.method = "fdr") %>% rownames_to_column("gene_symbol")
deg_list <- deg_results %>% filter(adj.P.Val <= 0.01, abs(logFC) >= 2)

# DMG - beta values must be M-value transformed before limma
# (beta is bounded [0,1], so |logFC|>=2 can never trigger on raw beta;
# M = log2(beta/(1-beta)) puts methylation on a comparable, unbounded scale)
common_genes_for_dmg <- intersect(rownames(methyl_aligned), rownames(expr_aligned))
beta_for_fit <- methyl_aligned[common_genes_for_dmg, ]
beta_clipped <- pmin(pmax(beta_for_fit, 1e-4), 1 - 1e-4)  # avoid log(0)/log(Inf)
m_values <- log2(beta_clipped / (1 - beta_clipped))

fit_methyl <- eBayes(lmFit(m_values, design))
dmg_results <- topTable(fit_methyl, coef = "hpv_groupHPV+", number = Inf,
                         adjust.method = "fdr") %>% rownames_to_column("gene_symbol")
dmg_list <- dmg_results %>% filter(adj.P.Val <= 0.01, abs(logFC) >= 2)

# Intersect + correlate
doubly_selected <- intersect(deg_list$gene_symbol, dmg_list$gene_symbol)
cat("DEG:", nrow(deg_list), " DMG:", nrow(dmg_list),
    " Doubly-selected:", length(doubly_selected), "\n")

if (length(doubly_selected) > 0) {
  correlation_results <- lapply(doubly_selected, function(gene) {
    expr_vals   <- expr_aligned[gene, ]
    methyl_vals <- methyl_aligned[gene, ]  # correlate against beta, not M-value - more interpretable
    hpv_pos <- hpv_group == "HPV+"
    hpv_neg <- hpv_group == "HPV-"
    data.frame(
      gene_symbol = gene,
      rho_hpv_pos = cor(expr_vals[hpv_pos], methyl_vals[hpv_pos]),
      rho_hpv_neg = cor(expr_vals[hpv_neg], methyl_vals[hpv_neg])
    )
  }) %>% bind_rows()

  print(correlation_results)

  paper_genes <- c("SYCP2", "MEI1", "UGT8", "ZFR2", "SOX30", "FLRT3",
                    "PITX2", "SPRR2G", "GJB6", "MMP3", "CCNA1")
  print(correlation_results %>% filter(gene_symbol %in% paper_genes))
} else {
  cat("No doubly-selected genes - check DMG threshold before proceeding.\n")
  correlation_results <- data.frame()
}

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(deg_list, "results/tables/deg_list.csv", row.names = FALSE)
write.csv(dmg_list, "results/tables/dmg_list.csv", row.names = FALSE)
write.csv(correlation_results, "results/tables/deg_dmg_correlation.csv", row.names = FALSE)
