# 07_wgcna_modules.R — see docs/07_wgcna_modules.md

library(WGCNA)
library(dplyr)
library(tibble)

enableWGCNAThreads()

expression_processed <- readRDS("data/processed/expression_normalized.rds")
log_cpm     <- expression_processed$log_cpm
sample_meta <- expression_processed$sample_meta

# Top 8000 most variant genes by MAD
gene_mad <- apply(log_cpm, 1, mad)
top_genes <- names(sort(gene_mad, decreasing = TRUE))[1:8000]
datExpr <- t(log_cpm[top_genes, ])  # WGCNA wants samples x genes

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
cat("Expression matrix for WGCNA:", nrow(datExpr), "samples x",
    ncol(datExpr), "genes\n")

# Soft-thresholding power selection
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 3)
softPower <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)
cat("Soft threshold power selected:", softPower, "\n")

# Module detection
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = 20,
  reassignThreshold = 0,
  mergeCutHeight = 0.45,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  maxBlockSize = ncol(datExpr) + 100,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)
cat("Modules detected:", length(unique(moduleColors)), "\n")
print(table(moduleColors))

MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)

# Build trait data aligned to datExpr sample order
sample_order <- rownames(datExpr)
meta_aligned <- sample_meta[match(sample_order, sample_meta$sample_barcode), ]

trait_data <- data.frame(
  hpv_status = ifelse(meta_aligned$hpv_status == "HPV+", 1, 0)
)

# Best-effort: pull extra clinical traits if available, skip gracefully if not
clinical <- tryCatch(readRDS("data/raw/tcga/clinical/clinical_raw.rds"), error = function(e) NULL)
if (!is.null(clinical)) {
  clinical_aligned <- clinical[match(meta_aligned$patient_barcode,
                                      substr(clinical$submitter_id, 1, 12)), ]
  candidate_traits <- c("age_at_index", "gender", "ajcc_pathologic_stage")
  for (tr in intersect(candidate_traits, colnames(clinical_aligned))) {
    trait_data[[tr]] <- tryCatch(
      as.numeric(as.factor(clinical_aligned[[tr]])),
      error = function(e) NULL
    )
  }
}
cat("Traits used for module correlation:", paste(colnames(trait_data), collapse = ", "), "\n")

moduleTraitCor <- cor(MEs, trait_data, use = "p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

cat("\nModule-trait correlation with HPV status:\n")
print(round(data.frame(module = rownames(moduleTraitCor),
                        cor_hpv = moduleTraitCor[, "hpv_status"],
                        p_hpv = moduleTraitP[, "hpv_status"]), 4))

sig_modules <- rownames(moduleTraitCor)[abs(moduleTraitCor[, "hpv_status"]) >= 0.25 &
                                          moduleTraitP[, "hpv_status"] <= 0.001]
cat("\nModules significantly associated with HPV status (|cor|>=0.25, p<=0.001):\n")
print(sig_modules)

# Save heatmap
dir.create("results/figures/fig03_module_trait_heatmap", recursive = TRUE, showWarnings = FALSE)
png("results/figures/fig03_module_trait_heatmap/module_trait_heatmap.png",
    width = 900, height = 900)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(trait_data),
  yLabels = rownames(moduleTraitCor),
  ySymbols = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = paste0(round(moduleTraitCor, 2), "\n(", round(moduleTraitP, 3), ")"),
  setStdMargins = FALSE,
  cex.text = 0.6,
  main = "Module-trait relationships"
)
dev.off()

# Save outputs
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
module_assignments <- data.frame(gene = colnames(datExpr), module = moduleColors)
write.csv(module_assignments, "results/tables/wgcna_module_assignments.csv", row.names = FALSE)
write.csv(as.data.frame(moduleTraitCor) %>% rownames_to_column("module"),
          "results/tables/wgcna_module_trait_correlation.csv", row.names = FALSE)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(net = net, moduleColors = moduleColors, MEs = MEs,
             datExpr = datExpr, trait_data = trait_data,
             moduleTraitCor = moduleTraitCor, moduleTraitP = moduleTraitP,
             sig_modules = sig_modules),
        "data/processed/wgcna_results.rds")

cat("\nSaved module assignments, correlation table, heatmap, and WGCNA object.\n")
