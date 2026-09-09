# 01_download_tcga_data.R
#
# Purpose: Download TCGA HNSCC (head-and-neck squamous cell carcinoma)
# multi-omic data via TCGAbiolinks: RNA-seq expression, 450K methylation,
# somatic mutations (MAF), and clinical data.
#

library(TCGAbiolinks)
library(SummarizedExperiment)

GDC_DATA_DIR <- "data/raw/tcga"
dir.create(GDC_DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# STEP 1: Clinical data (patient-level info: age, gender, stage, etc.)

clinical <- GDCquery_clinic(project = "TCGA-HNSC", type = "clinical")

# Quick sanity check
cat("Clinical records pulled:", nrow(clinical), "\n")

saveRDS(clinical, file.path("data/raw/tcga/clinical", "clinical_raw.rds"))

# STEP 2: RNA-seq gene expression

query_expr <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# This downloads the raw files into GDC_DATA_DIR 
GDCdownload(query_expr, directory = GDC_DATA_DIR)

expr_data <- GDCprepare(query_expr, directory = GDC_DATA_DIR)

saveRDS(expr_data, file.path("data/raw/tcga/expression", "expression_raw.rds"))
cat("Expression samples pulled:", ncol(expr_data), "\n")

# STEP 3: 450K methylation

query_methyl <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value"
)

GDCdownload(query_methyl, directory = GDC_DATA_DIR, method="api", files.per.chunk = 5)

methyl_data <- GDCprepare(query_methyl, directory = GDC_DATA_DIR)

saveRDS(methyl_data, file.path("data/raw/tcga/methylation", "methylation_raw.rds"))
cat("Methylation samples pulled:", ncol(methyl_data), "\n")

# STEP 4: Somatic mutations (MAF)

query_maf <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(
  query_maf,
  method = "api",
  files.per.chunk = 5,
  directory = GDC_DATA_DIR
)

maf_data <- GDCprepare(
  query_maf,
  directory = GDC_DATA_DIR
)

saveRDS(maf_data, file.path("data/raw/tcga/mutation", "mutation_raw.rds"))
cat("Mutation records pulled:", nrow(maf_data), "\n")

cat("\nDone. Raw .rds files saved under data/raw/tcga/{expression,methylation,mutation,clinical}/\n")
