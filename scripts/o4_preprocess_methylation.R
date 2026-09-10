
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(impute)
library(dplyr)
library(tidyr)

# STEP 1: Load manifest cohort + raw methylation data

manifest <- read.csv("data/processed/sample_manifest.csv", stringsAsFactors = FALSE)
cohort <- manifest %>% filter(complete_case)

methylation <- readRDS("data/raw/tcga/methylation/methylation_raw.rds")


# STEP 2: Restrict to primary tumor samples in our cohort (same logic as
# the expression script - sample type "01", dedupe to one per patient)


barcodes <- colnames(methylation)
patient_ids <- substr(barcodes, 1, 12)
sample_type <- substr(barcodes, 14, 15)

is_primary_tumor <- sample_type == "01"
is_cohort_patient <- patient_ids %in% cohort$patient_barcode
keep_cols <- which(is_primary_tumor & is_cohort_patient)

cat("Primary tumor samples matching cohort:", length(keep_cols), "\n")

methylation <- methylation[, keep_cols]

dedup_patient_ids <- substr(colnames(methylation), 1, 12)
methylation <- methylation[, !duplicated(dedup_patient_ids)]

cat("Samples after deduplication:", ncol(methylation), "\n")

beta <- assay(methylation)  # probes x samples beta-value matrix


# STEP 3: Probe-level QC filtering


ann450k <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

# Restrict annotation to probes actually present in our beta matrix
ann450k <- ann450k[rownames(ann450k) %in% rownames(beta), ]

# Remove sex chromosome probes
autosomal_probes <- rownames(ann450k)[!(ann450k$chr %in% c("chrX", "chrY"))]
beta <- beta[rownames(beta) %in% autosomal_probes, ]

cat("Probes remaining after removing chrX/chrY:", nrow(beta), "\n")

# Remove probes with more than 5% missing values across samples
missing_frac <- rowMeans(is.na(beta))
beta <- beta[missing_frac <= 0.05, ]

cat("Probes remaining after >5% missing filter:", nrow(beta), "\n")

# STEP 4: Impute remaining missing values (kNN, k=10, matching the paper)

# impute.knn expects a plain numeric matrix and can be memory-heavy on
# large inputs - this step may take a few minutes.

if (any(is.na(beta))) {
  imputed <- impute.knn(as.matrix(beta), k = 10)
  beta <- imputed$data
} else {
  cat("No missing values remained - skipping imputation.\n")
}

# STEP 5: Map probes to gene promoter regions (TSS200 > 1stExon > TSS1500)

ann_subset <- ann450k[rownames(ann450k) %in% rownames(beta),
                       c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]
ann_subset$probe_id <- rownames(ann_subset)

# Expand semicolon-delimited multi-gene annotations into long format:
# one row per (probe, gene, region) combination.
probe_gene_map <- ann_subset %>%
  filter(UCSC_RefGene_Name != "") %>%
  mutate(
    gene_list   = strsplit(UCSC_RefGene_Name, ";"),
    region_list = strsplit(UCSC_RefGene_Group, ";")
  ) %>%
  rowwise() %>%
  mutate(pairs = list(data.frame(gene = gene_list, region = region_list))) %>%
  ungroup() %>%
  select(probe_id, pairs) %>%
  tidyr::unnest(pairs) %>%
  distinct()

cat("Probe-gene-region annotations expanded:", nrow(probe_gene_map), "rows\n")

# For each gene, apply the priority: TSS200 > 1stExon > TSS1500
assign_promoter_probes <- function(gene_df) {
  for (priority_region in c("TSS200", "1stExon", "TSS1500")) {
    matched <- gene_df$probe_id[gene_df$region == priority_region]
    if (length(matched) > 0) {
      return(data.frame(probe_id = matched, region_used = priority_region))
    }
  }
  return(data.frame(probe_id = character(0), region_used = character(0)))
}

gene_promoter_probes <- probe_gene_map %>%
  group_by(gene) %>%
  group_modify(~ assign_promoter_probes(.x)) %>%
  ungroup()

cat("Genes with an assigned promoter probe set:",
    length(unique(gene_promoter_probes$gene)), "\n")

# STEP 6: Average beta values per gene's promoter probe set -> gene-level
# promoter methylation matrix (genes x samples)

gene_methylation <- gene_promoter_probes %>%
  left_join(
    as.data.frame(beta) %>%
      tibble::rownames_to_column("probe_id"),
    by = "probe_id"
  ) %>%
  select(-region_used, -probe_id) %>%
  group_by(gene) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

cat("Final gene-level promoter methylation matrix:",
    nrow(gene_methylation), "genes x", ncol(gene_methylation), "samples\n")

# STEP 7: Attach sample metadata and save

sample_meta <- data.frame(
  sample_barcode = colnames(gene_methylation),
  patient_barcode = substr(colnames(gene_methylation), 1, 12)
) %>%
  left_join(cohort, by = "patient_barcode")

stopifnot(all(sample_meta$patient_barcode ==
                substr(colnames(gene_methylation), 1, 12)))

methylation_processed <- list(
  gene_methylation = gene_methylation,
  sample_meta = sample_meta
)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(methylation_processed, "data/processed/methylation_normalized.rds")

cat("\nSaved to data/processed/methylation_normalized.rds\n")

rm(methylation, beta, ann450k, ann_subset, probe_gene_map, gene_promoter_probes)
gc()
