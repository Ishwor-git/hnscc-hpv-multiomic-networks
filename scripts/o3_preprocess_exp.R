
library(SummarizedExperiment)
library(edgeR)
library(dplyr)


# STEP 1: Load the manifest and restrict to the complete-case cohort

manifest <- read.csv("data/processed/sample_manifest.csv", stringsAsFactors = FALSE)
cohort <- manifest %>% filter(complete_case)

cat("Cohort size from manifest:", nrow(cohort), "\n")
cat("HPV+ :", sum(cohort$hpv_status == "HPV+"), "\n")
cat("HPV- :", sum(cohort$hpv_status == "HPV-"), "\n")

# STEP 2: Load raw expression data

expr_data <- readRDS("data/raw/tcga/expression/expression_raw.rds")


# STEP 3: Restrict to primary tumor samples for our cohort patients

barcodes <- colnames(expr_data)
patient_ids <- substr(barcodes, 1, 12)
sample_type <- substr(barcodes, 14, 15)

is_primary_tumor <- sample_type == "01"
is_cohort_patient <- patient_ids %in% cohort$patient_barcode

keep_cols <- which(is_primary_tumor & is_cohort_patient)

cat("Primary tumor samples matching cohort:", length(keep_cols), "\n")

expr_data <- expr_data[, keep_cols]

# Deduplicate: if any patient still has >1 sample after the above filter,
# keep only the first occurrence.
dedup_patient_ids <- substr(colnames(expr_data), 1, 12)
first_occurrence <- !duplicated(dedup_patient_ids)
expr_data <- expr_data[, first_occurrence]

cat("Samples after deduplication:", ncol(expr_data), "\n")


# STEP 4: Extract the raw count matrix and remove zero-variance genes


counts <- assay(expr_data, "unstranded")  # raw counts; STAR-Counts output
# has multiple count columns (unstranded/stranded_first/stranded_second) -
# "unstranded" is the standard choice unless your library prep was
# specifically stranded, which HNSC TCGA samples generally were not.

gene_sd <- apply(counts, 1, sd)
counts <- counts[gene_sd > 0, ]

cat("Genes remaining after removing zero-variance genes:", nrow(counts), "\n")


# STEP 5: CPM normalization + log transform

dge <- DGEList(counts = counts)

# Build a group vector aligned to the (deduplicated, filtered) sample
# columns, needed for filterByExpr's group-aware filtering.
sample_patient_ids <- substr(colnames(counts), 1, 12)
group <- cohort$hpv_status[match(sample_patient_ids, cohort$patient_barcode)]

keep_genes <- filterByExpr(dge, group = group)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

cat("Genes remaining after filterByExpr:", nrow(dge), "\n")

dge <- calcNormFactors(dge)  # TMM normalization - accounts for
# compositional differences between libraries, standard before CPM

# log2(CPM + 1)-equivalent transform (edgeR's prior.count handles the +1
# smoothing internally, avoiding log(0) issues)
log_cpm <- cpm(dge, log = TRUE, prior.count = 1)

cat("Final expression matrix dimensions:", nrow(log_cpm), "genes x",
    ncol(log_cpm), "samples\n")


# STEP 6: Attach sample metadata (HPV status, stage, etc.) and save


sample_meta <- data.frame(
  sample_barcode = colnames(log_cpm),
  patient_barcode = substr(colnames(log_cpm), 1, 12)
) %>%
  left_join(cohort, by = "patient_barcode")

# Sanity check - this should be TRUE. If not, something went wrong with
# the ordering/matching above and needs fixing before moving on.
stopifnot(all(sample_meta$patient_barcode == substr(colnames(log_cpm), 1, 12)))

expression_processed <- list(
  log_cpm = log_cpm,
  sample_meta = sample_meta
)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(expression_processed, "data/processed/expression_normalized.rds")

cat("\nSaved to data/processed/expression_normalized.rds\n")

# Free the large raw object now that we're done with it
rm(expr_data, counts, dge)
gc()
