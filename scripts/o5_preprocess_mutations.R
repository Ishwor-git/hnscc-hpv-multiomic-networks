# 05_preprocess_mutation.R
#
# Purpose: Restrict the somatic mutation (MAF) data to the complete-case
# cohort, build a proper maftools MAF object (needed for Fisher's exact
# test comparisons later), and produce a simple gene x patient binary
# mutation matrix for quick lookups (e.g. TP53/CDKN2A status per patient,
# used when annotating the co-expression networks in a later script).
#
# This is a lighter preprocessing step than expression/methylation - MAF
# data doesn't need normalization, just filtering and restructuring.

library(maftools)
library(dplyr)

# ---------------------------------------------------------------------------
# STEP 1: Load manifest cohort + raw MAF data
# ---------------------------------------------------------------------------

manifest <- read.csv("data/processed/sample_manifest.csv", stringsAsFactors = FALSE)
cohort <- manifest %>% filter(complete_case)

maf_data <- readRDS("data/raw/tcga/mutation/mutation_raw.rds")

cat("Raw MAF rows (mutation records):", nrow(maf_data), "\n")

# ---------------------------------------------------------------------------
# STEP 2: Restrict to primary tumor samples in our cohort
# ---------------------------------------------------------------------------
# Same barcode logic as the other two preprocessing scripts - sample type
# "01" (primary tumor), matched against cohort patient IDs.

sample_type <- substr(maf_data$Tumor_Sample_Barcode, 14, 15)
patient_ids <- substr(maf_data$Tumor_Sample_Barcode, 1, 12)

is_primary_tumor <- sample_type == "01"
is_cohort_patient <- patient_ids %in% cohort$patient_barcode

maf_filtered <- maf_data[is_primary_tumor & is_cohort_patient, ]

cat("MAF rows after filtering to cohort primary tumors:",
    nrow(maf_filtered), "\n")
cat("Unique patients represented in filtered MAF:",
    length(unique(substr(maf_filtered$Tumor_Sample_Barcode, 1, 12))), "\n")

# Note: unlike expression/methylation, it's normal and expected for the
# patient count here to be LOWER than 277 - not every patient necessarily
# has a mutation record for every gene, and some patients may have zero
# reported mutations after MAF filtering (silent tumors, low purity, etc).
# This is different from a data-matching bug, so don't chase this number
# down to exactly 277 the way we did for expression/methylation.

# ---------------------------------------------------------------------------
# STEP 3: Build a proper maftools MAF object
# ---------------------------------------------------------------------------
# This gives access to maftools' built-in summary/plotting functions and
# is the expected input format for the differential mutation testing
# we'll do in a later script (Fisher's exact test, HPV+ vs HPV-).
#
# We also attach clinical data (HPV status) directly to the MAF object so
# maftools functions can group/color by it automatically.

clinical_for_maf <- cohort %>%
  transmute(Tumor_Sample_Barcode = patient_barcode, hpv_status)

# maftools matches on the patient-level barcode by default if the MAF's
# Tumor_Sample_Barcode column is trimmed the same way - we trim ours here
# for the clinical join, but keep the full barcode in the MAF data itself
# (maftools handles this internally via its own barcode truncation logic).

maf_obj <- read.maf(maf = maf_filtered, clinicalData = clinical_for_maf)

cat("\nmaftools MAF object summary:\n")
print(maf_obj)

# ---------------------------------------------------------------------------
# STEP 4: Build a simple gene x patient binary mutation matrix
# ---------------------------------------------------------------------------
# Useful later for quickly checking "is gene X mutated in patient Y"
# without querying the full MAF object each time - e.g. for annotating
# TP53/CDKN2A status onto the co-expression network visualizations.

mutation_binary <- maf_filtered %>%
  transmute(
    patient_barcode = substr(Tumor_Sample_Barcode, 1, 12),
    gene = Hugo_Symbol
  ) %>%
  distinct() %>%
  mutate(mutated = 1) %>%
  tidyr::pivot_wider(names_from = gene, values_from = mutated, values_fill = 0)

cat("Gene x patient mutation matrix:", nrow(mutation_binary), "patients x",
    ncol(mutation_binary) - 1, "genes\n")

# ---------------------------------------------------------------------------
# STEP 5: Save
# ---------------------------------------------------------------------------

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

mutation_processed <- list(
  maf_obj = maf_obj,
  mutation_binary = mutation_binary
)

saveRDS(mutation_processed, "data/processed/mutation_normalized.rds")

cat("\nSaved to data/processed/mutation_normalized.rds\n")

rm(maf_data, maf_filtered)
gc()
