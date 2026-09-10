
library(SummarizedExperiment)
library(httr)
library(jsonlite)
library(dplyr)


# STEP 1: Load each raw object

clinical    <- readRDS("data/raw/tcga/clinical/clinical_raw.rds")
expr_data   <- readRDS("data/raw/tcga/expression/expression_raw.rds")
methyl_data <- readRDS("data/raw/tcga/methylation/methylation_raw.rds")
maf_data    <- readRDS("data/raw/tcga/mutation/mutation_raw.rds")

# STEP 2: Extract patient-level barcodes (first 12 characters) from each

patients_clinical <- unique(substr(clinical$submitter_id, 1, 12))

patients_expr <- unique(substr(colnames(expr_data), 1, 12))

patients_methyl <- unique(substr(colnames(methyl_data), 1, 12))

patients_mutation <- unique(substr(maf_data$Tumor_Sample_Barcode, 1, 12))

# Sanity check - print counts before merging
cat("Patients with clinical data:   ", length(patients_clinical), "\n")
cat("Patients with expression data: ", length(patients_expr), "\n")
cat("Patients with methylation data:", length(patients_methyl), "\n")
cat("Patients with mutation data:   ", length(patients_mutation), "\n")

# STEP 3: Pull HPV status from cBioPortal (same source confirmed earlier)

resp <- GET("https://www.cbioportal.org/api/studies/hnsc_tcga_pub/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY")
clin_cbio <- fromJSON(content(resp, "text", encoding = "UTF-8"))

hpv_status <- clin_cbio %>%
  filter(clinicalAttributeId == "HPV_STATUS") %>%
  select(patient_barcode = patientId, hpv_status = value) %>%
  distinct()

cat("Patients with HPV status from cBioPortal:", nrow(hpv_status), "\n")
cat("HPV status breakdown:\n")
print(table(hpv_status$hpv_status, useNA = "always"))

# STEP 4: Build the union of all patient IDs, then flag presence in each source with TRUE/FALSE columns

all_patients <- unique(c(patients_clinical, patients_expr,
                          patients_methyl, patients_mutation,
                          hpv_status$patient_barcode))

manifest <- data.frame(patient_barcode = sort(all_patients)) %>%
  mutate(
    has_clinical    = patient_barcode %in% patients_clinical,
    has_expression  = patient_barcode %in% patients_expr,
    has_methylation = patient_barcode %in% patients_methyl,
    has_mutation    = patient_barcode %in% patients_mutation
  ) %>%
  left_join(hpv_status, by = "patient_barcode")

# STEP 5: Define the analysis-ready cohort

manifest <- manifest %>%
  mutate(
    complete_case = has_expression & has_methylation & has_mutation &
                     hpv_status %in% c("HPV+", "HPV-")
  )

cat("\nTotal unique patients across all sources:", nrow(manifest), "\n")
cat("Patients with complete multi-omic + HPV data:",
    sum(manifest$complete_case), "\n")
cat("Of those, HPV+ :", sum(manifest$complete_case & manifest$hpv_status == "HPV+"), "\n")
cat("Of those, HPV- :", sum(manifest$complete_case & manifest$hpv_status == "HPV-"), "\n")

# STEP 6: Save

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write.csv(manifest, "data/processed/sample_manifest.csv", row.names = FALSE)

cat("\nManifest saved to data/processed/sample_manifest.csv\n")
