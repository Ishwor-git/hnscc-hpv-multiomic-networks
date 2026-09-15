# Reproduction instructions

This document covers environment setup, run order, and what to expect (runtime, disk, output files) when reproducing this pipeline. For what each step does and how it compares to the original paper, see [`METHODOLOGY.md`](METHODOLOGY.md).

## Environment

- **R ≥ 4.x** (developed and verified with R 4.6.1).
- Packages (CRAN/Bioconductor):

```
TCGAbiolinks  SummarizedExperiment  edgeR  limma  WGCNA
igraph        Hmisc                 GEOquery         clusterProfiler
org.Hs.eg.db  maftools              httr    jsonlite
dplyr         tibble  tidyr         IlluminaHumanMethylation450kanno.ilmn12.hg19
impute
```

- Optional but recommended for full output:
  - `dorothea` — transcription-factor annotation (o8). Without it, `is_TF` is skipped.
  - `STRINGdb` — PPI edges (o8). Needs internet. Without it, PPI edges are skipped.
  - `msigdbr` — MSigDB Hallmark enrichment (o9). Without it only GO/KEGG run.

- **Internet access** is required for:
  - o1: TCGA GDC (TCGAbiolinks `GDCdownload`) and cBioPortal API
  - o8: STRING (if STRINGdb used)
  - o10: GEO (GEOquery)

## Run order

```bash
cd hnscc-hpv-multiomic-networks
Rscript scripts/o1_download_tcga_data.R    # download TCGA-HNSC (large)
Rscript scripts/o2_build_sample_manifest.R # -> data/processed/sample_manifest.csv
Rscript scripts/o3_preprocess_exp.R        # -> data/processed/expression_normalized.rds
Rscript scripts/o4_preprocess_methylation.R# -> data/processed/methylation_normalized.rds
Rscript scripts/o5_preprocess_mutations.R  # -> data/processed/mutation_normalized.rds
Rscript scripts/o6_deg_dmg_analysis.R      # -> results/tables/deg_list.csv, dmg_list.csv, deg_dmg_correlation.csv
Rscript scripts/o7_wgcma_modules.R         # -> results/tables/wgcna_*.csv, fig03 heatmap, data/processed/wgcna_results.rds
Rscript scripts/o8_refine_network.R        # -> results/networks/*.graphml, network_connection_metrics.csv
Rscript scripts/o9_enrichment_analysis.R   # -> results/tables/enrichment_results.csv
Rscript scripts/o10_external_valid.R       # -> results/tables/methylation_geo_validation.csv
```

Steps o3–o10 are order-dependent; each reads the processed R objects written by the previous steps. o1/o2 need only be re-run if the raw data should be re-downloaded.

## Expected inputs, outputs, and runtime caveats

- **Disk:** the raw TCGA download (o1) is on the order of several GB across the four omic layers, plus ~380 MB of GEO soft/series files for o10. The two `Wed_Sep__9_…*.tar.gz` archives in the repo root are leftover earlier-stage downloads and are not required by the pipeline.
- **Runtime:** the slowest steps are o4 (kNN imputation + probe→gene collapse over ~400k probes) and o7 (WGCNA on 8,000 genes), each typically taking several minutes on a multicore machine (`enableWGCNAThreads()` is called in o7). o10's GEO downloads dominate when run cold.
- **Reproducibility caveats:**
  - Cohort count (277 = 241 HPV− + 36 HPV+) reflects the state of the GDC and cBioPortal `hnsc_tcga_pub` records at download time; if those sources change, counts may shift slightly.
  - WGCNA auto-selects the soft-threshold power (fallback 6); a different power changes the module partition. Expect the *module colours* to differ across runs/environments even though the biology is stable.
  - Network edge counts depend on STRING availability (PPI edges) and on which DEG list precedes the network step.
  - `o10` uses a brittle auto-detection of HPV columns in GEO pheno data; inspect `pData(eset)` for GSE6791/GSE38266 if you need those validations (see [Unintentional gaps in METHODOLOGY.md](METHODOLOGY.md#unintentional-gaps-and-blockers)).

Expected output file inventory (as produced for this run):

```
results/tables/
  deg_dmg_correlation.csv           # 10 doubly-selected genes, rho per HPV status
  deg_list.csv                      # 628 DEG (FDR<=0.01, |logFC|>=2)
  dmg_list.csv                      #  39 DMG (FDR<=0.01, |logFC|>=2, M-values)
  enrichment_results.csv            # 1,162 FDR<=0.05 terms over the 6 networks
  methylation_geo_validation.csv    # 42 GSE95036 per-gene t-tests (none FDR-significant)
  network_connection_metrics.csv    # node/edge counts per module x status
  wgcna_module_assignments.csv      # 8,000 genes -> 8 modules
  wgcna_module_trait_correlation.csv
results/networks/                   # blue/yellow/turquoise x hpv_{pos,neg} graphml
results/figures/fig03_module_trait_heatmap/module_trait_heatmap.png
data/processed/*                    # intermediate RDS + sample_manifest.csv
```
