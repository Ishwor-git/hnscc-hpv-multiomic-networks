# Distinct co-expression networks from multi-omic TCGA-HNSC data (HPV+ vs HPV−) — a replication of Costa, Boroni & Soares (2018)

This repository contains a from-scratch replication of the analysis described in:

> **Costa R.L., Boroni M. & Soares M.A.** — *"Distinct co-expression networks using multi-omic data reveal novel interventional targets in HPV-positive and negative head-and-neck squamous cell cancer."*
> *Scientific Reports* **8**, 15254 (2018). DOI: [10.1038/s41598-018-33498-5](https://doi.org/10.1038/s41598-018-33498-5)

The original paper is stored at [`tmp/s41598-018-33498-5.pdf`](tmp/s41598-018-33498-5.pdf).

The goal of the paper is to identify gene co-expression networks that differ between **HPV-positive (HPV+) and HPV-negative (HPV−)** head-and-neck squamous cell carcinoma (HNSCC), and to integrate multiple omic layers — gene expression, DNA promoter methylation, and somatic mutations — to highlight genes/pathways that could act as interventional targets. This project re-implements that pipeline on (nominally) the same TCGA-HNSC cohort, from raw data download to network reconstruction, enrichment analysis and external validation, and compares every step and every number against the published results.

**Documentation map:**

| Document | What's in it |
|---|---|
| `README.md` (this file) | What this is, repo layout, headline status |
| [`METHODOLOGY.md`](METHODOLOGY.md) | What the paper did, what was implemented step-by-step, which results reproduce the paper, and every deliberate/unintentional deviation with rationale |
| [`REPRODUCE.md`](REPRODUCE.md) | Environment setup, run order, runtime/disk caveats, expected output inventory |

---

## Repository structure

```
hnscc-hpv-multiomic-networks/
├── README.md                     <- this document
├── METHODOLOGY.md                <- paper vs. implementation, results, deviations
├── REPRODUCE.md                  <- environment + run instructions
├── hnscc-hpv-multiomic-networks.Rproj
├── data/
│   ├── raw/tcga/                 <- downloaded TCGA data: clinical, expression,
│   │                                methylation (450K), mutation (MAF), organised
│   │                                by GDC file-naming conventions
│   ├── raw/geo/                  <- external GEO datasets (GSE6791, GSE38266, GSE95036)
│   ├── external_refs/            <- empty (reserved; STRING/GO/etc. pulled on the fly)
│   └── processed/                <- intermediate R objects produced by scripts o2–o7
│                                    (sample_manifest.csv, *normalized.rds, wgcna_results.rds)
├── scripts/
│   ├── o1_download_tcga_data.R    <- download TCGA-HNSC multi-omic + clinical data
│   ├── o2_build_sample_manifest.R <- complete-case cohort (277 patients)
│   ├── o3_preprocess_exp.R        <- RNA-seq -> TMM log-CPM matrix (23,976 genes x 277)
│   ├── o4_preprocess_methylation.R<- 450K -> per-gene promoter beta matrix (19,160 genes x 277)
│   ├── o5_preprocess_mutations.R  <- MAF -> maftools MAF object + gene x patient matrix
│   ├── o6_deg_dmg_analysis.R      <- limma DEG/DMG, intersection + correlation (Fig. 2)
│   ├── o7_wgcma_modules.R         <- WGCNA modules + module–trait association (Fig. 3)
│   ├── o8_refine_network.R        <- HPV-status networks per module (Fig. 4) + metrics (Table 1)
│   ├── o9_enrichment_analysis.R   <- GO/KEGG/MSigDB enrichment (Table 2)
│   └── o10_external_valid.R       <- GEO validation experiments (Table 3)
├── results/
│   ├── tables/                    <- all quantitative outputs (CSV)
│   ├── networks/                  <- per-module HPV+ / HPV− co-expression graphs (graphml)
│   └── figures/                   <- generated figures (fig03 module-trait heatmap)
```

## Status at a glance

The full pipeline (`o1`–`o10`) runs end to end on the current TCGA-HNSC cohort and reproduces several of the paper's core biological conclusions:

- **Direction of differential expression** for every overlapping key gene (`ZFR2`, `SYCP2`, `CCNA1`, etc.) matches the paper.
- **Negative expression–methylation correlation** holds for all 10 doubly-selected DEG+DMG genes in this replication, in both HPV groups — the paper's central claim.
- **Module membership of key hub genes** (e.g. `SYCP2` in the Blue module) agrees with the paper wherever a module exists in both studies.
- **`SYCP2` as a Blue-module hub in HPV+**, densely connected and largely absent from HPV−, mirrors the paper's Fig. 4A.
- **Yellow-module keratinocyte/epidermis enrichment** matches the paper's Table 2 in both HPV groups.
- **HPV+ networks are consistently denser than HPV−**, reproducing the paper's most general conclusion.

It also diverges from the paper in several places — some deliberate (e.g. a global DEG/DMG test instead of per-stage tests, because HPV+ sample counts per stage are too small to power a per-stage design), some unintentional (e.g. no BMIQ normalisation in the methylation step, GEO validation blocked by inconsistent HPV-status column naming across series). See [`METHODOLOGY.md`](METHODOLOGY.md) for the full account of what matches, what doesn't, and why.

For environment setup and how to run the pipeline yourself, see [`REPRODUCE.md`](REPRODUCE.md).

---

## References

- Costa R.L., Boroni M. & Soares M.A. (2018). *Distinct co-expression networks using multi-omic data reveal novel interventional targets in HPV-positive and negative head-and-neck squamous cell cancer.* Scientific Reports 8:15254. [DOI:10.1038/s41598-018-33498-5](https://doi.org/10.1038/s41598-018-33498-5)
- Lawrence M.S. et al. (2015). *Comprehensive genomic characterization of head and neck squamous cell carcinomas.* Nature 517:576–582 (source TCGA-HNSC cohort).
- Langfelder P. & Horvath S. (2008/2005). WGCNA / Dynamic Tree Cut.
- Ritchie M.E. et al. (2015). limma. *Nucleic Acids Research* 43:e47.
- Yu G. et al. (2012/2016). clusterProfiler. *OMICS* 16:284–287.
- Original data availability from the paper: interactive networks and Neo4j HPV+ graph database at https://github.com/quelopes/HNSCC-network.
