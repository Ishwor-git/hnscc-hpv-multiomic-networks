# Distinct co-expression networks from multi-omic TCGA-HNSC data (HPV+ vs HPV−) — a replication of Costa, Boroni & Soares (2018)

This repository contains a from-scratch replication of the analysis described in:

> **Costa R.L., Boroni M. & Soares M.A.** — *"Distinct co-expression networks using multi-omic data reveal novel interventional targets in HPV-positive and negative head-and-neck squamous cell cancer."*
> *Scientific Reports* **8**, 15254 (2018). DOI: [10.1038/s41598-018-33498-5](https://doi.org/10.1038/s41598-018-33498-5)

The original paper is stored at [`tmp/s41598-018-33498-5.pdf`](tmp/s41598-018-33498-5.pdf).

The goal of the paper is to identify gene co-expression networks that differ between **HPV-positive (HPV+) and HPV-negative (HPV−)** head-and-neck squamous cell carcinoma (HNSCC), and to integrate multiple omic layers — gene expression, DNA promoter methylation, and somatic mutations — to highlight genes/pathways that could act as interventional targets. This project re-implements that pipeline on (nominally) the same TCGA-HNSC cohort, from raw data download to network reconstruction, enrichment analysis and external validation, and compares every step and every number against the published results.

---

## Table of contents

1. [Repository structure](#repository-structure)
2. [What the paper did](#what-the-paper-did)
3. [What was implemented](#what-was-implemented)
4. [Step-by-step comparison with the paper](#step-by-step-comparison-with-the-paper)
5. [Results that reproduce the paper](#results-that-reproduce-the-paper)
6. [Deliberate differences and their rationale](#deliberate-differences-and-their-rationale)
7. [Unintentional gaps and blockers](#unintentional-gaps-and-blockers)
8. [Comparison of achieved results vs the paper](#comparison-of-achieved-results-vs-the-paper)
9. [Reproduction instructions](#reproduction-instructions)
10. [References](#references)

---

## Repository structure

```
hnscc-hpv-multiomic-networks/
├── README.md                     <- this document
├── MANIFEST.txt                  <- legacy tcga data manifest (GDC helper output)
├── hnscc-hpv-multiomic-networks.Rproj
├── tmp/
│   └── s41598-018-33498-5.pdf    <- original paper
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
├── docs/  notebooks/  R/          <- empty structure reserved for future work
└── df.rds, results.rds            <- exploratory artifacts from data familiarisation;
                                      NOT used by the pipeline
```

Notes on layout:

- Several script headers (o6–o10) say *"see docs/06_…md"* etc. Those `docs/*.md` files were never written; **this README is the authoritative documentation**.
- `o7_wgcma_modules.R` is a filename typo for the WGCNA step (header says `07_wgcna_modules.R`). Keep the name as-is for numbering continuity.
- The two large `Wed_Sep__9_…*.tar.gz` archives at the repository root are legacy downloads of TCGA files from an earlier run; the pipeline reads the already-extracted files under `data/raw/tcga`, so they can be deleted (see *Reproduction instructions*).

---

## What the paper did

The 2018 paper proceeds through five main analyses, all on the TCGA-HNSC cohort defined in the 2015 TCGA HNSCC characterization (Lawrence et al., *Nature* 517:576):

1. **Differential expression and methylation (DEG/DMG).** Using limma moderated t-tests, HPV+ vs HPV− samples were compared *within each disease stage* (stage I vs I, II vs II, III vs III, IV vs IV). Genes significant in **at least one stage** were kept: 223 DEG (`|logFC| ≥ 4`, FDR ≤ 0.01) and 359 DMG (`|logFC| ≥ 2`, FDR ≤ 0.01). The intersection (14 genes; 12 retained after excluding `PNLDC1`/`CNTN1`) was checked for negative expression–methylation correlation, supporting that promoter methylation represses expression.
2. **WGCNA.** On the 8,000 most variant genes (median absolute deviation), blockwise WGCNA (`minModuleSize = 20`, `mergeCutHeight = 0.45`, soft-threshold power `β = 4`) identified **7 modules** (Supplementary Fig. 2). Module eigengenes were correlated with 7 clinical traits. Three modules — **Blue, Yellow, Grey** — passed `|cor| ≥ 0.25` and `p ≤ 0.001` for HPV status.
3. **Network refinement.** For each significant module, samples were split into HPV+ / HPV−; within each group, Spearman correlations were computed between all module genes; edges were kept at `rho ≥ 0.65`, `p ≤ 0.01`; a node was kept if it was a **DEG or directly connected to a DEG**. Signed off against transcription factors (TFcheckpoint), doubly DEG+DMG genes, DMG, and significantly mutated genes linked by STRING PPI (confidence ≥ 0.7). The paper found HPV+ networks consistently **more densely connected** than HPV− (Table 1).
4. **Enrichment.** GO Biological Process, KEGG, and MSigDB (v6.0) over-representation, hypergeometric tests with FDR adjustment (Table 2).
5. **External validation.** GSE6791 (microarray expression, Pyeon et al.), GSE38266 (450K methylation, Lechner et al.), GSE95036 (450K methylation FFPE, Esposti et al.). Methylation datasets processed with the same promoter-mapping methodology as the TCGA data (`|logFC| ≥ 1.5`, FDR ≤ 0.05).

### Paper's headline numbers

| Output | Paper value |
|---|---|
| Cohort | 240 HPV− + 36 HPV+ = 276 tumors analysed |
| Expression genes | 20,502 |
| Methylation genes | 14,861 |
| DEG / DMG | 223 / 359 (per-stage union) |
| Doubly selected DEG+DMG | 14 (12 after excluding `PNLDC1`, `CNTN1`); all negatively correlated |
| WGCNA | top 8,000 MAD genes, `β = 4`, 7 modules |
| Significant modules | Blue, Yellow, Grey (`|cor| ≥ 0.25`, `p ≤ 0.001`) |
| Blue network | HPV+ 539 nodes / 2,633 edges; HPV− 112 / 114 |
| Yellow network | HPV+ 145 / 640; HPV− 36 / 34 |
| Grey network | HPV+ 127 / 111; HPV− 8 / 7 |
| Validation | GSE6791: 38 concordant DE genes; GSE38266 + GSE95036: DMG congruence for `SYCP2`, `PITX2`, `GJB6`, `HSF4`, `MYO15B`, `SERINC4` |

---

## What was implemented

The pipeline is fully reproduced as 10 numbered R scripts (`scripts/o1_…R` … `o10_external_valid.R`), which must be run in order.

### o1 — download TCGA data
TCGAbiolinks is used to pull TCGA-HNSC:
- clinical (patient-level) data via `GDCquery_clinic`
- gene expression quantification (workflow `STAR - Counts`)
- 450K methylation beta values (`Illumina Human Methylation 450`)
- somatic mutations (`Masked Somatic Mutation`, aliquot-ensemble MAF)
- HPV-status clinical attribute cross-checked from cBioPortal (`hnsc_tcga_pub`)

`GDCdownload`/`GDCprepare` save raw `SummarizedExperiment`/data-frame objects under `data/raw/tcga/{clinical,expression,methylation,mutation}/`.

### o2 — build the sample manifest
Constructs a cohort from the union of all patients across the four omic sources plus cBioPortal HPV status. Each patient is flagged per-omic; a **complete case** requires expression + methylation + mutation data plus a known HPV status (`HPV+`/`HPV-`). At the time of this run this produced **277 complete cases (241 HPV−, 36 HPV+)**, stored in `data/processed/sample_manifest.csv`.

### o3 — expression preprocessing
- Restrict to **primary tumor** samples (barcode sample-type `01`), dedupe to one sample per patient.
- Remove zero-variance genes (`sd == 0`).
- `edgeR::filterByExpr` (group-aware, using HPV status) → TMM normalization (`calcNormFactors`) → `log2(CPM + 1)` via `cpm(log = TRUE, prior.count = 1)`.
- Collapse Ensembl gene IDs to unique symbols (mean of duplicate entries).
- Result: **23,976 genes × 277 samples** → `data/processed/expression_normalized.rds`.

### o4 — methylation preprocessing
- Restrict to primary tumors, dedupe (same barcode logic as expression).
- Probe QC: drop **sex-chromosome probes** (chrX/chrY) and probes with **>5% missing** values across samples; missing values imputed with **kNN** (`impute::impute.knn`, `k = 10`) exactly as in the paper.
- Probe → gene promoter mapping exactly as Jiao et al.: for each gene use TSS200 probes; if none, 1st-exon probes; if none, TSS1500 probes.
- Average beta values across a gene's assigned promoter probes.
- Result: **19,160 genes × 277 samples** → `data/processed/methylation_normalized.rds`.

> **Note.** Unlike the paper, this step does **not** apply BMIQ normalization or the cross-reactive / non-specific / SNP / zero-coordinate probe black-list filters (see [Unintentional gaps](#unintentional-gaps-and-blockers)).

### o5 — mutation preprocessing
- Filter the MAF to cohort primary-tumor samples.
- Build a `maftools::read.maf` object with clinical (HPV) data attached.
- Build a gene × patient binary mutation matrix.
- → `data/processed/mutation_normalized.rds`.

### o6 — DEG/DMG analysis (paper Fig. 2)
- Align expression and methylation to the common patient set.
- **DEG**: one limma design `~ hpv_group` (levels `HPV-`, `HPV+`); empirical-Bayes moderated test; keep `adj.P.Val ≤ 0.01` and `|logFC| ≥ 2` → **628 DEG**.
- **DMG**: beta values are transformed to **M values** (`log2(beta/(1-beta))`, clipped to avoid `log(0)`) before limma — necessary because `|logFC| ≥ 2` is meaningless on the bounded beta scale; same thresholds → **39 DMG**.
- Intersect → **10 doubly-selected genes**; for each, the **Pearson correlation** between expression and promoter beta is computed separately in HPV+ and HPV− samples.
- Outputs: `results/tables/deg_list.csv`, `dmg_list.csv`, `deg_dmg_correlation.csv`.

### o7 — WGCNA (paper Fig. 3 + Supplementary Figs. 2–3)
- Take the **top 8,000 most variant genes by MAD** from the log-CPM matrix (`t()`, samples × genes as WGCNA expects).
- `goodSamplesGenes` QC; `pickSoftThreshold(powerVector = 1:20)`; the selected power is piped into `blockwiseModules` with `TOMType = "unsigned"`, `minModuleSize = 20`, `reassignThreshold = 0`, `mergeCutHeight = 0.45`, single block (`maxBlockSize = ncol + 100`).
  - *If the soft-threshold estimate is `NA`, the script falls back to `power = 6`.* (The paper fixed `β = 4`.)
- Module eigengenes (`orderMEs`), correlated against the available clinical traits with `corPvalueStudent`.
- Significant modules: `|cor| ≥ 0.25` and `p ≤ 0.001` — **same thresholds as the paper**.
- Outputs: `results/tables/wgcna_module_assignments.csv`, `wgcna_module_trait_correlation.csv`, `results/figures/fig03_module_trait_heatmap/module_trait_heatmap.png`, `data/processed/wgcna_results.rds`.

### o8 — network refinement (paper Fig. 4 + Table 1)
For each significant module × HPV status:
- Map Ensembl IDs → symbols (keep the highest-MAD duplicate when a symbol maps from several Ensembl IDs).
- Restrict to the module's genes; compute the full **Spearman correlation** matrix among module genes in that HPV subgroup (`Hmisc::rcorr`).
- Keep edges with `|rho| ≥ 0.65` and `p ≤ 0.01` (the paper's thresholds).
- Keep only **DEG nodes plus nodes directly connected to a DEG**.
- Add **STRING PPI edges** (confidence ≥ 700, STRING DB v12.0) among the retained genes.
- Annotate nodes: `is_DEG`, `is_DMG`, `is_TF` (from **dorothea** high-confidence regulons, A/B — a substitute for TFcheckpoint), `is_significantly_mutated` (Fisher's exact test of mutation status vs HPV status, raw `p ≤ 0.05`).
- Write GraphML + connection metrics.

Outputs: `results/networks/{blue,yellow,turquoise}_{hpv_pos,hpv_neg}.graphml`, `results/tables/network_connection_metrics.csv`.

### o9 — enrichment (paper Table 2)
For each refined network (node list):
- `clusterProfiler::enrichGO` (Biological Process), `enrichKEGG` (`hsa`), and `enricher` against the **MSigDB Hallmark** gene sets (`msigdbr`, category `H`).
- Keep `p.adjust ≤ 0.05`.

Output: `results/tables/enrichment_results.csv` (1,162 significant terms across the six networks).

### o10 — external validation (paper Table 3 + Supplementary Figs. 6–7)
- **GSE6791** (Pyeon et al. microarray): GEOquery → limma moderated t-test of HPV+ vs HPV− (`FDR ≤ 0.05`, `|logFC| ≥ 1`), cervical samples excluded, probes mapped to symbols via platform annotation.
- **GSE38266** (Lechner et al.) and **GSE95036** (Esposti et al.) 450K methylation: per-gene promoter-beta t-test in HPV+ vs HPV− for all DMG genes plus the paper's key genes.

Output: `results/tables/methylation_geo_validation.csv` (see [Caveats in the validation step](#caveats-in-the-validation-step)).

---

## Step-by-step comparison with the paper

| Step | Paper (Costa et al. 2018) | This replication | Match? |
|---|---|---|---|
| **Cohort** | TCGA-HNSC, 279 primary tumours; 3 removed (no stage) → **240 HPV− + 36 HPV+ = 276** | TCGA-HNSC, cBioPortal HPV status; **241 HPV− + 36 HPV+ = 277** complete cases | ✓ essentially identical |
| **HPV status source** | 2015 TCGA HNSC characterization (Lawrence et al.) | cBioPortal `hnsc_tcga_pub` (same underlying study) | mostly ✓ |
| **Expression data** | Illumina HiSeq 2000, level-3 **RNAseqV2 (RSEM)** counts | current GDC **STAR-Counts** raw unstranded counts | ✗ different pipeline |
| **Expression norm.** | `log(1 + p)` on RSEM counts; remove zero-SD genes | `edgeR::filterByExpr` → **TMM** → `log2(CPM+1)`; remove zero-variance; mean-collapse duplicates | ✗ |
| **Expression genes** | 20,502 | 23,976 | ✗ |
| **Methylation QC** | remove cross-reactive, non-specific, SNP, chrX/Y, zero-coordinate probes; >5% missing; **BMIQ** (ChAMP) | remove chrX/Y and >5% missing only; kNN impute; **no BMIQ** | ✗ partial |
| **Methylation genes** | 14,861 | 19,160 | ✗ |
| **Promoter mapping** | average of TSS200 probes → 1st exon → TSS1500 (Jiao et al.) | same priority rule | ✓ |
| **DEG test** | limma, **per-stage** (HPV+ vs HPV− within stages I–IV), **`|logFC| ≥ 4`**, FDR ≤ 0.01, union over stages | limma, **one global HPV+ vs HPV− test**, **`|logFC| ≥ 2`**, FDR ≤ 0.01 | ✗ |
| **DEG count** | **223** | **628** | ✗ |
| **DMG test** | limma on M values, **per-stage**, `|logFC| ≥ 2`, FDR ≤ 0.01, union over stages | limma on M values, **global**, `|logFC| ≥ 2`, FDR ≤ 0.01 | ~ (same thresholds, different design) |
| **DMG count** | **359** | **39** | ✗ |
| **Doubly selected** | 14 (12 kept) | 10 | ✗ list differs, ✓ all negative correlations |
| **WGCNA input** | top **8,000** MAD-variant genes | top **8,000** MAD-variant genes | ✓ |
| **Soft-threshold** | fixed `β = 4` | auto `pickSoftThreshold` (fallback 6) | ✗ |
| **WGCNA params** | `minModuleSize = 20`, `mergeCutHeight = 0.45` | same | ✓ |
| **Modules** | **7** | **8** (black, blue, brown, green, grey, red, turquoise, yellow) | ✗ |
| **Traits correlated** | HPV status, staging, age, gender, alcohol, smoked, anatomical site (7) | HPV status, age, AJCC pathologic stage (3) | ✗ |
| **Sig. modules** (`|cor|≥0.25, p≤0.001`) | **Blue, Yellow, Grey** | **Yellow, Turquoise, Blue** (Grey: `−0.24`, just below the 0.25 cut-off) | ✗ |
| **Network rule** | Spearman within module×status; keep `rho ≥ 0.65`, `p ≤ 0.01`; keep DEG + DEG-connected nodes; add TF/DMG/mutated PPI | identical core rule | ✓ |
| **Network scale** | Blue 539/2633 (HPV+) vs 112/114 (HPV−); Yellow 145/640 vs 36/34; Grey 127/111 vs 8/7 | see [comparison tables](#comparison-of-achieved-results-vs-the-paper) | ✗ counts differ (but direction reproduced) |
| **TF source** | TFcheckpoint database | dorothea (confidence A/B regulons) | ✗ substitute |
| **PPI source** | STRING v9.1, confidence ≥ 0.7 | STRING v12.0, `score_threshold = 700` | ~ same semantics, newer build |
| **Mutated-gene test** | `maftools` Fisher's exact test, **FDR-adjusted** ≤ 0.05 | Fisher's exact test, **raw** p ≤ 0.05 (in o8) | ✗ |
| **Enrichment DBs** | GO BP, KEGG, **MSigDB v6.0 (incl. C6)** | GO BP, KEGG, **MSigDB Hallmark (H) only** | ✗ C6 omitted |
| **Enrichment FDR** | ≤ 0.05 (methods text says ≤ 0.001 — internal inconsistency in the paper) | ≤ 0.05 | ✓ |
| **Validation GSE6791** | 56 HNSCC samples, 16 HPV+; MAS5 + limma; 38 concordant DEG (Table 3) | **not produced** — GEO series matrix carries no HPV annotation that the script auto-detects | ✗ blocked |
| **Validation GSE38266** | 11 samples (6 HPV+), fresh frozen; DMG congruence | **not produced** — HPV status embedded in value strings, not detected | ✗ blocked |
| **Validation GSE95036** | 42 samples (21 HPV+), FFPE; DMG, significance lost after adjustment | 42 rows, t-test; **none significant after BH-FDR** | ~ consistent with the paper's "lost significance" note |

---

## Results that reproduce the paper

Several of the paper's core biological conclusions are reproduced despite the methodological differences above.

### 1. Differential-expression direction for the paper's key genes
Results in `results/tables/deg_list.csv` confirm the direction of the paper's Fig. 2 for every overlapping gene (logFC = HPV+ minus HPV−):

| Gene | logFC | Direction | Paper's finding |
|---|---|---|---|
| `ZFR2` | +6.74 | overexpressed in HPV+ | ✓ overexpressed in HPV+ |
| `SYCP2` | +5.23 | overexpressed in HPV+ | ✓ |
| `SOX30` | +4.26 | overexpressed in HPV+ | ✓ |
| `MEI1` | +3.34 | overexpressed in HPV+ | ✓ |
| `CCNA1` | −3.91 | underexpressed in HPV+ | ✓ |
| `SPRR2G` | −4.75 | underexpressed in HPV+ | ✓ |
| `GJB6` | −3.14 | underexpressed in HPV+ | ✓ |
| `MMP3` | −3.30 | underexpressed in HPV+ | ✓ (Yellow network) |
| `FLRT3` | −3.40 | underexpressed in HPV+ | ✓ |
| `PITX2` | −2.45 | underexpressed in HPV+ | ✓ |

### 2. Negative expression–methylation correlation
All **10** doubly-selected DEG+DMG genes of this replication show **negative Pearson correlations** between expression and promoter methylation in *both* HPV groups (`results/tables/deg_dmg_correlation.csv`), matching the paper's central claim that promoter methylation represses expression:

```
 gene     rho_hpv_pos   rho_hpv_neg
 SYCP2      -0.697        -0.181
 DDX25      -0.749        -0.143
 MEI1       -0.660        -0.386
 MSX2       -0.658        -0.493
 SLC35F3    -0.418        -0.362
 CCNA1      -0.532        -0.562
 POU4F1     -0.345        -0.783
 WNK2       -0.166        -0.530
 INA        -0.781        -0.653
 NEFL       -0.716        -0.482
```

Three of the paper's doubly-selected genes are shared with this replication: **`SYCP2`, `MEI1`, `CCNA1`**.

### 3. Module membership of key hub genes
From `results/tables/wgcna_module_assignments.csv` (mapped Ensembl → symbol):

| Module | Replication genes | Paper's assignment |
|---|---|---|
| **blue** | `SYCP2`, `HSF4`, `MYO15B`, `EYA2`, `UGT8`, `KRT19` | Blue ✓ (SYCP2/HSF4/MYO15B/EYA2/UGT8 in Blue; KRT19 in Yellow in the paper) |
| **yellow** | `FLRT3`, `GJB6`, `HOXC13`, `MMP3`, `KRT14`, `SPRR2G`, `DMRTA2`, `YBX2`, `MYO3A` | Yellow ✓ for FLRT3/GJB6/HOXC13/MMP3/KRT14; DMRTA2/YBX2/MYO3A were Blue in the paper |
| **grey** | `CCNA1`, `CDKN2A`, `PITX2`, `CTSE`, `PAX1` | Grey ✓ for all five |
| **turquoise** | `MEI1`, `SOX30`, `ZFR2`, `COL4A6` | Turquoise ✓ (paper Table 3 lists MEI1/SOX30/ZFR2 as Turquoise; not refined in the paper) |

So for the modules that exist in **both** studies the hub-gene assignments agree well — in particular `SYCP2` (the paper's headline hub) lands in the blue module, connected with its DMG neighbours `HSF4` and `MYO15B` (see below).

### 4. SYCP2 hub in the Blue HPV+ network
The blue HPV+ graph (`results/networks/blue_hpv_pos.graphml`) contains `SYCP2`, `HSF4`, `MYO15B`, `EYA2`, `UGT8`; `SYCP2` has degree 38. In the blue HPV− graph none of these appear except `EYA2`. This mirrors the paper's Fig. 4A: `SYCP2` is densely connected in HPV+ and effectively absent from the HPV− network, and `HSF4`/`MYO15B` (DMG in both studies) are co-expressed with it.

### 5. Yellow module keratinocyte/epidermis biology
The yellow network is enriched in the same processes the paper reports in Table 2 — and with extreme significance:

- `epidermis development`, `skin development`, `keratinocyte differentiation`, `epidermal cell differentiation`, `keratinization`, `cornification` (GO BP)
- `Cornified envelope formation` (KEGG)
- in **both** the HPV+ and HPV− yellow networks

### 6. HPV+ networks are denser than HPV−
In every module that was refined, the HPV+ graph has more nodes and edges than the HPV− graph (see comparison table below) — reproducing the paper's most general conclusion that HPV+ tumours show more densely connected co-expression networks.

---

## Deliberate differences and their rationale

These decisions were made consciously during the implementation and are recorded here for transparency.

1. **Global DEG/DMG test instead of per-stage tests, and `|logFC| ≥ 2` for DEG instead of `≥ 4`.**
   The paper tests HPV+ vs HPV− *within* stage I–IV and unions the results. With just 2–5 HPV+ samples in stages I–III, per-stage tests are extremely under-powered. A single well-powered global comparison (36 vs 241) was chosen to obtain stable, reproducible differential lists. Because a global test finds smooth, large changes more easily than small-n stage tests, the stricter `|logFC| ≥ 4` cut-off would have almost emptied the DEG list (only ~40 genes pass), so the threshold was relaxed to `≥ 2`. Consequences: DEG = 628 (vs 223), DMG = 39 (vs 359) — the DMG drop reflects both the global (vs per-stage-union) design and the missing BMIQ/probe-QC (see gaps).

2. **STAR-Counts + TMM log-CPM instead of RNAseqV2 RSEM.**
   The current GDC no longer ships the 2015-era `RNAseqV2` RSEM level-3 files used by the paper; the default quantification for TCGA-HNSC is `STAR - Counts`. TMM log-CPM is the standard, best-practice normalisation for these raw counts and is a reasonable substitute for the paper's `log(1+p)` transform of RSEM counts.

3. **Expression rows collapsed Ensembl → symbol.**
   DEG/DMG, WGCNA and the networks all operate on gene symbols for consistency with the paper and with the methylation/mutation annotations. Ties are collapsed by mean; WGCNA keeps the highest-MAD duplicate — a pragmatic mapping rather than an official ID translation.

4. **TF list from dorothea instead of TFcheckpoint.**
   TFcheckpoint's downloads are no longer reliably accessible; dorothea's A/B-confidence human regulons provide a curated, current TF gene list. This is only used for node annotation (`is_TF`) and does not change network topology.

5. **STRING v12.0 confidence ≥ 700 instead of STRING v9.1 ≥ 0.7.**
   STRING v9.1 is long retired; v12.0 with the same 700 (0.7) score threshold implements the paper's "high confidence" criterion on current data.

6. **MSigDB Hallmark (H) only.**
   The annotation package `msigdbr` provides the H (Hallmark) set out of the box. The paper's C6 cancer-signature sets (`RICKMAN_HEAD_AND_NECK_CANCER_A`, `PYEON_HPV_POSITIVE_TUMORS_UP` in its Table 2) come from MSigDB v6.0 and require a separate licensed download, so they were omitted. This is why those two signature enrichments are absent here.

7. **ClusterProfiler used per-network on refined network nodes.**
   The paper reports enrichment of "genes within modules"; the replication runs enrichment on the *refined network node sets* (which are exactly the module genes passing the DEG/edge filter). In the paper's Table 2 the listed genes are all network nodes, so this is a faithful interpretation.

8. **FDR ≤ 0.05 everywhere for enrichment.**
   The paper's methods say `p ≤ 0.001` FDR-adjusted, yet Table 2 lists p-values up to 0.042; the replication consistently uses FDR ≤ 0.05.

---

## Unintentional gaps and blockers

These are places where the replication does **not** match the paper, and that should be kept in mind when interpreting the results. They are documented factually rather than "recommended to fix".

1. **No BMIQ normalisation and no full probe black-list.**
   `o4_preprocess_methylation.R` removes sex-chromosome and >5%-missing probes and kNN-imputes, but the paper additionally removes cross-reactive, non-specific, SNP-bearing and zero-coordinate probes and applies **BMIQ** normalisation (ChAMP pipeline). Skipping these inflates the methylation gene count (19,160 vs 14,861) and is the most likely driver of the much smaller DMG list (39 vs 359).

2. **No per-stage analysis.**
   The paper explicitly analyses disease progression (stage-by-stage DEG/DMG, Supplementary Table 1 and Supplementary Fig. 1). This replication uses a single global comparison, so the stage-specific results in the paper — e.g. "only six DEG selected in all stages; no DMG common across stages" — are not reproduced.

3. **Only 3 clinical traits in the module–trait analysis.**
   The paper correlates module eigengenes with 7 traits (HPV status, staging, age, gender, alcohol, smoked, anatomical site); the replication could only align `hpv_status`, `age_at_index` and `ajcc_pathologic_stage` from the downloaded clinical record.

4. **Grey module just misses significance.**
   The replication's grey module has `cor = −0.242` vs HPV status — below the paper's `|cor| ≥ 0.25` threshold — so it was not refined into networks, and **no grey networks exist** under `results/networks/`. The paper's Grey module findings (PITX2/CCNA1 sub-network, TP53/CDKN2A node in HPV−) therefore have no counterpart here. Instead the (statistically significant) turquoise module was refined — a module the paper detected but did not study.

5. **Enrichment uses only Hallmark, not C6.**
   Consequently the `RICKMAN_HEAD_AND_NECK_CANCER_A` and `PYEON_HPV_POSITIVE_TUMORS_UP` rows of the paper's Table 2 cannot be reproduced.

6. **Mutation calling simplified.**
   The paper calls significantly mutated genes genome-wide with `maftools` + FDR correction and integrates only those via STRING. Here, Fisher's exact test is applied only to genes already present in a network, with raw `p ≤ 0.05` and no multiple-testing correction.

7. **Caveats in the validation step (o10).**
   The script auto-detects the HPV-status column by searching for `"hpv"` in the **column names** of each GEO series matrix:
   - **GSE6791 (Pyeon et al.):** the series matrix contains **no HPV-status field at all** (case/site/gender/age/stage only), so the function bails out and no `gse6791_validation.csv` is produced. The paper's Table 3 congruence (38 genes) therefore could not be reproduced from GEO metadata alone.
   - **GSE38266 (Lechner et al.):** HPV status is present but embedded in a *value string* (`tissue: HPV- HNSCC tumor` / `tissue: HPV+ HNSCC tumor`), not in a column name — so it is not detected and no rows are produced. (Its 42 samples did download successfully.)
   - **GSE95036 (Esposti et al.):** a `hpv status:` column exists, so 42 genes were tested with per-gene t-tests on promoter-beta means. Raw p-values range down to 3.9e-3 (`MYO15B`), but **after Benjamini–Hochberg correction none of the 42 genes is significant** — consistent with the paper's own statement that significance was lost in this dataset once p-values were adjusted. The replication applies no `|logFC| ≥ 1.5` filter and tests raw beta (not M values), so it is a weaker check than the paper's.
   - **Dataset-size discrepancy:** the paper describes GSE38266 as 11 samples (6 HPV+) fresh-frozen and GSE95036 as 42 samples (21 HPV+) FFPE. In the *current* GEO records the sizes are swapped — GSE95036 is the 11-sample fresh-frozen set and GSE38266 holds 42 samples. This is documented here so readers are not misled by the paper's numbers.

8. **Figures not generated.**
   Only the module–trait heatmap (paper Fig. 3) is produced (`results/figures/fig03_module_trait_heatmap/module_trait_heatmap.png`). The per-gene correlation scatter plots (Fig. 2), the network diagrams (Fig. 4) and the validation plots (Fig. 5) were not rendered; the corresponding figure folders exist but are empty.

9. **Filename/docs inconsistencies.**
   `o7_wgcma_modules.R` is a typo for the WGCNA step; script headers reference missing `docs/0*.md` files (this README replaces them).

---

## Comparison of achieved results vs the paper

### Module–trait association (HPV status)

Replication (from `results/tables/wgcna_module_trait_correlation.csv`):

| Module | cor (HPV status) | p | Significant (`|cor|≥0.25, p≤0.001`)? |
|---|---|---|---|
| MEyellow | −0.442 | 1.2e-14 | ✓ |
| MEturquoise | +0.316 | 7.8e-8 | ✓ |
| MEblue | +0.267 | 6.8e-6 | ✓ |
| MEred | +0.243 | 4.2e-5 | borderline (|cor| < 0.25) |
| MEgrey | −0.242 | 4.8e-5 | borderline (|cor| < 0.25) |
| MEbrown | −0.139 | 2.0e-2 | ✗ |
| MEgreen | −0.110 | 6.9e-2 | ✗ |
| MEblack | −0.015 | 8.1e-1 | ✗ |

Paper: significant modules = **Blue, Yellow, Grey**. Replication = **Yellow, Turquoise, Blue** — with grey just missing the 0.25 cut-off.

Module sizes (replication): black 23, red 27, green 521, brown 1,061, blue 1,462, grey 1,948, yellow 955, turquoise 2,003 (8 modules vs the paper's 7).

### Network connection metrics — paper Table 1 vs replication

| Module | Status | Paper nodes / edges | Replication nodes / edges |
|---|---|---|---|
| Blue | HPV+ | 539 / 2,633 | **312 / 1,649** |
| Blue | HPV− | 112 / 114 | **132 / 418** |
| Yellow | HPV+ | 145 / 640 | **569 / 10,583** |
| Yellow | HPV− | 36 / 34 | **320 / 4,325** |
| Grey | HPV+ | 127 / 111 | **not built** (below significance) |
| Grey | HPV− | 8 / 7 | **not built** |
| Turquoise | HPV+ | not studied | 1,095 / 149,226 |
| Turquoise | HPV− | not studied | 767 / 80,664 |

Shared qualitative result: within each module the HPV+ network is denser (more nodes/edges) than HPV−. Absolute numbers are not reproduced (this replication keeps *all* pairwise correlations passing `rho ≥ 0.65, p ≤ 0.01` among DEG-connected genes, which with 36 HPV+ samples and large yellow/turquoise modules yields very dense graphs).

### Enrichment highlights — paper Table 2 vs replication

| Paper finding | Category | Replication status |
|---|---|---|
| Blue HPV+ — cell-fate specification, glycolysis/gluconeogenesis, glycerophospholipid metabolism, stem-cell pluripotency, RICKMAN/PYEON C6 signatures | GO/KEGG/MSigDB C6 | GO "cell fate **commitment**" found, but in the **turquoise** HPV+ network; no glycolysis/gluconeogenesis; C6 unavailable; **xenobiotic metabolism** instead dominates blue HPV+ |
| Yellow HPV+ — epidermis/skin development, keratinocyte differentiation, negative regulation of epithelial proliferation, EMT | GO/MSigDB H | **Reproduced**: epidermis development, skin development, keratinocyte differentiation, keratinization/cornification, Cornified-envelope formation in the yellow networks (both statuses); EMT hallmark found in **blue HPV−** instead of yellow HPV+ |
| Grey C1/C2 — dendritic spine morphogenesis, lysosome | GO/KEGG | **Not reproducible** — grey module not refined |
| Grey HPV− — cellular ketone metabolism, aging, xenobiotic metabolism | GO/MSigDB H | **Not reproducible** — grey module not refined (xenobiotic signals appear in the blue HPV+ network instead) |

### Validation summary

| Dataset | Paper outcome | Replication outcome |
|---|---|---|
| GSE6791 (microarray) | 38 concordant DE genes | **Blocked** — no HPV field in the GEO series matrix |
| GSE38266 (450K) | DMG congruence (6/12 double-selected genes) | **Blocked** — HPV status not in a detectable column name |
| GSE95036 (450K) | DMG congruence, lost on adjustment | Raw t-test p-values produced; **none survive FDR** — consistent with the paper |

---

## Reproduction instructions

### Environment

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

### Run order

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

### Expected inputs, outputs, and runtime caveats

- **Disk:** the raw TCGA download (o1) is on the order of several GB across the four omic layers, plus ~380 MB of GEO soft/series files for o10. The two `Wed_Sep__9_…*.tar.gz` archives in the repo root are leftover earlier-stage downloads and are not required by the pipeline.
- **Runtime:** the slowest steps are o4 (kNN imputation + probe→gene collapse over ~400k probes) and o7 (WGCNA on 8,000 genes), each typically taking several minutes on a multicore machine (`enableWGCNAThreads()` is called in o7). o10's GEO downloads dominate when run cold.
- **Reproducibility caveats:**
  - Cohort count (277 = 241 HPV− + 36 HPV+) reflects the state of the GDC and cBioPortal `hnsc_tcga_pub` records at download time; if those sources change, counts may shift slightly.
  - WGCNA auto-selects the soft-threshold power (fallback 6); a different power changes the module partition. Expect the *module colours* to differ across runs/environments even though the biology is stable.
  - Network edge counts depend on STRING availability (PPI edges) and on which DEG list precedes the network step.
  - `o10` uses a brittle auto-detection of HPV columns in GEO pheno data; inspect `pData(eset)` for GSE6791/GSE38266 if you need those validations (see [Unintentional gaps](#unintentional-gaps-and-blockers)).

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

---

## References

- Costa R.L., Boroni M. & Soares M.A. (2018). *Distinct co-expression networks using multi-omic data reveal novel interventional targets in HPV-positive and negative head-and-neck squamous cell cancer.* Scientific Reports 8:15254. [DOI:10.1038/s41598-018-33498-5](https://doi.org/10.1038/s41598-018-33498-5)
- Lawrence M.S. et al. (2015). *Comprehensive genomic characterization of head and neck squamous cell carcinomas.* Nature 517:576–582 (source TCGA-HNSC cohort).
- Langfelder P. & Horvath S. (2008/2005). WGCNA / Dynamic Tree Cut.
- Ritchie M.E. et al. (2015). limma. *Nucleic Acids Research* 43:e47.
- Yu G. et al. (2012/2016). clusterProfiler. *OMICS* 16:284–287.
- Original data availability from the paper: interactive networks and Neo4j HPV+ graph database at https://github.com/quelopes/HNSCC-network.