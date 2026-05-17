# Insomnia Immune–Stress Bioinformatics — Supplementary Code

[![DOI](https://zenodo.org/badge/DOI/REPLACE_WITH_DOI_AFTER_RELEASE.svg)](https://doi.org/REPLACE_WITH_DOI_AFTER_RELEASE)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0007--4575--1368-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0007-4575-1368)

This repository contains all analytical code used to produce the results of:

> **Oxidative Stress Shows the Largest Immune-Pathway Effect Size in Chronic Insomnia: Integrative Identification of SPTLC1, PTGES3, and HLA-G as Candidate Hub Genes**
> Xin, Z. *PeerJ* (under review), 2026.
> Manuscript DOI: TBA upon acceptance.

The analysis uses **GSE208668** (Illumina HumanHT-12 V4 microarray, PBMC, n = 42:
17 chronic insomnia patients vs. 25 controls), publicly available at the
[NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208668).

---

## What's in this repository

| Path | What it is |
|---|---|
| `Insomnia_Supplementary_Code_ALL.R` | **Single-file merged script** — all 12 analysis phases in execution order. The recommended entry point. |
| `00_setup_and_data.R` ... `11_Supplementary_Tables.R` | The same code split into 12 standalone scripts (alternative for users who prefer to run section-by-section). |
| `RUN_ALL.R` | Orchestration script that sources all 12 in order. |
| `README.md` | This file. |
| `LICENSE` | MIT license. |
| `CITATION.cff` | How to cite this code (machine-readable). |
| `sessionInfo.txt` | R + package versions used for the published results (regenerated at the end of `Insomnia_Supplementary_Code_ALL.R`). |

---

## Pipeline overview

```
GSE208668 (GEO)
      │
      ▼
[ Section 1 ]  GEO download + preprocessing  (14,539 genes × 42 samples)
      │
      ├──▶ [ Section 2 ]  Phase 1 — GSVA of 5 stress / immune gene sets
      │
      ├──▶ [ Section 3 ]  Phase 2 — WGCNA modules  (turquoise: r = −0.44, n = 2,178)
      │       │
      │       ▼
      │   [ Section 4 ]  Phase 3 — ML (LASSO + RF + Boruta + LOOCV)  →  3 core genes
      │
      ├──▶ [ Section 5 ]  Phase 3.5 — Directional consistency across 3 sleep cohorts
      │
      ├──▶ [ Section 6 ]  Phase 0 — eQTLGen cis-eQTL lookup
      │
      ├──▶ [ Section 7 ]  Phase 4 — CIBERSORTx immune deconvolution  (external web tool)
      │
      ├──▶ [ Section 8 ]  Phase 5 — SMR Mendelian randomization      (external CLI tool)
      │
      ├──▶ [ Section 9 ]  Phase 6 — DGIdb drug–target lookup         (external web tool)
      │
      ├──▶ [ Section 10 ] GO / KEGG enrichment of the turquoise module
      ├──▶ [ Section 11 ] Bootstrap + LOOCV diagnostics of the 3-gene classifier
      └──▶ [ Section 12 ] Export all supplementary tables (Table S2 etc.)
```

---

## Requirements

- **R ≥ 4.3**
- CRAN packages: `tidyverse`, `data.table`, `glmnet`, `randomForest`, `Boruta`, `pROC`, `VennDiagram`, `WGCNA`, `ggplot2`, `openxlsx`, `BiocManager`
- Bioconductor packages: `GEOquery`, `limma`, `GSVA`, `clusterProfiler`, `org.Hs.eg.db`

Install everything at once:

```r
install.packages(c(
  "tidyverse", "data.table", "glmnet", "randomForest", "Boruta",
  "pROC", "VennDiagram", "WGCNA", "ggplot2", "openxlsx", "BiocManager"
))
BiocManager::install(c(
  "GEOquery", "limma", "GSVA", "clusterProfiler", "org.Hs.eg.db"
))
```

Exact versions used in the published analysis are recorded in `sessionInfo.txt`.

---

## How to reproduce

1. Open `Insomnia_Supplementary_Code_ALL.R` in RStudio.
2. Edit the `PROJECT_DIR` line near the top to point to a working directory on your machine.
3. Source the file:
   ```r
   source("Insomnia_Supplementary_Code_ALL.R", encoding = "UTF-8")
   ```
4. Three sections (7, 8, 9) involve external tools that cannot be invoked from R:
   - **Section 7** — upload the matrix to [CIBERSORTx](https://cibersortx.stanford.edu) (free academic account required)
   - **Section 8** — run the [SMR command-line tool](https://yanglab.westlake.edu.cn/software/smr/) (Linux/Mac/Windows binaries available)
   - **Section 9** — paste 9 gene symbols into [DGIdb](https://dgidb.org)
   Each section contains the exact instructions for the external step and the parsing code for the returned files.

---

## External data sources

| Resource | URL | Access date |
|---|---|---|
| GSE208668 (microarray, PBMC) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208668 | April 2026 |
| GSE39445, GSE98566, GSE56931 (cross-cohort) | https://www.ncbi.nlm.nih.gov/geo | April 2026 |
| eQTLGen Phase I cis-eQTL | https://www.eqtlgen.org/ | April 2026 |
| Jansen 2019 Insomnia GWAS | https://ctg.cncr.nl/software/summary_statistics | April 2026 |
| CIBERSORTx LM22 reference | https://cibersortx.stanford.edu | April 2026 |
| DGIdb v5.x | https://dgidb.org | April 2026 |
| 1000 Genomes EUR LD reference (for SMR) | https://yanglab.westlake.edu.cn/software/smr/ | April 2026 |

> **Note on third-party databases:** DGIdb and eQTLGen are versioned databases that may update over time. The drug-gene interactions and cis-eQTL associations reported in this manuscript reflect query results retrieved on the dates listed above.

---

## Important caveats (as stated in the manuscript)

- The 3-gene classifier (SPTLC1, PTGES3, HLA-G) achieved **AUC = 1.000** under leave-one-out cross-validation, but this is a **discovery-stage estimate** in a small cohort and **requires external validation** in independent samples.
- Cross-cohort analyses (Section 5) test **directional consistency**, not formal replication.
- SMR results (Section 8) are reported as evidence **"consistent with"** a causal effect, not as proof of causality.
- CIBERSORTx-derived immune cell fractions (Section 7) are **relative estimates**, not absolute cell counts.

---

## How to cite this code

If you use this code, please cite both the manuscript and the archived code snapshot:

```
Xin, Z. (2026). Insomnia Immune–Stress Bioinformatics — Supplementary Code (v1.0.0)
[Software]. Zenodo. https://doi.org/REPLACE_WITH_DOI_AFTER_RELEASE
```

Machine-readable citation info is in `CITATION.cff`.

---

## License

MIT License — see `LICENSE` file.

This means: you are free to use, modify, and redistribute this code,
provided that the copyright notice and license text are preserved.

---

## Contact

For questions about the code or methods, open an issue on this repository,
or contact:

- **Zhipeng Xin** — School of Basic Medical Sciences, Hubei University of Chinese Medicine
- Email: jin85956@stmail.hbucm.edu.cn
- ORCID: [0009-0007-4575-1368](https://orcid.org/0009-0007-4575-1368)
