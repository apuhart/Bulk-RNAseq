# Bulk RNA-seq Pipeline

Differential expression and pathway enrichment pipeline for bulk RNA-seq data (Salmon → DESeq2 → fgsea), developed for iPSC-derived neuronal models of **Friedreich's Ataxia (FA)**.

---

## Table of contents

- [Overview](#overview)
- [Repository structure](#repository-structure)
- [Dependencies](#dependencies)
- [Configuration](#configuration)
- [Usage](#usage)
- [Outputs](#outputs)
- [Pipeline details](#pipeline-details)

---

## Overview

This pipeline takes **Salmon quantification outputs** (`quant.sf`) as input and runs:

1. **Transcript-to-gene aggregation** via `tximport` + Ensembl BioMart
2. **Differential expression** via `DESeq2` — one model per contrast, with sex as covariate
3. **Visualization** — PCA (2D/3D), volcano plots, heatmaps (top DEGs + genes of interest)
4. **Gene set enrichment** via `fgsea` with MSigDB collections
5. **Cross-dataset comparison** — scatter plots and z-scored heatmaps across two datasets (e.g. iPSC vs iNeurons)

Designed for multi-dataset, multi-contrast experiments with automatic contrast detection (control vs FA, control vs FA_edited, FA_edited vs FA).

---

## Repository structure

```
Bulk-RNAseq/
│
├── mainPatchSave.R                  # Pipeline orchestrator — run this
│
├── config/
│   ├── config.R                     # All user-defined parameters (paths, samples, thresholds)
│   ├── packages.R                   # Package installation and loading
│   └── biomart.R                    # BioMart tx2gene helper
│
├── R/
│   ├── deseq2_helpers.R             # DESeq2 wrapper (tximport → DESeq2 → results)
│   ├── fgsea_helpers.R              # fgsea wrappers + MSigDB collection loading
│   ├── plots.R                      # PCA, volcano, heatmap plotting functions
│   └── cross_dataset_comparison.R   # Cross-dataset scatter + z-scored heatmaps
│
└── results/                         # Generated outputs (gitignored)
    └── <dataset>/
        └── <contrast>/
```

---

## Dependencies

### R packages

| Package | Role |
|---|---|
| `DESeq2` | Differential expression |
| `tximport` | Salmon → gene counts |
| `biomaRt` | tx2gene mapping (Ensembl) |
| `fgsea` | Gene set enrichment |
| `msigdbr` | MSigDB gene set collections |
| `pheatmap` | Heatmaps |
| `plotly` | Interactive 3D PCA |
| `ggplot2` | Plots |
| `ggrepel` | Volcano labels |
| `htmlwidgets` | Save interactive HTML plots |
| `dplyr` | Data wrangling |

Install all dependencies via `config/packages.R`:

```r
source("config/packages.R")
```

### External requirement

**Salmon** quantification must be run upstream. Each sample needs a `quant.sf` file at:
```
<base_path>/<sample_name>/quant.sf
```

---

## Configuration

All parameters are set in **`config/config.R`**. Key variables:

```r
# Output directory
BASE_OUT <- "results/"

# Ensembl mirror for BioMart
ENSEMBL_MIRROR <- "https://www.ensembl.org"

# DESeq2 thresholds
MIN_COUNTS     <- 10       # minimum count filter
LFC_THRESHOLD  <- 1        # |log2FC| cutoff for significance
PADJ_THRESHOLD <- 0.05     # adjusted p-value cutoff

# Heatmap
HEATMAP_TOP_N  <- 50       # number of top DEGs in heatmap

# fgsea
FGSEA_COLLECTIONS <- c("H", "C2", "C5")   # MSigDB collections
FGSEA_NPERM       <- 1000
FGSEA_TOP_N       <- 20
FGSEA_PADJ_CUT    <- 0.05
NEURONAL_KEYWORDS <- c("NEURON", "AXON", "SYNAP", "DOPAMIN", "SEROTON")

# Genes of interest (for focused volcano + heatmap)
GENES_OF_INTEREST <- c("FXN", "ISCU", "NFS1", "ACO2", "SOD2", "GPX4")

# Dataset definitions
DATASETS <- list(
  iPSC = list(
    path       = "/path/to/iPSC/salmon/",
    samples    = c("iPSC_Ctrl_1", "iPSC_Ctrl_2", "iPSC_FA_1", "iPSC_FA_2"),
    conditions = c("control", "control", "FA", "FA"),
    sex        = c("Male", "Female", "Male", "Female")
  ),
  iNeurons = list(
    path       = "/path/to/iNeurons/salmon/",
    samples    = c("iN_Ctrl_1", "iN_Ctrl_2", "iN_FA_1", "iN_FA_2"),
    conditions = c("control", "control", "FA", "FA"),
    sex        = c("Male", "Female", "Male", "Female")
  )
)
```

**Contrasts are detected automatically** based on condition labels:

| Conditions present | Contrast run |
|---|---|
| `control` + `FA` | `control_vs_FA` |
| `control` + `FA_edited` | `control_vs_FA_edited` |
| `FA_edited` + `FA` | `FA_edited_vs_FA` |

---

## Usage

```r
# From R or RStudio — adjust paths in config/config.R first
source("mainPatchSave.R")
```

Or from the terminal:

```bash
Rscript mainPatchSave.R
```

---

## Outputs

For each dataset and each contrast, the following files are generated under `results/<dataset>/<contrast>/`:

| File | Description |
|---|---|
| `<contrast>_PCA_2D.pdf` | 2D PCA colored by condition and sex |
| `<contrast>_PCA_3D.html` | Interactive 3D PCA (plotly) |
| `<contrast>_volcano.pdf` | Volcano plot — all genes |
| `<contrast>_volcano_GOI.pdf` | Volcano plot — genes of interest highlighted |
| `<contrast>_heatmap_topDEGs.pdf` | Heatmap of top N DEGs (VSD z-scored) |
| `<contrast>_heatmap_GOI.pdf` | Heatmap restricted to genes of interest |
| `<contrast>_fgsea_*.pdf` | fgsea enrichment plots per MSigDB collection |
| `<contrast>_fgsea_*.csv` | fgsea full results table |

Cross-dataset outputs are saved under `results/SCATTER/`:

| File | Description |
|---|---|
| `scatter_<A>_vs_<B>.pdf` | log2FC scatter across two datasets |
| `heatmap_cross_<A>_<B>.pdf` | z-scored expression heatmap — shared DEGs |

---

## Pipeline details

### 1. tx2gene mapping

`config/biomart.R` connects to Ensembl via `biomaRt` to retrieve transcript-to-gene mappings (Ensembl transcript ID → Ensembl gene ID → HGNC symbol). This mapping is used by `tximport` to aggregate transcript-level Salmon quantifications to gene level.

### 2. DESeq2 model

For each contrast, the design formula is:

```
~ sex + condition   (if sex is not confounded with condition)
~ condition         (otherwise)
```

Results are shrunk with `lfcShrink` (apeglm) and filtered by `MIN_COUNTS` before testing.

### 3. fgsea

Gene ranks are computed as `sign(log2FC) × -log10(pvalue)`. Three MSigDB collections are tested by default: **H** (Hallmarks), **C2** (curated gene sets), **C5** (GO terms). Neuronal pathways are highlighted separately using keyword filtering.

### 4. Cross-dataset comparison

`R/cross_dataset_comparison.R` aligns DEG results from two datasets (e.g. iPSC vs iNeurons) on shared genes and produces:
- A **log2FC scatter** to identify concordant vs discordant regulation
- A **z-scored expression heatmap** of shared significant DEGs using the VSD matrices

---

## Biological context

This pipeline was developed for the study of **Friedreich's Ataxia (FA)**, a neurodegenerative disease caused by frataxin (*FXN*) deficiency. Datasets include:

- **iPSC** — induced pluripotent stem cells (Ctrl / FA / FA_edited)
- **iNeurons** — iPSC-derived neurons at different differentiation timepoints (Div40 / Div80)

Key pathways of interest: iron-sulfur cluster biogenesis (FXN/ISCU/NFS1), mitochondrial dysfunction, and oxidative stress (SOD2, GPX4, NFE2L2).

---

## Author

**A-P. Uhart**
