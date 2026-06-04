# ============================================================
# R/functions/packages.R — Package installation & loading
# ============================================================

load_packages <- function() {

  bioc_pkgs <- c(
    "tximport", "DESeq2", "biomaRt", "apeglm",
    "fgsea", "limma", "edgeR"
  )

  cran_pkgs <- c(
    "dplyr", "ggplot2", "ggrepel", "pheatmap",
    "matrixStats", "msigdbr", "WGCNA",
    "plotly", "htmlwidgets", "tibble"
  )

  # Install BiocManager if needed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  # Install missing Bioconductor packages
  missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE)
  }

  # Install missing CRAN packages
  missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_cran) > 0) {
    message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran)
  }

  # Load all
  all_pkgs <- c(bioc_pkgs, cran_pkgs)
  invisible(lapply(all_pkgs, library, character.only = TRUE))
  message("All packages loaded.")
}
