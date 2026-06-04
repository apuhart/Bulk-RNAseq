options(repos = c(CRAN = "https://cloud.r-project.org"))

# ============================================================
# main.R — RNA-seq pipeline orchestrator
# ============================================================

# ----------------------------
# 0. Sources
# ----------------------------
source("/home/apuhart/Downloads/rna_seq_pip/config/packages.R")
source("/home/apuhart/Downloads/rna_seq_pip/config/config.R")
source("/home/apuhart/Downloads/rna_seq_pip/config/biomart.R")
source("/home/apuhart/Downloads/rna_seq_pip/R/deseq2_helpers.R")
source("/home/apuhart/Downloads/rna_seq_pip/R/fgsea_helpers.R")
source("/home/apuhart/Downloads/rna_seq_pip/R/plots.R")
source("/home/apuhart/Downloads/rna_seq_pip/R/cross_dataset_comparison.R")

# ----------------------------
# 1. Global setup
# ----------------------------
library(biomaRt)
library(tximport)
library(DESeq2)
library(plotly)
library(pheatmap)
library(msigdbr)
library(fgsea)

tx2gene <- get_tx2gene(ENSEMBL_MIRROR)

# gene mapping (ENSG <-> SYMBOL)
gene_map <- tx2gene %>%
  dplyr::select(GENEID, GENENAME) %>%
  dplyr::distinct()

# convert GOI symbols -> ENSG
goi_ensg <- gene_map$GENEID[
  gene_map$GENENAME %in% GENES_OF_INTEREST
]

dir.create(BASE_OUT, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2. Helper: run one dataset
# ----------------------------
run_dataset <- function(dataset_name, dataset) {
  
  message("\n==============================")
  message("DATASET: ", dataset_name)
  message("==============================\n")
  
  out_dir <- file.path(BASE_OUT, dataset_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  samples    <- dataset$samples
  conditions <- dataset$conditions
  sex        <- dataset$sex
  base_path  <- dataset$path
  
  # ----------------------------
  # 3. Define contrasts
  # ----------------------------
  cond_levels <- unique(conditions)
  contrasts   <- list()
  
  if (all(c("control", "FA") %in% cond_levels))
    contrasts[["control_vs_FA"]] <- c("control", "FA")
  
  if (all(c("control", "FA_edited") %in% cond_levels))
    contrasts[["control_vs_FA_edited"]] <- c("control", "FA_edited")
  
  if (all(c("FA_edited", "FA") %in% cond_levels))
    contrasts[["FA_edited_vs_FA"]] <- c("FA_edited", "FA")
  
  all_results <- list()
  vsd_store   <- list()   # <-- collects vsd per contrast
  
  # ----------------------------
  # 4. Run DESeq2 per contrast
  # ----------------------------
  for (cn in names(contrasts)) {
    
    contrast_dir <- file.path(out_dir, cn)
    dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)
    
    message("\n--- Contrast: ", cn, " ---")
    
    condA <- contrasts[[cn]][1]
    condB <- contrasts[[cn]][2]
    idx   <- which(conditions %in% c(condA, condB))
    
    files <- file.path(base_path, samples[idx], "quant.sf")
    names(files) <- samples[idx]
    
    res <- run_deseq2_contrast(
      files_pair  = files,
      samples_pair = samples[idx],
      conds_pair   = conditions[idx],
      donor_pair   = sex[idx],
      tx2gene      = tx2gene,
      condA        = condA,
      condB        = condB,
      min_counts   = MIN_COUNTS
    )
    
    all_results[[cn]] <- res$res_df
    vsd_store[[cn]]   <- res$vsd   # save vsd object for cross-dataset heatmap
    
    # ----------------------------
    # 5. PCA
    # ----------------------------
    if (!is.null(res$vsd)) {
      pdf(file.path(contrast_dir, paste0(cn, "_PCA_2D.pdf")))
      print(plot_pca_2d(res$vsd, dataset_name, cn))
      dev.off()
      
      p3d <- plot_pca_3d(res$vsd, dataset_name, cn)
      htmlwidgets::saveWidget(
        p3d,
        file.path(contrast_dir, paste0(cn, "_PCA_3D.html")),
        selfcontained = TRUE
      )
    }
    
    # ----------------------------
    # 6. Volcano plots
    # ----------------------------
    pdf(file.path(contrast_dir, paste0(cn, "_volcano.pdf")), width = 15, height = 12)
    print(plot_volcano(
      res$res_df, dataset_name, cn,
      direction_label = paste(condB, "vs", condA),
      lfc_thr  = LFC_THRESHOLD,
      padj_thr = PADJ_THRESHOLD
    ))
    dev.off()
    
    pdf(file.path(contrast_dir, paste0(cn, "_volcano_GOI.pdf")), width = 15, height = 12)
    print(plot_volcano_focus(
      res$res_df, dataset_name, cn,
      direction_label   = paste(condB, "vs", condA),
      genes_of_interest = GENES_OF_INTEREST,
      lfc_thr  = LFC_THRESHOLD,
      padj_thr = PADJ_THRESHOLD
    ))
    dev.off()
    
    # ----------------------------
    # 7. Heatmaps
    # ----------------------------
    sample_order <- rownames(res$colData)
    
    pdf(file.path(contrast_dir, paste0(cn, "_heatmap_topDEGs.pdf")), width = 10, height = 20)
    print(plot_heatmap_top_degs(
      vsd_mat      = assay(res$vsd),
      res_df       = res$res_df,
      colData      = res$colData,
      sample_order = sample_order,
      dataset_name = dataset_name,
      contrast_name = cn,
      top_n        = HEATMAP_TOP_N,
      padj_thr     = PADJ_THRESHOLD
    ))
    dev.off()
    
    pdf(file.path(contrast_dir, paste0(cn, "_heatmap_GOI.pdf")), width = 10, height = 20)
    print(plot_heatmap_goi(
      vsd_mat           = assay(res$vsd),
      colData           = res$colData,
      sample_order      = sample_order,
      gene_map          = gene_map,
      genes_of_interest = goi_ensg,
      dataset_name      = dataset_name,
      contrast_name     = cn
    ))
    dev.off()
    
    # ----------------------------
    # 8. fgsea
    # ----------------------------
    ranks <- make_ranks(res$res_df)
    
    run_all_fgsea(
      ranks             = ranks,
      out_dir           = contrast_dir,
      label             = cn,
      collections       = FGSEA_COLLECTIONS,
      neuronal_keywords = NEURONAL_KEYWORDS,
      direction_label   = paste(condB, "vs", condA),
      nperm             = FGSEA_NPERM,
      top_n             = FGSEA_TOP_N,
      padj_cut          = FGSEA_PADJ_CUT
    )
  }
  
  # ----------------------------
  # 9. Cross-contrast scatter
  # ----------------------------
  cnames <- names(all_results)
  
  if (length(cnames) >= 2) {
    for (i in 1:(length(cnames) - 1)) {
      for (j in (i + 1):length(cnames)) {
        run_scatter(
          comp1       = cnames[i],
          comp2       = cnames[j],
          all_results = all_results,
          out_dir     = file.path(out_dir, "SCATTER"),
          dataset_name = dataset_name
        )
      }
    }
  }
  
  message("\nFinished dataset: ", dataset_name)
  
  # Return both DEG results and vsd matrices
  return(invisible(list(deg = all_results, vsd = vsd_store)))
}

# ----------------------------
# 10. Run all datasets
# ----------------------------
all_dataset_results <- list()
all_vsd_results     <- list()

for (ds_name in names(DATASETS)) {
  out <- run_dataset(ds_name, DATASETS[[ds_name]])
  all_dataset_results[[ds_name]] <- out$deg
  all_vsd_results[[ds_name]]     <- lapply(out$vsd, function(v) assay(v))
}

# ----------------------------
# 11. Cross-dataset comparison
# ----------------------------
DATASET_A <- "iPSC"       # must match keys in DATASETS (config.R)
DATASET_B <- "iNeurons"   # must match keys in DATASETS (config.R)

run_cross_dataset_comparison(
  all_dataset_results = all_dataset_results,
  base_out            = BASE_OUT,
  name_a              = DATASET_A,
  name_b              = DATASET_B,
  padj_thr            = PADJ_THRESHOLD,
  lfc_thr             = LFC_THRESHOLD,
  genes_highlight     = GENES_OF_INTEREST,
  vsd_list            = all_vsd_results    # enables z-scored expression heatmaps
)

message("\n✔ ALL ANALYSES COMPLETED (including cross-dataset comparison)")