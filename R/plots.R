# ============================================================
# R/functions/plots.R — All plotting helpers
# ============================================================

# ----------------------------
# Volcano plot
# ----------------------------

#' Standard volcano plot (all significant genes labelled)
plot_volcano <- function(
    res_df, dataset_name, contrast_name, direction_label,
    lfc_thr = 0.3, padj_thr = 0.05
) {
  df <- res_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(sig = padj < padj_thr & abs(log2FoldChange) > lfc_thr)

  ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = sig), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    coord_cartesian(
      xlim = c(min(df$log2FoldChange) - 3, max(df$log2FoldChange) + 3)
    ) +
    geom_text_repel(
      data = subset(df, padj < 0.1 & abs(log2FoldChange) > lfc_thr),
      aes(label = gene), size = 5, max.overlaps = 50
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x     = "log2 Fold Change",
      y     = "-log10 adjusted p-value",
      color = "Significant",
      title = paste("Volcano -", dataset_name, contrast_name, "\n", direction_label)
    )
}


#' Focused volcano plot — labels only genes of interest
plot_volcano_focus <- function(
    res_df, dataset_name, contrast_name, direction_label,
    genes_of_interest,
    lfc_thr = 0.3, padj_thr = 0.05
) {
  df <- res_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(
      sig         = padj < padj_thr & abs(log2FoldChange) > lfc_thr,
      is_interest = gene %in% genes_of_interest
    )

  ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = sig), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    geom_text_repel(
      data     = subset(df, is_interest),
      aes(label = gene),
      fontface = "bold", size = 4, max.overlaps = 50
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x     = "log2 Fold Change",
      y     = "-log10 adjusted p-value",
      color = "Significant",
      title = paste("Volcano (genes of interest) -", dataset_name, contrast_name, "\n", direction_label)
    )
}


# ----------------------------
# PCA plots
# ----------------------------

#' 2D PCA with condition colour + sex shape
plot_pca_2d <- function(vsd, dataset_name, contrast_name) {
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pca_data$sex <- colData(vsd)$sex[match(pca_data$name, rownames(colData(vsd)))]
  percentVar   <- round(100 * attr(pca_data, "percentVar"))

  ggplot(pca_data, aes(PC1, PC2)) +
    geom_point(aes(color = condition, shape = sex), size = 4) +
    geom_text_repel(aes(label = name), size = 4) +
    scale_color_manual(values = c(
      "control"   = "steelblue",
      "FA"        = "firebrick",
      "FA_edited" = "darkorange"
    )) +
    scale_shape_manual(values = c("M" = 17, "F" = 16)) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) +
    ggtitle(paste("PCA -", dataset_name, contrast_name)) +
    theme_minimal(base_size = 14)
}


#' 3D interactive PCA (plotly) — returns a plotly object
plot_pca_3d <- function(vsd, dataset_name, contrast_name) {
  mat      <- t(assay(vsd))
  pca      <- prcomp(mat, scale. = TRUE)
  pca_data <- as.data.frame(pca$x[, 1:3])
  pca_data$condition <- colData(vsd)$condition
  pca_data$sex       <- colData(vsd)$sex
  pca_data$name      <- rownames(pca_data)
  pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))

  color_map <- c("control" = "steelblue", "FA" = "firebrick", "FA_edited" = "darkorange")
  shape_map <- c("M" = "triangle-up", "F" = "circle")

  p <- plot_ly(
    data         = pca_data,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color        = ~condition, colors  = color_map,
    symbol       = ~sex,       symbols = shape_map,
    text         = ~name,
    type         = "scatter3d",
    mode         = "markers+text",
    textposition = "top center",
    marker       = list(size = 6)
  ) %>% plotly::layout(
    scene = list(
      xaxis = list(title = paste0("PC1: ", pct[1], "%")),
      yaxis = list(title = paste0("PC2: ", pct[2], "%")),
      zaxis = list(title = paste0("PC3: ", pct[3], "%"))
    ),
    title = paste("3D PCA -", dataset_name, contrast_name)
  )
  p
}


# ----------------------------
# Heatmaps
# ----------------------------

#' Heatmap of top N DEGs (by padj)
plot_heatmap_top_degs <- function(
    vsd_mat, res_df, colData, sample_order,
    dataset_name, contrast_name,
    top_n = 100, padj_thr = 0.05
) {
  top_genes <- res_df %>%
    dplyr::filter(padj < padj_thr) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(GENEID)

    # Keep only genes present in vsd_mat
  top_genes <- top_genes[top_genes %in% rownames(vsd_mat)]

  if (length(top_genes) == 0) {
    message("No significant DEGs to plot heatmap for ", contrast_name)
    return(invisible(NULL))
  }

  pheatmap(
    vsd_mat[top_genes, sample_order, drop = FALSE],
    scale          = "row",
    annotation_col = as.data.frame(colData[sample_order, ]),
    main           = paste("Top", top_n, "DEGs -", dataset_name, contrast_name),
    fontsize_row   = 16,
    fontsize_col   = 14,
    fontsize       = 16,
    cluster_cols   = FALSE
  )
}




plot_heatmap_goi <- function(
    vsd_mat, colData, sample_order,
    genes_of_interest,
    gene_map,
    dataset_name, contrast_name
) {
  
  # convert SYMBOL → ENSG si nécessaire
  if (!all(genes_of_interest %in% rownames(vsd_mat))) {
    genes_of_interest <- gene_map$GENEID[
      gene_map$GENENAME %in% genes_of_interest
    ]
  }
  
  goi_present <- genes_of_interest[
    genes_of_interest %in% rownames(vsd_mat)
  ]
  
  if (length(goi_present) == 0) {
    message("No genes of interest found in vsd for ", contrast_name)
    return(invisible(NULL))
  }
  
  pheatmap(
    vsd_mat[goi_present, sample_order, drop = FALSE],
    scale = "row",
    annotation_col = as.data.frame(colData[sample_order, ]),
    main = paste("Genes of Interest -", dataset_name, contrast_name),
    fontsize_row = 16,
    fontsize_col = 14,
    cluster_cols = FALSE,
    cluster_rows = FALSE
  )
}
#' 
#' #' Heatmap of genes of interest
#' plot_heatmap_goi <- function(
#'     vsd_mat, colData, sample_order,
#'     genes_of_interest,
#'     dataset_name, contrast_name
#' ) {
#'   
#'   browser()
#'   goi_present <- genes_of_interest[genes_of_interest %in% rownames(vsd_mat)]
#' 
#'   if (length(goi_present) == 0) {
#'     message("No genes of interest found in vsd for ", contrast_name)
#'     return(invisible(NULL))
#'   }
#' 
#'   pheatmap(
#'     vsd_mat[goi_present, sample_order, drop = FALSE],
#'     scale          = "row",
#'     annotation_col = as.data.frame(colData[sample_order, ]),
#'     main           = paste("Genes of Interest -", dataset_name, contrast_name),
#'     fontsize_row   = 16,
#'     fontsize_col   = 14,
#'     fontsize       = 16,
#'     cluster_cols   = FALSE,
#'     cluster_rows   = FALSE
#'   )
#' }


# ----------------------------
# Scatter (cross-contrast)
# ----------------------------

#' Build and save scatter plot comparing two contrasts
#'
#' @param comp1,comp2   Names of contrasts in all_results
#' @param all_results   Named list of DESeq2 result data.frames
#' @param out_dir       Root output directory for this scatter comparison
#' @param dataset_name  Dataset label for plot titles
run_scatter <- function(comp1, comp2, all_results, out_dir, dataset_name) {

  if (!all(c(comp1, comp2) %in% names(all_results))) return(invisible(NULL))

  # ---- directory layout ----
  comp_out      <- file.path(out_dir, paste0(comp1, "_VS_", comp2))
  full_dir      <- file.path(comp_out, "FULL")
  sig_dir       <- file.path(comp_out, "SIGNIFICANT")
  quad_full_dir <- file.path(full_dir,  "QUAD")
  quad_sig_dir  <- file.path(sig_dir,   "QUAD")
  for (d in c(quad_full_dir, quad_sig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  diag_slope <- ifelse(comp2 == "FA_vs_FA_edited", -1, 1)

  fc1 <- paste0("log2FC_", comp1)
  fc2 <- paste0("log2FC_", comp2)
  p1  <- paste0("padj_",   comp1)
  p2  <- paste0("padj_",   comp2)

  # ---- merge ----
  merged <- dplyr::full_join(
    all_results[[comp1]] %>%
      dplyr::select(GENEID, gene_1 = gene, log2FoldChange, padj) %>%
      dplyr::rename(!!fc1 := log2FoldChange, !!p1 := padj),
    all_results[[comp2]] %>%
      dplyr::select(GENEID, gene_2 = gene, log2FoldChange, padj) %>%
      dplyr::rename(!!fc2 := log2FoldChange, !!p2 := padj),
    by = "GENEID"
  ) %>%
    dplyr::mutate(
      gene        = dplyr::coalesce(gene_1, gene_2),
      gene        = ifelse(is.na(gene) | gene == "", GENEID, gene),
      dist_diag   = abs(.data[[fc2]] - diag_slope * .data[[fc1]]) / sqrt(1 + diag_slope^2),
      significant = (.data[[p1]] < 0.05 | .data[[p2]] < 0.05)
    )

  write.csv(
    merged,
    file.path(comp_out, paste0(dataset_name, "_", comp1, "_VS_", comp2, "_scatter_table.csv")),
    row.names = FALSE
  )

  # ---- inner helpers ----
  make_plot <- function(df, filename, save_dir, label_type = "none") {
    p <- ggplot(df, aes(x = .data[[fc1]], y = .data[[fc2]])) +
      geom_point(
        aes(fill = dist_diag, color = significant),
        shape = 21, size = 2.3, alpha = 0.85, stroke = 0.4
      ) +
      scale_fill_viridis_c(name = "Distance to diagonal") +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "transparent"), guide = "none") +
      geom_abline(slope = diag_slope, intercept = 0, linetype = "dashed", color = "grey60") +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_vline(xintercept = 0, linetype = "dotted") +
      theme_minimal(base_size = 14) +
      labs(
        x     = paste("Log2FC", comp1),
        y     = paste("Log2FC", comp2),
        title = paste(dataset_name, "-", comp1, "vs", comp2)
      )

    if (label_type != "none") {
      label_col <- ifelse(label_type == "ENSG", "GENEID", "gene")
      label_df  <- df %>%
        dplyr::mutate(score = abs(.data[[fc1]]) + abs(.data[[fc2]])) %>%
        dplyr::arrange(dplyr::desc(score)) %>%
        dplyr::slice_head(n = 50)
      p <- p + geom_text_repel(data = label_df, aes(label = .data[[label_col]]),
                               size = 3, max.overlaps = 50)
    }
    ggplot2::ggsave(file.path(save_dir, filename), p, width = 10, height = 8, bg = "white")
  }

  plot_set <- function(df, suffix, save_dir) {
    base <- paste0(dataset_name, "_", comp1, "_VS_", comp2, "_", suffix)
    make_plot(df, paste0(base, ".png"),       save_dir, "none")
    make_plot(df, paste0(base, "_ENSG.png"),  save_dir, "ENSG")
    make_plot(df, paste0(base, "_genes.png"), save_dir, "gene")
  }

  get_quadrants <- function(df) {
    if (diag_slope == 1) {
      list(
        hd = df %>% dplyr::filter(.data[[fc1]] > 0 & .data[[fc2]] > 0),
        bg = df %>% dplyr::filter(.data[[fc1]] < 0 & .data[[fc2]] < 0)
      )
    } else {
      list(
        hd = df %>% dplyr::filter(.data[[fc1]] < 0 & .data[[fc2]] > 0),
        bg = df %>% dplyr::filter(.data[[fc1]] > 0 & .data[[fc2]] < 0)
      )
    }
  }

  # ---- run plots ----
  plot_set(merged, "FULL", full_dir)
  quads <- get_quadrants(merged)
  plot_set(quads$hd, "quad_HD", quad_full_dir)
  plot_set(quads$bg, "quad_BG", quad_full_dir)

  sig_df    <- merged %>% dplyr::filter(.data[[p1]] < 0.05 | .data[[p2]] < 0.05)
  plot_set(sig_df, "SIGNIFICANT", sig_dir)
  quads_sig <- get_quadrants(sig_df)
  plot_set(quads_sig$hd, "quad_HD", quad_sig_dir)
  plot_set(quads_sig$bg, "quad_BG", quad_sig_dir)

  invisible(merged)
}
