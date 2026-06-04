# ============================================================
# R/functions/deseq2_helpers.R — DESeq2 wrappers
# ============================================================

#' Run DESeq2 for a pairwise contrast
#'
#' @param files_pair   Named character vector of quant.sf paths
#' @param samples_pair Character vector of sample names
#' @param conds_pair   Character vector of conditions (same order as samples)
#' @param donor_pair   Character vector of sex/donor covariates
#' @param tx2gene      tx2gene data.frame
#' @param condA        Reference condition (denominator)
#' @param condB        Numerator condition
#' @param min_counts   Minimum rowSums(counts) to keep a gene
#' @return list: dds, res, resLFC, res_df, vsd
run_deseq2_contrast <- function(
    files_pair, samples_pair, conds_pair, donor_pair,
    tx2gene, condA, condB,
    min_counts = 3
) {

  # --- tximport ---
  txi <- tximport(
    files_pair,
    type       = "salmon",
    tx2gene    = tx2gene[, c("TXNAME", "GENEID")],
    ignoreTxVersion = TRUE
  )

  # --- colData ---
  colData <- data.frame(
    sample    = samples_pair,
    sex       = factor(donor_pair),
    condition = factor(conds_pair, levels = c(condA, condB)),
    row.names = samples_pair
  )

  # --- DESeq2 ---
  dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > min_counts, ]
  dds <- DESeq(dds)
  dds <- estimateSizeFactors(dds)

  # --- Results ---
  res <- results(dds, contrast = c("condition", condB, condA))

  resLFC <- tryCatch(
    lfcShrink(dds, contrast = c("condition", condB, condA), type = "apeglm"),
    error = function(e) {
      message("apeglm lfcShrink failed, returning unshrunken results: ", e$message)
      res
    }
  )

  # --- Annotate ---
  res_df <- res %>%
    as.data.frame() %>%
    tibble::rownames_to_column("GENEID") %>%
    dplyr::mutate(GENEID = sub("\\..*", "", GENEID)) %>%
    annotate_results(tx2gene) %>%
    dplyr::arrange(padj) %>%
    dplyr::distinct(GENEID, .keep_all = TRUE) %>%
    dplyr::select(GENEID, GENENAME, gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

  # --- VST ---
  vsd <- tryCatch(
    vst(dds),
    error = function(e) { message("vst failed: ", e$message); NULL }
  )

  list(dds = dds, res = res, resLFC = resLFC, res_df = res_df, vsd = vsd, colData = colData)
}


#' Build ranked gene list for fgsea from a DESeq2 results data.frame
#'
#' @param res_df  data.frame with columns gene and log2FoldChange
#' @return Named numeric vector sorted decreasingly
make_ranks <- function(res_df) {
  ranks <- res_df$log2FoldChange
  names(ranks) <- res_df$gene
  ranks <- ranks[!is.na(names(ranks)) & !is.na(ranks)]
  sort(ranks, decreasing = TRUE)
}


#' Check which genes of interest are present/absent in DESeq2 results
#'
#' @param res_df            data.frame with gene column
#' @param genes_of_interest character vector of gene symbols
#' @return list: present, missing
check_genes_of_interest <- function(res_df, genes_of_interest) {
  present <- intersect(genes_of_interest, res_df$gene)
  missing <- setdiff(genes_of_interest, res_df$gene)
  cat("Genes of interest — total:", length(genes_of_interest),
      "| present:", length(present),
      "| missing:", length(missing), "\n")
  list(present = present, missing = missing)
}


#' Extract sex-chromosome marker expression for QC
#'
#' @param dds  DESeqDataSet after estimateSizeFactors
#' @return matrix (or empty matrix if none found)
check_sex_expression <- function(dds) {
  sex_genes <- c(
    "ENSG00000229807",  # XIST
    "ENSG00000067048",  # RPS4Y1
    "ENSG00000012817",  # KDM5D
    "ENSG00000129824",  # EIF1AY
    "ENSG00000198692",  # UTY
    "ENSG00000000003"   # DDX3Y
  )
  norm_counts <- counts(dds, normalized = TRUE)
  rownames(norm_counts) <- sub("\\..*", "", rownames(norm_counts))
  norm_counts[rownames(norm_counts) %in% sex_genes, , drop = FALSE]
}
