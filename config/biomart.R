# ============================================================
# R/functions/biomart.R — Transcript-to-gene mapping via BioMart
# ============================================================

#' Fetch tx2gene table from Ensembl
#'
#' @param mirror  Ensembl mirror to use (default: "useast")
#' @return data.frame with columns TXNAME, GENEID, GENENAME
get_tx2gene <- function(mirror = "useast") {

  message("Connecting to Ensembl BioMart (mirror: ", mirror, ") ...")

  ensembl <- useEnsembl(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    mirror  = mirror
  )

  tx2gene <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
    mart       = ensembl
  )

  colnames(tx2gene) <- c("TXNAME", "GENEID", "GENENAME")

  # Strip version suffixes (e.g. ENSG00000001234.5 -> ENSG00000001234)
  tx2gene$TXNAME  <- sub("\\..*", "", tx2gene$TXNAME)
  tx2gene$GENEID  <- sub("\\..*", "", tx2gene$GENEID)

  message("tx2gene fetched: ", nrow(tx2gene), " transcripts, ",
          length(unique(tx2gene$GENEID)), " genes.")
  return(tx2gene)
}


#' Map Ensembl gene IDs to gene symbols in a results data.frame
#'
#' @param res_df   data.frame with a GENEID column
#' @param tx2gene  tx2gene table (output of get_tx2gene)
#' @return res_df with GENENAME and gene columns added
annotate_results <- function(res_df, tx2gene) {

  gene_map <- tx2gene %>%
    dplyr::select(GENEID, GENENAME) %>%
    dplyr::distinct()

  res_df %>%
    dplyr::left_join(gene_map, by = "GENEID") %>%
    dplyr::mutate(
      gene = dplyr::coalesce(GENENAME, GENEID),
      gene = ifelse(is.na(gene) | gene == "", GENEID, gene)
    )
}
