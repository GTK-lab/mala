#' Aggregate GRanges-represented sleuth results using TF to loci mappings
#'
#' @param gr GRanges representation of sleuth results obtained using
#'   \code{sr2gr}.
#' @param tf2loci List of GRanges representing transcription factor to loci
#'   bindings, obtained using \code{get_tf2loci}.
#' @param weight_fn Unary function to produce weights for Lancaster p-value
#'   aggregation. This function will receive \code{mcols(gr)}, and should return
#'   a numeric vector of the same length as \code{length(gr)}. See details.
#'
#' @details This function will add a column "tfs_overlapped" to \code{mcols(gr)}
#'   which is the number of transcription factors which each TSS in \code{gr}
#'   maps to, if that column does not yet exists. So, you can use that column
#'   in \code{weight_fn}.
#'
#' @return A tibble of results.
#'
#' @references Yi, L., Pimentel, H., Bray, N.L. et al. Gene-level
#'   differential analysis at transcript-level resolution. Genome Biol
#'   19, 53 (2018). https://doi.org/10.1186/s13059-018-1419-z
#'
#' @importFrom GenomicRanges findOverlaps mcols mcols<- reduce
#' @importFrom S4Vectors queryHits
#' @importFrom aggregation lancaster
#' @importFrom assertthat assert_that not_empty has_name
#' @importFrom dplyr bind_rows mutate select group_by summarise n arrange
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust
aggregate_gr <- function(
    gr, tf2loci,
    weight_fn=function(mcols)
      (1/mcols$tfs_overlapped) / sum(1/mcols$tfs_overlapped)) {

  assert_that(not_empty(gr))
  assert_that(not_empty(tf2loci))

  if (has_name(mcols(gr), "tfs_overlapped")) {
    message("Column tfs_overlapped already exists in gr, not overwriting.")
  } else {
    message("Creating column tfs_overlapped.")
    # GenomicRanges::reduces merges loci to avoid double-counting
    hits <- suppressWarnings(findOverlaps(gr, reduce(tf2loci))) %>%
      as_tibble() %>%
      group_by(queryHits) %>%
      summarise(tfs_overlapped=n())
    tfs_overlapped <- rep(0, length(gr))
    tfs_overlapped[hits$queryHits] <- hits$tfs_overlapped
    mcols(gr)$tfs_overlapped <- tfs_overlapped
  }

  message(sprintf(
    "Calculating. Each dot represents a transcription factor (total: %d).",
    length(tf2loci)))
  # /FOR EACH/ transcription factor to loci binding
  results <- lapply(tf2loci, function(tfloci) {
      message(".", appendLF=FALSE)
      # Different BED files include different 'esoteric' sequences (chromosomes)
      # and this generates a warning from `c` when combining GRanges.
      gr_sub <- suppressWarnings(findOverlaps(gr, tfloci)) %>%
        queryHits() %>%
        unique() %>%
        gr[., ]
      pval <- lancaster(gr_sub$qval, weight_fn(mcols(gr_sub)))
      tibble(pval=pval)
    }) %>%
    bind_rows() %>%
    mutate(TF=names(tf2loci)) %>%
    select(.data$TF, .data$pval) %>%  # reorder columns
    mutate(qval=p.adjust(.data$pval, "BH"))
  message()  # just for the newline

  results %>% arrange(.data$pval)
}
