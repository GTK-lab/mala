#' Aggregate GRanges-represented sleuth results using TF to loci mappings
#'
#' @param gr GRanges representation of sleuth results obtained using
#'   \code{sr2gr}.
#' @param tf2loci List of GRanges representing transcription factor to loci
#'   bindings, obtained using \code{get_tf2loci}.
#' @param weight_fn Unary function to produce weights for Lancaster p-value
#'   aggregation. This function should take a length-N numeric vector of
#'   mean_obs, where N is the number of loci for that transcription factor.
#'
#' @return A tibble of results.
#'
#' @references Yi, L., Pimentel, H., Bray, N.L. et al. Gene-level
#'   differential analysis at transcript-level resolution. Genome Biol
#'   19, 53 (2018). https://doi.org/10.1186/s13059-018-1419-z
#'
#' @importFrom magrittr %>%
#' @importFrom assertthat assert_that not_empty
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom aggregation lancaster
#' @importFrom dplyr bind_rows mutate select
aggregate_gr <- function(
    gr, tf2loci, weight_fn=function(obs) log(obs+1)/sum(log(obs+1))) {
  # TODO custom weight function

  assert_that(not_empty(gr))
  assert_that(not_empty(tf2loci))

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
        gr[., ]
      pval <- lancaster(gr_sub$qval, weight_fn(gr_sub$mean_obs))
      tibble(pval=.data$pval)
    }) %>%
    bind_rows() %>%
    mutate(TF=names(tf2loci)) %>%
    select(.data$TF, .data$pval)  # reorder columns
  message()  # just for the newline

  results %>% arrange(.data$pval)
}
