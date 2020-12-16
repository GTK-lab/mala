#' Aggregate GRanges-represented sleuth results using TF to loci mappings
#'
#' @param gr GRanges representation of sleuth results obtained using
#'   \code{sr2gr}.
#' @param tf2loci List of GRanges representing transcription factor to loci
#'   bindings, obtained using \code{get_tf2loci}.
#' @param intermediates If this function should return a list containing
#'   intermediate data structures used in the computation.
#' @param weight_fn Unary function to produce weights for Lancaster p-value
#'   aggregation. This function will receive \code{mcols(gr)}, and should return
#'   a numeric vector of the same length as \code{length(gr)}. See details.
#' @param filter_mapping_fn Unary function to filter out TSS to TF mappings
#'   which seem poorly supported by \code{tf2loci}. This function will receive
#'   a numeric vector which represents the number of rows in tf2loci which
#'   support each TSS to TF mapping, for each TSS (that is, each numeric vector
#'   contains TSS to TF mappings for one TSS only). Then, the function should
#'   return a logical vector representing which TSS to TF mappings should be
#'   used. To keep all TSS to TF mappings without filtering, use
#'   \code{filter_mapping_fn=function(.) TRUE}.
#'
#' @details This function will add a column "tfs_overlapped" to \code{mcols(gr)}
#'   which is the number of transcription factors which each TSS in \code{gr}
#'   maps to, if that column does not yet exists. So, you can use that column
#'   in \code{weight_fn}.
#'
#' @return If \code{intermediates=FALSE}, a tibble of results, one row per TF.
#'   If \code{intermediates=TRUE}, a list with named elements:
#'   \describe{
#'     \item{gr}{The GRanges of the argument \code{gr}, possibly with an added
#'     metadata columns "tfs_overlapped" --- see details.}
#'     \item{grl}{A list of GRanges. Each GRanges is derived from overlapping
#'       the argument \code{gr} with a transcription factor from the argument
#'       \code{tf2loci}.}
#'     \item{tfs}{The tibble of results which would have been returned by itself
#'       if \code{intermediates=FALSE}}.
#'   }
#'
#' @references Yi, L., Pimentel, H., Bray, N.L. et al. Gene-level
#'   differential analysis at transcript-level resolution. Genome Biol
#'   19, 53 (2018). https://doi.org/10.1186/s13059-018-1419-z
#'
#' @importFrom GenomicRanges findOverlaps mcols mcols<- reduce
#' @importFrom aggregation lancaster
#' @importFrom assertthat assert_that not_empty has_name
#' @importFrom dplyr bind_rows mutate select group_by summarise n arrange pull
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats p.adjust median
#' @importFrom tibble as_tibble
aggregate_gr <- function(
    gr, tf2loci, intermediates=FALSE,
    weight_fn=function(mcols)
      (1/mcols$tfs_overlapped) / sum(1/mcols$tfs_overlapped),
    filter_mapping_fn=function(n_rows)
      n_rows > median(n_rows)) {

  assert_that(not_empty(gr))
  assert_that(not_empty(tf2loci))

  ## Preprocess to generate "tfs_overlapped" metadata column
  if (has_name(mcols(gr), "tfs_overlapped")) {
    message("Column tfs_overlapped already exists in gr, not overwriting.")
  } else {
    message("Creating column tfs_overlapped.")
    # For each TSS, how many transcription factors do they map?
    # (Below, GenomicRanges::reduces merges loci to avoid double-counting. I
    # think these overlapping tf2loci regions come about because UniBind bulk
    # mappings draw from many different experiments.)
    hits <- suppressWarnings(findOverlaps(gr, reduce(tf2loci))) %>%
      as_tibble() %>%
      group_by(.data$queryHits) %>%
      summarise(tfs_overlapped=n())
    tfs_overlapped <- rep(0, length(gr))
    tfs_overlapped[hits$queryHits] <- hits$tfs_overlapped
    mcols(gr)$tfs_overlapped <- tfs_overlapped
  }

  ## For each transcription factor, find its overlaps with gr, filtering by the
  ## number of overlaps (keep only TSS in gr which overlap at least some number
  ## of times with UniBind BED records).
  message(sprintf(
    "Calculating. Each dot represents a transcription factor (total: %d).",
    length(tf2loci)))
  # /FOR EACH/ transcription factor
  grl <- lapply(tf2loci, function(tfloci) {
      message(".", appendLF=FALSE)
      gr_sub <- filter_by_num_overlaps(gr, tfloci, filter_mapping_fn)
      gr_sub
    })

  ## For each TSS-TF filtered overlap, aggregate their p-values.
  tfs <- grl %>%
    lapply(function(gr_sub) {
        # NAs are ignored by lancaster
        pval <- lancaster(gr_sub$qval, weight_fn(mcols(gr_sub)))
        tibble(pval=pval)
      }) %>%
    bind_rows() %>%
    mutate(TF=names(tf2loci)) %>%
    select(.data$TF, .data$pval) %>%  # reorder columns
    mutate(qval=p.adjust(.data$pval, "BH")) %>%
    arrange(.data$pval)
  message()  # just for the newline

  if (intermediates) list(gr=gr, grl=grl,tfs=tfs) else tfs
}
