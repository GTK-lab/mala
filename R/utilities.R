globalVariables(".")  # Oh, the joy of programming in R.

#' Resolve file or URL
#'
#' @param path File path or URL.
#'
#' @details If \code{path} is a local file path, assert that it
#'   exists, then just return it. If \code{path} is a URL, see if it
#'   has been downloaded before by checking if the file
#'   \code{basename(path)} exists. If it exists, just return
#'   \code{basename(path)}. Otherwise, download the file at
#'   \code{path} to \code{basename(path)} and return
#'   \code{basename(path)}.
#'
#' @importFrom assertthat assert_that is.string
#' @importFrom urltools scheme
#' @importFrom utils download.file
download_if_url <- function(path) {
  assert_that(is.string(path))

  if (is.na(scheme(path)) || scheme(path) == "file") {
    assert_that(file.exists(path))
    return(path)
  }

  # Otherwise, first check if it had already been downloaded.
  if (file.exists(basename(path))) {
    message(sprintf("Skipping download: file %s already exists", basename(path)))
    return(basename(path))
  }

  # If not already downloaded, download it.
  message(sprintf("Downloading to %s...", basename(path)))
  download.file(path, basename(path), "auto")

  basename(path)
}

#' Decompress a gzipped tarball
#'
#' @param path File path to gzipped tarball.
#'
#' @details The contents of the tarball are extracted into a directory which has
#'   directory name equals to the basename of the \code{path}. For example,
#'   \code{a.tar.gz} containing \code{b.txt} will generate directory \code{./a/}
#'   and file \code{./a/b.txt}.
#'
#' @importFrom assertthat assert_that is.string is.dir
#' @importFrom stringr str_extract
#' @importFrom utils untar
extract_if_targz <- function(path) {
  assert_that(is.string(path))

  is_targz <- grepl("\\.tar\\.gz$", path)
  dirn <- ifelse(is_targz, str_extract(path, ".*(?=\\.tar\\.gz$)"), path)

  if (is_targz & !dir.exists(dirn)) {
    message(sprintf("Extracting %s to %s...", path, dirn))
    untar(path, exdir=dirn)
  } else if (is_targz & dir.exists(dirn)) {
    message(sprintf("Skip extracting %s, because %s exists", path, dirn))
  } else { }  # !is_targz, so nothing to do.

  assert_that(is.dir(dirn))
  dirn
}

#' Convert sleuth results data.frame to GRanges
#'
#' @param sr \code{data.frame} from \code{sleuth_results},
#' @param mean_obs_col_name Column name to find the quantity representing
#'   mean log observations. For sleuth \code{gene_mode=TRUE}, this is
#'   \code{"mean_obs"}, for sleuth \code{pval_aggregate=TRUE}, this is
#'   \code{"sum_mean_obs_count"}.
#'
#' @return GRanges with qval, mean_obs, and target_id (usually a unique
#'   TSS-associated ID) as metadata.
#'
#' @details Rows with \code{NA} qvals are dropped. Whatever values are in the
#'   column with the name \code{mean_obs_col_name} will be saved in the mcol
#'   column named "mean_obs".
#'
#' @importFrom assertthat assert_that has_name not_empty
#' @importFrom tibble as_tibble
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @importFrom tidyr separate
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
sr2gr <- function(sr, mean_obs_col_name="mean_obs") {

  assert_that(has_name(sr, "target_id"))
  assert_that(has_name(sr, "qval"))
  assert_that(has_name(sr, mean_obs_col_name))
  assert_that(not_empty(sr))

  sr2 <- sr %>%
    as_tibble() %>%
    select(.data$target_id, .data$qval, .env$mean_obs_col_name) %>%
    filter(!is.na(.data$qval)) %>%
    # loci information is encoded in tss_id, which is in the target_id column
    separate(
      .data$target_id, into=c("seqnames", "start", "end", "strand"), sep=",",
      remove=FALSE, convert=TRUE)

  # The piping breaks here, because it is difficult to continue with a column
  # that can have different names.
  gr <- GRanges(
    seqnames=sr2$seqnames,
    ranges=IRanges(sr2$start, sr2$end),
    strand=sr2$strand,
    tss_id=sr2$target_id,
    qval=sr2$qval,
    mean_obs=pull(sr2[, mean_obs_col_name]))

  message(sprintf("Removed %d rows with NA qvals", nrow(sr) - length(gr)))

  gr
}

#' Convert tx2tss tibble to GRanges
#'
#' @param tx2tss tibble from get_tx2tss
#'
#' @return GRanges object with metadata column target_id.
#' @importFrom assertthat assert_that not_empty
#' @importFrom rlang .data
#' @importFrom tidyr separate
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
tx2tss2gr <- function(tx2tss) {
  assert_that(not_empty(tx2tss))
  tx2tss %>%
    separate(
      .data$tss_id, into=c("seqnames", "start", "end", "strand"), sep=",",
      convert=TRUE) %>%
    { GRanges(
      seqnames=.$seqnames, ranges=IRanges(start=.$start, end=.$end),
      strand=.$strand, target_id=.$target_id) }
}

#' Filter a GRanges object based on its number of overlaps with another
#'
#' @param gr_to_filter a GRanges object to filter
#' @param gr_to_overlap a GRanges object to overlap \code{to_filter} with
#' @param filter_fn A unary function which takes a numeric vector of length
#'   \code{length(to_filter)}, where the i'th element of that vector is the
#'   number of overlaps of the i'th row of \code{to_filter} with
#'   \code{gr_to_overlap}. This function should then return a logical vector
#'   (usually of length also equals \code{length(to_filter)}) with which to
#'   subset (index) \code{to_filter} with.
#'
#' @return GRanges subset of \code{gr_to_filter}.
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom dplyr group_by summarise n pull
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tibble as_tibble
filter_by_num_overlaps <- function(gr_to_filter, gr_to_overlap, filter_fn) {
  overlaps <- suppressWarnings(findOverlaps(gr_to_filter, gr_to_overlap)) %>%
    as_tibble() %>%
    # Filtering is implicit --- if a row in gr_to_filter does not overlap with
    # gr_to_overlap, then its row index will not be in .data$queryHits.
    group_by(.data$queryHits) %>%
    summarise(num_overlaps=n())

  if (nrow(overlaps) == 0) {
    return(GRanges())  # no overlaps to report
  } else {  # otherwise...
    overlaps_idxs_to_keep <- overlaps %>%
      pull(.data$num_overlaps) %>%
      filter_fn()
    # NB overlaps_idxs_to_keep is not the same as gr_to_filter_idxs_to_keep.
    # e.g.        overlaps$queryHits == c(11, 29, 37),
    #             overlaps$num_overlaps == c( 5,  1,  3),
    #          overlaps_idxs_to_keep == c( T,  F,  F),
    #      gr_to_filter_idxs_to_keep == c(11).
    # Also, gr_to_filter_idxs_to_keep is guaranteed unique by virtue of group_by
    gr_to_filter_idxs_to_keep <- overlaps$queryHits[overlaps_idxs_to_keep]
    gr_to_filter[gr_to_filter_idxs_to_keep, ]
  }
}
