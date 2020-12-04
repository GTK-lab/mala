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
#'
#' @return GRanges with qval, mean_obs, and target_id (usually a unique
#'   TSS-associated ID) as metadata.
#'
#' @details Rows with \code{NA} qvals are dropped.
#'
#' @importFrom assertthat assert_that has_name not_empty
#' @importFrom tibble as_tibble
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @importFrom tidyr separate
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
sr2gr <- function(sr) {

  assert_that(has_name(sr, "target_id"))
  assert_that(has_name(sr, "qval"))
  assert_that(has_name(sr, "mean_obs"))
  assert_that(not_empty(sr))

  gr <- sr %>%
    as_tibble() %>%
    select(.data$target_id, .data$qval, .data$mean_obs) %>%
    filter(!is.na(.data$qval)) %>%
    # loci information is encoded in tss_id, which is in the target_id column
    separate(
      .data$target_id, into=c("seqnames", "start", "end", "strand"), sep=",",
      remove=FALSE, convert=TRUE) %>%
    { GRanges(
        seqnames=.$seqnames,
        ranges=IRanges(.$start, .$end),
        strand=.$strand,
        tss_id=.$target_id,
        qval=.$qval,
        mean_obs=.$mean_obs) }

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
