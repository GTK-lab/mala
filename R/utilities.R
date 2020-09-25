globalVariables(".")  # Oh, the joy of programming in R.

#' Resolve file or URL
#' @keywords internal
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
#' @importFrom RCurl url.exists
#' @importFrom utils download.file
download_if_url <- function(path) {
  assert_that(is.string(path))

  # If path is a local file path...
  if (!url.exists(path)) {
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
#' @keywords internal
#'
#' @param path File path to gzipped tarball.
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
#' @return GRanges with qva, mean_obs, and target_id (usually a unique
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
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
sr2gr <- function(sr) {

  assert_that(has_name(sr, "target_id"))
  assert_that(has_name(sr, "qval"))
  assert_that(has_name(sr, "mean_obs"))
  assert_that(not_empty(sr))

  sr %>%
    as_tibble() %>%
    select(.data$target_id, .data$qval, .data$mean_obs) %>%
    filter(!is.na(.data$qval)) %>%
    # loci information is encoded in tss_id, which is in the target_id column
    separate(
      .data$target_id, into=c("seqnames", "start", "end", "strand"), sep=",",
      remove=FALSE, convert=TRUE) %>%
    { GRanges(
        seqnames=Rle(.$seqnames),
        ranges=IRanges(.$start, .$end),
        strand=.$strand,
        tss_id=.$target_id,
        qval=.$qval,
        mean_obs=.$mean_obs) }
}
