#' Get a mapping of transcription factors to genomic loci
#'
#' @param unibind_bed_dir File path or URL to a UniBind "TFBSs per TF (in BED
#'   format)" directory or gzipped tarball (ending in \code{tar.gz}).
#' @param save_to_cache Should the raw archive downloaded from UniBind be cached?
#'   \code{overwrite_cache} modifies this behaviour.
#' @param overwrite_cache If there was a previous download from the same
#'   \code{unibind_bed_dir} with \code{save_to_cache}, should this new function
#'   call overwrite that cached download? Note that different
#'   \code{unibind_bed_dir} files are cached differently (they are
#'   version-aware).
#' @param read_from_cache If there was a previous download from the same
#'   \code{unibind_bed_dir}, should this new function call just read from that,
#'   rather than try to download and parse it again?
#'
#' @return A GRangesList. The names of the list are transcription factor
#'   names, and each element of the list is a GRanges recording the loci that
#'   that transcription factor maps to, according to the \code{unibind_bed_dir}.
#'
#' @details \code{unibind_bed_dir} will be used directly if it is a local
#'   directory. Otherwise, if it is a local gzipped tarball (determined by
#'   checking for the file extension \code{tar.gz}), it will be extracted to a
#'   local directory. If it is a remote (part of a URL) gzipped tarball, it will
#'   downloaded then extracted to a local directory.
#'
#'   All UniBind TFBSs per TF BED files are available at their website,
#'   https://unibind.uio.no/downloads/.
#'
#' @importFrom GenomicRanges GRangesList
#' @importFrom assertthat assert_that is.string is.flag
#' @importFrom digest digest
#' @importFrom fs is_dir path dir_create file_exists
#' @importFrom magrittr %>%
#' @importFrom readr read_rds write_rds
#' @importFrom rtracklayer import
#' @importFrom tools R_user_dir
#' @importFrom utils packageName
get_tf2loci <- function(
    unibind_bed_dir="https://unibind.uio.no/static/data/20210421/bulk_Robust/Homo_sapiens/damo_hg38_TFBS_per_TF.tar.gz",
    save_to_cache=TRUE,
    overwrite_cache=FALSE,
    read_from_cache=TRUE
  ) {

  assert_that(is.string(unibind_bed_dir))
  assert_that(is.flag(save_to_cache))
  assert_that(is.flag(overwrite_cache))
  assert_that(is.flag(read_from_cache))

  # Does cache file already exist? Do we read_from_cache?
  cache_dir <- R_user_dir(packageName(), "cache") %>% dir_create()
  cache_fn <- path(
    cache_dir,
    paste(
      "tf2loci-",
      digest(unibind_bed_dir, "md5", serialize=TRUE) %>% substr(1, 8),
      ".rds.gz",
      sep=""))
  if (read_from_cache && file_exists(cache_fn)) {
    message("Using saved tf2loci at ", cache_fn, "...")
    return(read_rds(cache_fn))
  }

  # The download_if_url and extract_if_targz functions perform no operation if
  # the file appears to already have been downloaded or extracted, respectively.
  # This logic is redundant (and indeed, was written before) the caching logic.
  unibind_bed_dir <- unibind_bed_dir %>%
    download_if_url() %>%
    extract_if_targz()

  # CAREFUL: different versions of the UniBind archives are structured differently.
  # For most of 2020, it was     damo_hg38_TFBS_per_TF/TFBS_per_TF/*/*.bed
  # For 20200918, it was         damo_hg38_TFBS_per_TF/Homo_sapiens/TFBS_per_TF/*/*.bed
  # For 20210421, it reverted to damo_hg38_TFBS_per_TF/TFBS_per_TF/*/*.bed
  dir_table <- dir(unibind_bed_dir)
  dir_of_tfs <- path(unibind_bed_dir, dir_table)
  tfs <- dir(dir_of_tfs) %>%
    Filter(function(x) is_dir(path(dir_of_tfs, x)), .)  # "remove" file TFs.txt

  message(sprintf(
    "Loading. Each dot represents a transcription factor (total: %d).",
    length(tfs)))
  # /FOR EACH/ transcription factor
  list_of_granges <- sapply(tfs, function(tf) {
    message(".", appendLF=FALSE)
    dirp <- path(dir_of_tfs, tf)

    if (length(dir(dirp)) == 0) {
      warning(dirp, " is empty. Skipping")
      return(NULL)  # NULLs are explicitly filtered out later.
    }

    # /FOR EACH/ BED file, convert to GRanges. BED is 0-based start-inclusive
    # but end-exclusive, GRanges in 1-based inclusive, but rtracklayer takes
    # care of the mess.
    granges <- lapply(dir(dirp), function(bedfn) {
        granges_sub <- import(path(dirp, bedfn))
      }) %>%
      # Different BED files include different 'esoteric' sequences (chromosomes)
      # and this generates a warning from `c` when combining GRanges.
      { suppressWarnings(do.call(c, .)) }

    granges
  }, simplify=FALSE)
  message()  # just for the newline

  not_nulls <- list_of_granges %>%
    sapply(is.null) %>%
    `!`

  granges_list <- GRangesList(list_of_granges[not_nulls])
  names(granges_list) <- tfs[not_nulls]

  # not save_to_cache => don't save
  if (!save_to_cache) return(granges_list)
  # save_to_cache but cache_fn exists and not overwrite_cache => don't save
  if (save_to_cache && file_exists(cache_fn) && !overwrite_cache) {
    message(
      "You indicated save_to_cache=TRUE, but a cache file already exists. ",
      "To overwrite, set overwrite_cache=TRUE.")
    return(granges_list)
  }
  # in any other case, save to cache

  message("Caching tf2loci to ", cache_fn, "...")
  write_rds(granges_list, cache_fn)  # returns granges_list invisibly
}
