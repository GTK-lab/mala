#' Get a mapping of transcription factors to genomic loci
#'
#' @param unibind_bed_dir File path or URL to a UniBind "TFBs per TF (in BED
#'   format)" directory or gzipped tarball (ending in \code{tar.gz}).
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
#' @importFrom magrittr %>%
#' @importFrom assertthat assert_that is.string
#' @importFrom fs is_dir path
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRangesList
get_tf2loci <- function(
    unibind_bed_dir="https://unibind.uio.no/static/data/20200918/bulk_Robust/Homo_sapiens/damo_hg38_TFBS_per_TF.tar.gz"
  ) {
  assert_that(is.string(unibind_bed_dir))

  unibind_bed_dir <- unibind_bed_dir %>%
    download_if_url() %>%
    extract_if_targz()

  # extract_if_targz dumps BED files like:
  # damo_hg38_TFBS_per_TF/Homo_sapiens/TFBS_per_TF/*.bed
  dir_species <- dir(unibind_bed_dir)
  dir_table <- dir(path(unibind_bed_dir, dir_species))
  dir_of_tfs <- path(unibind_bed_dir, dir_species, dir_table)
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

  granges_list
}
