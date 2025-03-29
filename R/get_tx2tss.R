#' Get a mapping of transcripts to transcription start sites and promoters 
#' (as defined by the user)
#'
#' @param transcripts_ranges A GRanges object for transcripts. These will 
#'   be used to build a map of distinct TSSs and promoters.
#' @param upstream Nucleotide distance for how far upstream of a TSS
#'   should be considered as part of its promoter.
#' @param downstream Nucleotide distance for how far downstream of a TSS
#'   should be considered as part of its promoter.
#' @param save_to_cache Should the results be cached? (\code{overwrite_cache}
#'   modifies this behaviour.)
#' @param overwrite_cache If there was a pre-existing cache from the same
#'   set of transcripts and parameters with \code{save_to_cache}, should this new
#'   function call overwrite that cached result? Note that different
#'   \code{annotation_gtf_file} files are cached differently (they are
#'   version-aware).
#' @param read_from_cache If there was a pre-existing cache from the same
#'   \code{transcript_ranges}, should this new function call just read from
#'   that, rather than try to regenerate it?
#'
#' @return tibble with three columns: \code{transcript_id} which contains
#'   transcript accessions,  \code{tss_id} which is a unique ID string for
#'   each transcription start site (TSS), and \code{promoter_id} for each
#'   promoter region.
#'
#'   This tibble can be used as \code{sleuth_prep}'s \code{target_mapping}. Just
#'   remember to set \code{aggregation_column = "promoter_id"}!
#'
#' @importFrom GenomicRanges GRanges resize promoters reduce findOverlaps mcols
#'   seqnames start end strand
#' @importFrom S4Vectors DataFrame queryHits subjectHits
#' @importFrom assertthat assert_that is.string is.count is.flag
#' @importFrom dplyr filter select mutate rename arrange
#' @importFrom plyranges select mutate
#' @importFrom fs dir_create path file_exists
#' @importFrom options opt
#' @importFrom readr read_rds write_rds
#' @importFrom rlang .data
#' @importFrom rtracklayer import
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr unite unnest
#' @importFrom tools R_user_dir
#' @importFrom utils packageName
#' @export
get_tx2tss <- function(
    transcripts_ranges,
    upstream=500, downstream=100,
    save_to_cache=TRUE,
    overwrite_cache=FALSE,
    read_from_cache=TRUE
  ) {

  assert_that(class(transcripts_ranges) == "GRanges")
  assert_that(is.count(upstream) && is.count(downstream))
  assert_that(is.flag(save_to_cache))
  assert_that(is.flag(overwrite_cache))
  assert_that(is.flag(read_from_cache))

  bfc=BiocFileCache()
  # Does a cache file already exist? Do we read_from_cache?

  cache_dir <- R_user_dir(ifelse(is.character(packageName()),packageName(),"mala"), "cache") |> dir_create()
  cache_fn <- path(
    cache_dir,
    paste(
      "tx2tss-u", upstream, "-d", downstream, "-",
      digest(transcripts_ranges, "md5", serialize=TRUE),
      ".rds.gz",
      sep=""))
  if (read_from_cache && file_exists(cache_fn)) {
    message("Using saved tx2tss at ", cache_fn, "...")
    return(read_rds(cache_fn))
  }

    tss <- plyranges::select(resize(transcripts_ranges,width=1),transcript_id) |>
      plyranges::mutate(tss_id=sprintf("%s:%d:%d:%s",seqnames,start,end,strand))
    # Extend the TSS into (sort-of) promoter regions
    tss_tibble <- dplyr::select(as_tibble(tss),tss_id,transcript_id)
    promoter_regions <- GenomicRanges::reduce(promoters(tss,upstream=upstream, downstream=downstream),
                                              with.revmap=TRUE) |>
      plyranges::mutate(promoter_id=sprintf("%s:%d:%d:%s",seqnames,start,end,strand))
    
    promoter_tibble <- dplyr::select(as_tibble(promoter_regions),promoter_id,revmap) |> 
      unnest(cols=c(revmap)) |> arrange(revmap) 
    
    results <-    as_tibble(cbind(tss_tibble,promoter_tibble)) |> 
      dplyr::select(promoter_id,tss_id,transcript_id)
    
  # not save_to_cache => don't save
  if (!save_to_cache) return(results)
  # save_to_cache but cache_fn exists and not overwrite_cache => don't save
  if (save_to_cache && file_exists(cache_fn) && !overwrite_cache) {
    message(
      "You indicated save_to_cache=TRUE, but a cache file already exists. ",
      "To overwrite, set overwrite_cache=TRUE.")
    return(results)
  }

  # in any other case, save to cache
  message("Caching tx2tss to ", cache_fn, "...")
  write_rds(results, cache_fn)  # returns results invisibly
}
