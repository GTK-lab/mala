#' Get a mapping of transcripts to transcription start sites.
#'
#' @param upstream Nucleotide distance for how far upstream of a transcript
#'   should be considered as part of its transcription start site.
#' @param downstream Nucleotide distance for how far downstream of a transcript
#'   should be considered as part of its transcription start site.
#' @param annotation_gtf_file File path or URL to a GTF annotation file in the
#'   style of ENSEMBL.
#' @param save_to_cache Should the results be cached? (\code{overwrite_cache}
#'   modifies this behaviour.)
#' @param overwrite_cache If there was a pre-existing cache from the same
#'   \code{annotation_gtf_file} with \code{save_to_cache}, should this new
#'   function call overwrite that cached result? Note that different
#'   \code{annotation_gtf_file} files are cached differently (they are
#'   version-aware).
#' @param read_from_cache If there was a pre-existing cache from the same
#'   \code{annotation_gtf_file}, should this new function call just read from
#'   that, rather than try to download and parse it again?
#'
#' @return tibble with two columns: \code{target_id} which contains ENST
#'   transcript accessions, and \code{tss_id} which is a unique ID string for
#'   each trascription start site (TSS). \code{tss_id} is a comma separated
#'   string with chromosome name, start of TSS window, end of TSS window, and
#'   strand (\code{"+"} or \code{"-"}). Note that \code{start <= end} regardless
#'   of strand.
#'
#'   This tibble can be used as \code{sleuth_prep}'s \code{target_mapping}. Just
#'   remebmer to set \code{aggregation_column = "tss_id"}!
#'
#' @importFrom GenomicRanges GRanges resize promoters reduce findOverlaps mcols
#'   seqnames start end strand
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame queryHits subjectHits
#' @importFrom assertthat assert_that is.string is.count is.flag
#' @importFrom dplyr filter select mutate rename
#' @importFrom fs dir_create path file_exists
#' @importFrom magrittr |>
#' @importFrom readr read_rds write_rds
#' @importFrom rlang .data
#' @importFrom rtracklayer import
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr unite
#' @importFrom tools R_user_dir
#' @importFrom utils packageName
get_tx2tss <- function(
    upstream=500, downstream=100,
    annotation_gtf_file="http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz",
    save_to_cache=TRUE,
    overwrite_cache=FALSE,
    read_from_cache=TRUE
  ) {

  assert_that(is.count(upstream) && is.count(downstream))
  assert_that(is.string(annotation_gtf_file))
  assert_that(is.flag(save_to_cache))
  assert_that(is.flag(overwrite_cache))
  assert_that(is.flag(read_from_cache))

  # Does a cache file already exist? Do we read_from_cache?
  cache_dir <- R_user_dir(packageName(), "cache") |> dir_create()
  cache_fn <- path(
    cache_dir,
    paste(
      "tx2tss-u", upstream, "-d", downstream, "-",
      digest(annotation_gtf_file, "md5", serialize=TRUE) |> substr(1, 8),
      ".rds.gz",
      sep=""))
  if (read_from_cache && file_exists(cache_fn)) {
    message("Using saved tx2tss at ", cache_fn, "...")
    return(read_rds(cache_fn))
  }

  annotation_gtf_file <- download_if_url(annotation_gtf_file)

  # From annotation_gtf_file, which has about 3 billion records and 27 fields,
  # filter for only transcript records, and only the fields we care about. We
  # end up with about 230,000 records.
  tx <- annotation_gtf_file |>  # GTF is 1-based and inclusive at both ends.
    import() |>
    as_tibble() |>  # also 1-based and inclusive at both ends.
    filter(.data$type == "transcript") |>
    select(
      .data$seqnames, .data$start, .data$end, .data$width, .data$strand,
      .data$transcript_id) |>
    { GRanges(
        seqnames=.$seqnames, ranges=IRanges(start=.$start, end=.$end),
        strand=.$strand, mcols=DataFrame(transcript_id=.$transcript_id)) } |>
    # Get only the 1nt TSS. GenomicRanges::resize is strand-aware.
    resize(1) |>
    # Extend the TSS into (sort-of) promoter regions
    promoters(upstream=upstream, downstream=downstream)

  # GenomicRanges::reduce merges!
  tss <- reduce(tx)

  # Now, we are ready to create the target mapping table. For each transcript,
  # we encode its TSS by using an ID which looks like "chr1,69055,70108,+"
  tx2tss_hits <- findOverlaps(tx, tss)
  tss_ordered <- tss[subjectHits(tx2tss_hits)]

  results <- tx[queryHits(tx2tss_hits)] |>
    mcols |>
    as_tibble() |>
    rename(target_id=.data$mcols.transcript_id) |>
    mutate(
      seqnames=paste0("chr", as.character(seqnames(tss_ordered))),
      start=start(tss_ordered), end=end(tss_ordered),
      strand=as.character(strand(tss_ordered))) |>
    unite("tss_id", "seqnames", "start", "end", "strand", sep=",")

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
