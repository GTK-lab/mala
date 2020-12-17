#' Get a mapping of transcripts to transcription start sites.
#'
#' @param upstream Nucleotide distance for how far upstream of a transcript
#'   should be considered as part of its transcription start site.
#' @param downstream Nucleotide distance for how far downstream of a transcript
#'   should be considered as part of its transcription start site.
#' @param annotation_gtf_file File path or URL to a GTF annotation file in the
#'   style of ENSEMBL.
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
#' @importFrom assertthat assert_that is.string is.count
#' @importFrom dplyr filter select mutate rename
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom rtracklayer import
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr unite
get_tx2tss <- function(
    upstream=500, downstream=100,
    annotation_gtf_file="http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz") {

  assert_that(is.count(upstream) && is.count(downstream))
  assert_that(is.string(annotation_gtf_file))
  annotation_gtf_file <- download_if_url(annotation_gtf_file)

  # From annotation_gtf_file, which has about 3 billion records and 27 fields,
  # filter for only transcript records, and only the fields we care about. We
  # end up with about 230,000 records.
  tx <- annotation_gtf_file %>%  # GTF is 1-based and inclusive at both ends.
    import() %>%
    as_tibble() %>%  # also 1-based and inclusive at both ends.
    filter(.data$type == "transcript") %>%
    select(
      .data$seqnames, .data$start, .data$end, .data$width, .data$strand,
      .data$transcript_id) %>%
    { GRanges(
        seqnames=.$seqnames, ranges=IRanges(start=.$start, end=.$end),
        strand=.$strand, mcols=DataFrame(transcript_id=.$transcript_id)) } %>%
    # Get only the 1nt TSS. GenomicRanges::resize is strand-aware.
    resize(1) %>%
    # Extend the TSS into (sort-of) promoter regions
    promoters(upstream=upstream, downstream=downstream)

  # GenomicRanges::reduce merges!
  tss <- reduce(tx)

  # Now, we are ready to create the target mapping table. For each transcript,
  # we encode its TSS by using an ID which looks like "chr1,69055,70108,+"
  tx2tss_hits <- findOverlaps(tx, tss)
  tss_ordered <- tss[subjectHits(tx2tss_hits)]

  tx[queryHits(tx2tss_hits)] %>%
    mcols %>%
    as_tibble() %>%
    rename(target_id=.data$mcols.transcript_id) %>%
    mutate(
      seqnames=paste0("chr", as.character(seqnames(tss_ordered))),
      start=start(tss_ordered), end=end(tss_ordered),
      strand=as.character(strand(tss_ordered))) %>%
    unite("tss_id", "seqnames", "start", "end", "strand", sep=",")
}
