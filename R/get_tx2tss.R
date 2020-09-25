#' Get a mapping of transcripts to transcription start sites.
#'
#' @param threshold If two transcripts are within \code{threshold} nucleotides
#'   of each other, they are grouped as having the same transcription start sites.
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
#' @importFrom assertthat assert_that is.string is.count
#' @importFrom magrittr %>%
#' @importFrom rtracklayer import
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr filter select arrange mutate transmute lag rename left_join
#'   bind_rows
#' @importFrom rlang .data .env
get_tx2tss <- function(
    threshold=500,
    annotation_gtf_file="http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz") {

  assert_that(is.count(threshold))
  assert_that(is.string(annotation_gtf_file))
  annotation_gtf_file <- download_if_url(annotation_gtf_file)

  # From annotation_gtf_file, which has about 3 billion records and 27 fields,
  # filter for only transcript records, and only the fields we care about. We
  # end up with about 230,000 records.
  gtf <- annotation_gtf_file %>%  # GTF is 1-based and inclusive at both ends.
    import() %>%
    as_tibble() %>%  # also 1-based and inclusive at both ends.
    filter(.data$type == "transcript") %>%
    select(
      .data$seqnames, .data$start, .data$end, .data$width, .data$strand,
      .data$transcript_id)

  # /FOR EACH/ chromosome (indicated by seqnames), group transcripts into
  # TSS using the value of threshold.
  message("Finding TSS for seqnames ", appendLF=FALSE)
  tx2tss <- lapply(unique(gtf$seqnames), function(seqname) {
      message(seqname, " ", appendLF=FALSE)
      gtf_chr <- filter(gtf, .data$seqnames == seqname)

      # /FOR EACH/ strand (we must be strand-aware!)...
      lapply(unique(gtf_chr$strand), function(strand) {
          assert_that(strand %in% c("+", "-"))
          # First, group tx into tss simply with integer ids e.g. 1, 2, 3...
          tb_simple_id <- filter(gtf_chr, .data$strand == .env$strand) %>%
            arrange(.data$start) %>%
            mutate(
              # start pos of tx just upstream
              start_prev = lag(.data$start, default=0),
              # dist of start from tx just upstream
              dist = .data$start - .data$start_prev,
              is_new_tss = .data$dist > threshold,
              tss_id = cumsum(.data$is_new_tss))

          # /FOR EACH/ TSS, we must assign a name which will be unique across
          # all strands of all chromosomes.
          new_ids <- lapply(unique(tb_simple_id$tss_id), function(tss_id) {
              tss <- filter(tb_simple_id, .data$tss_id == .env$tss_id)
              tibble(tss_id=tss_id, new_tss_id=sprintf(
              "chr%s,%d,%d,%s", seqname, min(tss$start), max(tss$end), strand))
            }) %>%
            do.call(bind_rows, .)  # tibble of tss_id, new_tss_id, for this strand

          # Finally, replace the original and simple ids in tb_simple_id with the
          # new, unique ones. Also, just keep the two columns we need:
          # transcript_id, and the new_tss_id.
          tb_simple_id %>%
            left_join(new_ids, by="tss_id") %>%
            transmute(transcript_id=.data$transcript_id, tss_id=.data$new_tss_id)

        }) %>%
        do.call(bind_rows, .)  # tibble of transcript_id, tss_id, for this chr

    }) %>%
    do.call(bind_rows, .)  # tibble of transcript_id, tss_id.
  message()  # just for the newline

  # for compatibility with sleuth
  tx2tss %>% rename(target_id=.data$transcript_id)
}
