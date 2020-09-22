#' Group kallisto results by promoter
#'
#' @param sample_to_variates A \code{data.frame} in the same style as sleuth_prep.
#' @param transcripts_gtf_file File path or URL to a transcript GTF file in the
#'   style of ENSEMBL.
#' @param bedtools_genome_index_file File path or URL to a "genome file" as
#'   specified by bedtools.
#'
#' @importFrom assertthat assert_that is.string has_name
#' @importFrom RCurl url.exists
group_results_by_promoter <- function(
    samples_to_covariates,
    transcripts_gtf_file="http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz",
    bedtools_genome_index_file="https://raw.githubusercontent.com/arq5x/bedtools2/master/genomes/human.hg38.genome") {

  assert_that(is.data.frame(samples_to_covariates))
  assert_that(has_name(samples_to_covariates, "sample"))
  assert_that(has_name(samples_to_covariates, "path"))
  assert_that(is.string(transcripts_gtf_file))
  assert_that(is.string(bedtools_genome_index_file))

  transcripts_gtf_file <- download_if_url(transcripts_gtf_file)
  bedtools_genome_index_file <- download_if_url(bedtools_genome_index_file)
}
