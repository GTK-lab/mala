% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_gr.R
\name{aggregate_gr}
\alias{aggregate_gr}
\title{Aggregate GRanges-represented sleuth results using TF to loci mappings}
\usage{
aggregate_gr(
  gr,
  tf2loci,
  intermediates = FALSE,
  weight_fn = function(mcols) (1/mcols$tfs_overlapped)/sum(1/mcols$tfs_overlapped),
  filter_mapping_fn = function(n_rows) TRUE
)
}
\arguments{
\item{gr}{GRanges representation of sleuth results obtained using
\code{sr2gr}.}

\item{tf2loci}{List of GRanges representing transcription factor to loci
bindings, obtained using \code{get_tf2loci}.}

\item{intermediates}{If this function should return a list containing
intermediate data structures used in the computation.}

\item{weight_fn}{Unary function to produce weights for Lancaster p-value
aggregation. This function will receive \code{mcols(gr)}, and should return
a numeric vector of the same length as \code{length(gr)}. mcols(gr) can be
an empty DataFrame if a TF overlaps nothing in \code{gr}! Also see details.}

\item{filter_mapping_fn}{Unary function to filter out TSS to TF mappings
which seem poorly supported by \code{tf2loci}. This function will receive
a numeric vector which represents the number of rows in tf2loci which
support each TSS to TF mapping, for each TSS (that is, each numeric vector
contains TSS to TF mappings for one TSS only). Then, the function should
return a logical vector representing which TSS to TF mappings should be
used. E.g. to keep all TSS to TF mappings with more than a median number of
supporting BED entries in the UniBind database, use
\code{filter_mapping_fn=function(x) x > median(x)}.}
}
\value{
If \code{intermediates=FALSE}, a tibble of results, one row per TF.
If \code{intermediates=TRUE}, a list with named elements:
\describe{
\item{gr}{The GRanges of the argument \code{gr}, possibly with an added
metadata columns "tfs_overlapped" --- see details.}
\item{grl}{A list of GRanges. Each GRanges is derived from overlapping
the argument \code{gr} with a transcription factor from the argument
\code{tf2loci}.}
\item{tfs}{The tibble of results which would have been returned by itself
if \code{intermediates=FALSE}}.
}
}
\description{
Aggregate GRanges-represented sleuth results using TF to loci mappings
}
\details{
This function will add a column "tfs_overlapped" to \code{mcols(gr)}
which is the number of transcription factors which each TSS in \code{gr}
maps to, if that column does not yet exists. So, you can use that column
in \code{weight_fn}.
}
\references{
Yi, L., Pimentel, H., Bray, N.L. et al. Gene-level
differential analysis at transcript-level resolution. Genome Biol
19, 53 (2018). https://doi.org/10.1186/s13059-018-1419-z
}
