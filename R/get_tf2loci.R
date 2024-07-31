#' Get a mapping of transcription factors to genomic loci
#'
#' @param unibind_fpath File path or URL to a UniBind "TFBSs per TF (in BED format) gzipped tarball (ending in \code{tar.gz}).
#' @param unibind_bed_dir File path to UniBind "TFBSs per TF (in BED
#'   format)" directory. Relative to cache directory if cache is used.
#' @param save_to_cache Should the results be cached? (\code{overwrite_cache}
#'   modifies this behaviour.)
#' @param overwrite_cache If there was a pre-existing cache from the same
#'   \code{unibind_bed_dir} with \code{save_to_cache}, should this new function
#'   call overwrite that cached result? Note that different
#'   \code{unibind_bed_dir} files are cached differently (they are
#'   version-aware).
#' @param read_from_cache If there was a pre-existing cache from the same
#'   \code{unibind_bed_dir}, should this new function call just read from that,
#'   rather than try to download and parse it again?
#'
#' @return A GRangesList. The names of the list are transcription factor
#'   names, and each element of the list is a GRanges recording the loci that
#'   that transcription factor maps to, according to the \code{unibind_bed_dir}.
#'
#' @details \code{unibind_bed_dir} will be used directly if it is a local
#'   directory and exists. Otherwise, if it is a local gzipped tarball (determined by
#'   checking for the file extension \code{tar.gz}), it will be extracted to a
#'   local directory. If it is a remote (part of a URL) gzipped tarball, it will
#'   downloaded then extracted to a local directory.
#'
#'   All UniBind TFBSs per TF BED files are available at their website,
#'   https://unibind.uio.no/downloads/.
#'
#' @importFrom BiocFileCache bfcadd bfcquery bfcneedsupdate BiocFileCache
#' @importFrom GenomicRanges GRangesList
#' @importFrom assertthat assert_that is.string is.flag
#' @importFrom digest digest
#' @importFrom fs is_dir path dir_create file_exists file_delete
#' @importFrom glue glue
#' @importFrom options opt
#' @importFrom purrr keep
#' @importFrom readr read_rds write_rds
#' @importFrom rtracklayer import
#' @importFrom stringr str_ends
#' @importFrom tools R_user_dir
#' @importFrom utils packageName

get_tf2loci <- function(
                        unibind_fpath=opt("unibind-fpath"),
                        unibind_rname=opt("unibind-rname"),
                        overwrite_cache=FALSE,
                        bfc=BiocFileCache()
                        ) {

    assert_that(is.string(unibind_fpath))
    assert_that(is.string(unibind_rname))
    assert_that(is.flag(overwrite_cache))

    rid <- bfcquery(bfc,unibind_fpath,field='fpath',exact=TRUE)$rid

    cache_ok <- (length(rid) == 1L) && (!bfcneedsupdate(bfc,rids=rid))
    if (cache_ok) {
        unibind_tarball <- as.character(bfc[[rid]])
    } else {
        unibind_tarball <- bfcadd(bfc,rname=unibind_rname,fpath=unibind_fpath)
    }
    cache_dir <- R_user_dir(ifelse(is.character(packageName()),packageName(),"mala"), "cache") |> dir_create()

    message(glue("cach_dir is {cache_dir}"))

    cache_fn <- fs::path(
        cache_dir,
        paste(
            "tf2loci-",
            digest(unibind_rname, "md5", serialize=TRUE) |> substr(1, 8),
            ".rds.gz",
            sep=""))

    message(glue("cache_fn is {cache_fn}"))
    

    if (cache_ok && file_exists(cache_fn)) {
        message("Using saved tf2loci at ", cache_fn, "...")
        return(read_rds(cache_fn))
    } else {
        if  (file_exists(cache_fn)) {
            file_delete(cache_fn)
        }
        exdir=tempdir()
        bedfiles <- keep(untar(unibind_tarball,list=TRUE),str_ends,".bed")
        success_p <- untar(unibind_tarball,exdir=exdir,verbose=TRUE)        
        bedfiles_by_tf <- split(bedfiles,sapply(bedfiles,str_split_i,"/",2))        

        message(sprintf(
            "Loading. Each dot represents a transcription factor (total: %d).",
            length(bedfiles_by_tf)))
                                        # /FOR EACH/ transcription factor
        granges_by_tf <- GRangesList(sapply(names(bedfiles_by_tf), function(tf) {
            message(glue("{tf}-"), appendLF=FALSE)
            beds <- keep(paste(exdir,bedfiles_by_tf[[tf]],sep="/"),file_exists)

            if (length(beds) == 0) {
                warning(glue("No bed files for {tf}. Skipping"))
                return(NULL)  # NULLs are explicitly filtered out later.
            }

            return(keepStandardChromosomes(
                reduce(unlist(GRangesList(lapply(beds, import)))),
                pruning.mode = 'coarse'))

        }, simplify=FALSE))
        message("Caching tf2loci to ", cache_fn, "...")
        write_rds(granges_by_tf, cache_fn,compress = "gz")
    }
    
}
