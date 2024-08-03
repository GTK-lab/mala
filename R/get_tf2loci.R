#' Get a mapping of transcription factors to genomic loci
#'
#' @param unibind_fpath File path or URL to a UniBind "TFBSs per TF (in BED format) gzipped tarball
#' (ending in \code{tar.gz}).
#' @param unibind_rname Resource name to be used for the Unibind or equivalent set
#' @param overwrite_cache If there was a pre-existing cache from the same
#'   \code{unibind_bed_dir} with \code{save_to_cache}, should this new function
#'   call overwrite that cached result? Note that different
#'   \code{unibind_bed_dir} files are cached differently (they are
#'   version-aware).
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
#'

get_tf2loci <- function(
                        unibind_fpath=opt("unibind-fpath"),
                        unibind_rname=opt("unibind-rname"),
                        overwrite_cache=FALSE
                        ) {

    assert_that(is.string(unibind_fpath))
    assert_that(is.string(unibind_rname))
    assert_that(is.flag(overwrite_cache))


    main_bfc <- BiocFileCache()
    mala_bfc <- R_user_dir(ifelse(is.character(packageName()),packageName(),"mala"), "cache") |> dir_create() |> BiocFileCache()

    unibind_rid <- bfcquery(main_bfc, unibind_fpath, field="fpath", exact=TRUE)$rid

    cache_ok <- (length(unibind_rid) == 1L) && (!bfcneedsupdate(main_bfc, rids=unibind_rid))
    if (cache_ok) {
        message(glue("using existing unibind cache file {main_bfc[[unibind_rid]]}"))
        unibind_tarball <- as.character(main_bfc[[unibind_rid]])
    } else {
        message(glue("caching unibind data from {unibind_fpath}"))
        unibind_tarball <- bfcadd(main_bfc,rname=unibind_rname,fpath=unibind_fpath)
    }

    message(glue("package cache directory is {bfccache(mala_bfc)}"))
    tf2loci_rid <- bfcquery(mala_bfc, unibind_rname, field="rname", exact=TRUE)$rid
    
    mala_cache_ok <- (length(tf2loci_rid) == 1L) 

    if (mala_cache_ok) {
        conn <- dbConnect(RSQLite::SQLite(),mala_bfc[[tf2loci_rid]])
        tf2loci_db <- UnibindDb(mala_bfc[[tf2loci_rid]])
    }  else {
        message("creating new tf2loci database")
        tf2loci_fn <- bfcnew(mala_bfc,unibind_rname,ext=".sqlite")
        tf2loci_db <- UnibindDb(tf2loci_fn)
        conn <- tf2loci_db@conn

        bedfiles <- keep(untar(unibind_tarball,list=TRUE),str_ends,".bed")
        bedfiles_by_tf <- split(bedfiles,sapply(bedfiles,str_split_i,"/",2))
        exdir=tempdir()

        untar(unibind_tarball,exdir=exdir,verbose=TRUE)
    sapply(names(bedfiles_by_tf), function(tf) {
        message(glue("{tf}-"), appendLF=FALSE)
        beds <- keep(paste(exdir,bedfiles_by_tf[[tf]],sep="/"),file_exists)
        grl <- lapply(beds, function(fn) {
            gr <-keepStandardChromosomes(import(fn),pruning.mode='coarse')
            gr$bedfile <- str_split_i(fn,"/",-1)
            gr  })
        tbl <- as.data.frame(GenomicRanges::reduce(unlist(GRangesList(grl))))
        dbWriteTable(conn,tf,tbl)
    })
    message("\ncompleted") ## for the new line

    }

    return(tf2loci_db)


        


    
}
