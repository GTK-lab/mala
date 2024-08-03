setClass("UnibindObj",
         contains="VIRTUAL",
         slots=list(
             objName="character",
             objTarget="character"     
         ))


setClass("UnibindDb",
         contains="UnibindObj",
         slots=list(
             dbfile="character",
             conn="SQLiteConnection"
         ))




UnibindDb <- function(dbfile) {
    new("UnibindDb",
        dbfile=dbfile,
        conn=DBI::dbConnect(RSQLite::SQLite(),dbfile),
        objName=basename(dbfile))

}

setGeneric("TFs",
           function(db) standardGeneric("TFs"))
setMethod("TFs","UnibindDb",
          function(db) dbListTables(db@conn))

setMethod("[[",c("UnibindDb","character","missing"),
          
          function(x,i,j) {
              if ((length(i) == 1L) &&  (i %in% TFs(x))) {
                  GenomicRanges::GRanges(dbReadTable(x@conn,i))
              } else {
                  message(glue("{i} is not a transcription factor in {ub@objName}"))
              }
          })

