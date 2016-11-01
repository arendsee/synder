#' @useDynLib synder
#' @importFrom Rcpp sourceCpp
NULL

# Changes data.frames to temporary files
df2file <- function(x) {
  if(!is.null(x) && 'data.frame' %in% class(x)){
    xfile <- tempfile()
    write.table(x, file=xfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    x <- xfile
    class(x) <- append(class(x), 'tmp')
  }
  x
}

wrapper <- function(FUN, x, y=NULL, ...) {
  x <- df2file(x)
  y <- df2file(y)

  if(file.exists(x) && is.null(y)){
    d <- FUN(x, ...)
  } else if (file.exists(x) && file.exists(y)) {
    d <- FUN(x, y, ...)
  } else {
    warning("Failed to open files")
    d <- NULL
  }

  # Remove temporary files
  if('tmp' %in% class(x)) file.remove(x)
  if('tmp' %in% class(y)) file.remove(y)

  d
}


#' print all blocks with contiguous set ids
#'
#' @param filename synteny map file name
#' @export
dump <- function(synfile) {
  wrapper(c_dump, synfile)
}

#' remove links that disagree with the synteny map
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
filter <- function(synfilename, intfilename) {
  d <- wrapper(c_filter, synfilename, intfilename)
  if(is.null(d)) return(NULL)

  d <- sub(pattern="\n", replacement="") %>%
    strsplit(split="\t")                 %>%
    do.call(what=rbind)                  %>%
    as.data.frame(stringsAsFactors=FALSE)
  names(d)[1:6] <- c("qseqid", "qstart", "qstop", "tseqid", "tstart", "tstop")
  d$qstart <- as.numeric(d$qstart)
  d$qstop  <- as.numeric(d$qstop)
  d$tstart <- as.numeric(d$tstart)
  d$tstop  <- as.numeric(d$tstop)
  d
}

#' trace intervals to target search intervals
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
search <- function(synfilename, gfffilename) {
  wrapper(c_search, synfilename, gfffilename)
}

#' trace intervals across genomes
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
map <- function(synfilename, gfffilename) {
  wrapper(c_map, synfilename, gfffilename)
}

#' count overlaps
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
count <- function(synfilename, gfffilename) {
  wrapper(c_count, synfilename, gfffilename)
}
