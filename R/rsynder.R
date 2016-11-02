#' @useDynLib synder
#' @importFrom Rcpp sourceCpp
NULL

# Changes data.frames to temporary files
df2file <- function(x) {
  if(!is.null(x) && 'data.frame' %in% class(x)){
    xfile <- tempfile()
    wite.table(
      x,
      file      = xfile,
      quote     = FALSE,
      sep       = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
    x <- xfile
    class(x) <- append(class(x), 'tmp')
  }
  x
}

wrapper <- function(FUN, x, y=NULL, ...) {
  x <- df2file(x)
  y <- df2file(y)

  if(is.null(x)){
    d <- NULL
    warning("x is NULL")
  } else if(file.exists(x) && is.null(y)){
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
dump <- function(synfile, swap=FALSE) {
  wrapper(c_dump, synfile, swap=swap)
}

#' remove links that disagree with the synteny map
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
filter <- function(synfilename, intfilename, swap=FALSE, k=0, r=0, trans="i") {
  d <- wrapper(
    c_filter,
    synfilename,
    intfilename,
    swap  = swap,
    k     = k,
    r     = r,
    trans = trans
  )
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
search <- function(
  synfilename,
  gfffilename,
  tclfile = "",
  qclfile = "",
  swap    = FALSE,
  k       = 0,
  r       = 0,
  trans   = "i"
) {
  wrapper(
    c_search,
    synfilename,
    gfffilename,
    tclfile = tclfile,
    qclfile = qclfile,
    swap    = swap,
    k       = k,
    r       = r,
    trans   = trans
  )
}

#' trace intervals across genomes
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
map <- function(
  synfilename,
  gfffilename,
  swap=FALSE
) {
  wrapper(c_map, synfilename, gfffilename, swap=swap)
}

#' count overlaps
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
count <- function(
  synfilename,
  gfffilename,
  swap=FALSE
) {
  wrapper(c_count, synfilename, gfffilename, swap=swap)
}
