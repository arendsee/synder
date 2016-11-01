#' @useDynLib synder
#' @importFrom Rcpp sourceCpp
NULL

#' print all blocks with contiguous set ids
#'
#' @param filename synteny map file name
#' @export
dump <- function(filename) {
  if(file.exists(filename)) {
    c_dump(filename)
  } else {
    warning("Failed to open file")
    NULL
  }
}

#' remove links that disagree with the synteny map
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
filter <- function(synfilename, intfilename) {
  if(file.exists(synfilename) && file.exists(intfilename)){
    d <- c_filter(synfilename, intfilename)     %>%
      sub(pattern="\n", replacement="")         %>%
      strsplit(split="\t")                      %>%
      do.call(what=rbind)                       %>%
      as.data.frame(stringsAsFactors=FALSE)
    names(d)[1:6] <- c("qseqid", "qstart", "qstop", "tseqid", "tstart", "tstop")
    d$qstart <- as.numeric(d$qstart)
    d$qstop  <- as.numeric(d$qstop)
    d$tstart <- as.numeric(d$tstart)
    d$tstop  <- as.numeric(d$tstop)
    d
  } else {
    warning("Failed to open file")
    NULL
  }
}

#' trace intervals to target search intervals
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
search <- function(synfilename, gfffilename) {
  if(file.exists(synfilename) && file.exists(gfffilename)){
    c_search(synfilename, gfffilename)
  } else {
    warning("Failed to open file")
    NULL
  }
}


#' trace intervals across genomes
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
map <- function(synfilename, gfffilename) {
  if(file.exists(synfilename) && file.exists(gfffilename)){
    c_map(synfilename, gfffilename)
  } else {
    warning("Failed to open file")
    NULL
  }
}

#' count overlaps
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
count <- function(synfilename, gfffilename) {
  if(file.exists(synfilename) && file.exists(gfffilename)){
    c_count(synfilename, gfffilename)
  } else {
    warning("Failed to open file")
    NULL
  }
}
