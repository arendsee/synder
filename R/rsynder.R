#' @useDynLib synder
#' @importFrom Rcpp sourceCpp
NULL

SYNMAP_COLS <- c(
  "qconid" = "character",
  "qstart" = "integer",
  "qstop"  = "integer",
  "tconid" = "character",
  "tstart" = "integer",
  "tstop"  = "integer",
  "score"  = "numeric",
  "strand" = "character"
)

GFF_COLS <- c(
  "conid"   = "character",
  "source"  = "character",
  "type"    = "character",
  "start"   = "integer",
  "stop"    = "integer",
  "score"   = "numeric",
  "strand"  = "character",
  "phase"   = "integer",
  "seqname" = "character"
)

CON_LENGTH <- c(
  "conid"  = "character",
  "length" = "integer"
)

defaults <- list(
  tclfile = "",
  qclfile = "",
  offsets = c(0,1,0,0,1,1),
  k       = 0,
  r       = 0.001,
  trans   = 'i',
  swap    = FALSE
)

# NOTE: Handling of offsets:
# --------------------------
# Start and stop base offsets (0 or 1) vary between between formats.  The main
# synteny build I use, Satsuma, is 0-based relative to start, and 1-based
# relative to stop. GFF files are required to be 0-based. bioconductor (and R
# in general) is 1-based. Internally, Synder is 0-based. There are 3 sets of
# start/stop offsets to consider: synteny map base, input GFF or hitmap base,
# and output base. I will let C-side synder handle input bases, and R-side
# Synder handle output base.

do_offsets <- function(d, offsets){
  d$qstart <- d$qstart + offsets[5]
  d$qstop  <- d$qstop  + offsets[6]
  d$tstart <- d$tstart + offsets[5]
  d$tstop  <- d$tstop  + offsets[6]
  d
}

#' Read a synteny map
#'
#' @param file Synteny map file name
#' @return A dataframe
#' @export
read_synmap <- function(file) {
  d <- readr::read_tsv(file, col_names=FALSE)

  # Assert the correct number of columns were read
  stopifnot(ncol(d) == 8)  

  # Set column names
  names(d) <- names(SYNMAP_COLS)

  # Set types (readr guesses the others correctly)
  d$score <- as.numeric(d$score)

  # Assert column types match expectations
  stopifnot(SYNMAP_COLS == lapply(d, class))

  # Assign class membership
  class(d) <- append('synmap', class(d))

  return(d)
}

#' Read a GFF file 
#'
#' @param file GFF file name
#' @return A dataframe
#' @export
read_gff <- function(file) {
  d <- readr::read_tsv(file, col_names=FALSE)
  class(d) <- 'data.frame'

  # Assert the correct number of columns were read
  stopifnot(ncol(d) == 9)

  # Set column names
  names(d) <- names(GFF_COLS)

  # Set types (readr guesses the others correctly)
  d$score <- as.numeric(d$score)
  d$phase <- as.integer(d$phase)

  # Assert column types match expectations
  stopifnot(GFF_COLS == lapply(d, class))

  # Assign class membership
  class(d) <- append('gff', class(d))

  return(d)
}

#' Read a hit file 
#'
#' A hit file is required to have the same first 6 columns as a synteny map,
#' i.e. qcon, qstart, qstop, tcon, tstart, tstop. There can be any number of
#' additional columns.
#'
#' @param file hit file name
#' @return A dataframe
#' @export
read_hitmap <- function(file) {
  d <- readr::read_tsv(file, col_names=FALSE)

  # Assert the correct number of columns were read
  stopifnot(ncol(d) >= 6)

  # Set column names
  names(d)[1:6] <- names(SYNMAP_COLS)[1:6]

  # Assert column types match expectations
  stopifnot(SYNMAP_COLS[1:6] == lapply(d[,1:6], class))

  # Assign class membership
  class(d) <- append('hitmap', class(d))

  return(d)
}

#' Read a contig length file 
#'
#' @param file contig length file
#' @return A dataframe
#' @export
read_conlen <- function(file) {
  d <- readr::read_tsv(file, col_names=FALSE)

  # Assert the correct number of columns were read
  stopifnot(ncol(d) == 2)

  # Set column names
  names(d) <- names(CON_LENGTH)

  # Assert column types match expectations
  stopifnot(CON_LENGTH == lapply(d, class))

  # Assign class membership
  class(d) <- append('conlen', class(d))

  return(d)
}

# Changes data.frames to temporary files
df2file <- function(x) {
  if(!is.null(x) && 'data.frame' %in% class(x)){
    xfile <- tempfile()
    readr::write_tsv(x, path = xfile, col_names = FALSE)
    x <- xfile
    class(x) <- append(class(x), 'tmp')
  }
  x
}

check_parameters <- function(
  offsets = NULL,
  k       = NULL,
  r       = NULL,
  swap    = NULL,
  trans   = NULL,
  ...
){
  stopifnot(is.null(offsets) || all(offsets %in% c(1,0)))
  stopifnot(is.null(k)       || is.numeric(k))
  stopifnot(is.null(r)       || is.numeric(r))
  stopifnot(is.null(swap)    || is.logical(swap))
  stopifnot(is.null(trans)   || trans %in% c('i', 'd', 'p', 'l'))
}

wrapper <- function(FUN, x, y=NULL, ...) {

  x <- df2file(x)
  y <- df2file(y)

  check_parameters(...)

  if(is.null(x)){
    d <- NULL
    warning("x is NULL")
  } else if(file.exists(x) && is.null(y)){
    d <- FUN(x, ...) %>% tibble::as_data_frame()
  } else if (file.exists(x) && file.exists(y)) {
    d <- FUN(x, y, ...) %>% tibble::as_data_frame()
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
dump <- function(
  synfile,
  swap    = defaults$swap,
  trans   = defaults$trans,
  offsets = defaults$offsets
) {
  syn <- wrapper(
    FUN     = c_dump,
    x       = synfile,
    swap    = swap,
    trans   = trans,
    offsets = offsets[1:4]
  )
  syn$qcon <- as.character(syn$qcon)
  syn$tcon <- as.character(syn$tcon)

  # Assign class and attributes
  class(syn) <- append('dump_result', class(syn))
  attributes(syn)$swap  = swap
  attributes(syn)$trans = trans

  syn <- do_offsets(syn, offsets)

  syn
}

#' remove links that disagree with the synteny map
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
filter <- function(
  synfilename,
  intfilename,
  swap    = defaults$swap,
  trans   = defaults$trans,
  k       = defaults$k,
  r       = defaults$r,
  offsets = defaults$offsets
) {
  d <- wrapper(
    FUN     = c_filter,
    x       = synfilename,
    y       = intfilename,
    swap    = swap,
    k       = k,
    r       = r,
    trans   = trans,
    offsets = offsets[1:4]
  )
  if(is.null(d)) return(NULL)

  d <- sub(d, pattern="\n", replacement="") %>%
    strsplit(split="\t")                    %>%
    do.call(what=rbind)                     %>%
    tibble::as_data_frame()
  names(d)[1:6] <- names(SYNMAP_COLS)[1:6]
  d$qstart <- as.numeric(d$qstart)
  d$qstop  <- as.numeric(d$qstop)
  d$tstart <- as.numeric(d$tstart)
  d$tstop  <- as.numeric(d$tstop)

  class(d) <- append('filter_result', class(d))
  attributes(d)$swap  <- swap
  attributes(d)$k     <- k
  attributes(d)$r     <- r
  attributes(d)$trans <- trans

  d <- do_offsets(d, offsets)

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
  tclfile = defaults$tclfile,
  qclfile = defaults$qclfile,
  swap    = defaults$swap,
  trans   = defaults$trans,
  k       = defaults$k,
  r       = defaults$r,
  offsets = defaults$offsets
) {
  d <- wrapper(
    FUN     = c_search,
    x       = synfilename,
    y       = gfffilename,
    tclfile = tclfile,
    qclfile = qclfile,
    swap    = swap,
    k       = k,
    r       = r,
    trans   = trans,
    offsets = offsets[1:4]
  )

  d$seqname <- as.character(d$seqname)
  d$qcon    <- as.character(d$qcon)
  d$tcon    <- as.character(d$tcon)
  d$strand  <- as.character(d$strand)

  class(d) <- append('search_result', class(d))
  attributes(d)$swap  <- swap
  attributes(d)$k     <- k
  attributes(d)$r     <- r
  attributes(d)$trans <- trans

  d <- do_offsets(d, offsets)

  d
}

#' trace intervals across genomes
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
map <- function(
  synfilename,
  gfffilename,
  swap    = defaults$swap,
  offsets = defaults$offsets
) {
  d <- wrapper(
    FUN     = c_map,
    x       = synfilename,
    y       = gfffilename,
    swap    = swap,
    offsets = offsets[1:4]
  )

  d$seqname <- as.character(d$seqname)
  d$qcon    <- as.character(d$qcon)
  d$tcon    <- as.character(d$tcon)
  d$strand  <- as.character(d$strand)

  class(d) <- append('map_result', class(d))
  attributes(d)$swap <- swap

  d <- do_offsets(d, offsets)

  d
}

#' count overlaps
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
count <- function(
  synfilename,
  gfffilename,
  swap    = defaults$swap,
  offsets = defaults$offsets
) {
  d <- wrapper(
    FUN     = c_count,
    x       = synfilename,
    y       = gfffilename,
    swap    = swap,
    offsets = offsets[1:4]
  )

  d$seqname = as.character(d$seqname)

  class(d) <- append('count_result', class(d))
  attributes(d)$swap <- swap

  d
}
