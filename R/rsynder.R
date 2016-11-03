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
  "strand" = "logical"
)

GFF_COLS <- c(
  "conid"  = "character",
  "source" = "character",
  "type"   = "character",
  "start"  = "integer",
  "stop"   = "integer",
  "score"  = "numeric",
  "strand" = "logical",
  "phase"  = "integer",
  "attr"   = "character"
)

CON_LENGTH <- c(
  "conid"  = "character",
  "length" = "integer"
)

# Converts a character vector of '+' and '-' to a logical vector
# where '+' -> TRUE, '-' -> FALSE, else -> NA
logical_strand <- function(x) {
  as.logical(c_logical_strand(x)) 
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
  d$strand <- logical_strand(d$strand)

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
  d$strand <- logical_strand(d$strand)

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

wrapper <- function(FUN, x, y=NULL, ...) {
  x <- df2file(x)
  y <- df2file(y)

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
dump <- function(synfile, swap=FALSE, trans="i") {
  syn <- wrapper(FUN=c_dump, x=synfile, swap=swap, trans=trans)
  syn$qseqid <- as.character(syn$qseqid)
  syn$tseqid <- as.character(syn$tseqid)

  # Cast strand from int [0,1,NA_INTEGER] to logical
  syn$strand <- as.logical(syn$strand)

  # Assign class and attributes
  class(syn) <- append('dump_result', class(syn))
  attributes(syn)$swap  = swap
  attributes(syn)$trans = trans

  syn
}

#' remove links that disagree with the synteny map
#'
#' @param synfilename synteny map file name
#' @param intfilename int file name
#' @export
filter <- function(synfilename, intfilename, swap=FALSE, k=0, r=0, trans="i") {
  d <- wrapper(
    FUN   = c_filter,
    x     = synfilename,
    y     = intfilename,
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
  d <- wrapper(
    FUN     = c_search,
    x       = synfilename,
    y       = gfffilename,
    tclfile = tclfile,
    qclfile = qclfile,
    swap    = swap,
    k       = k,
    r       = r,
    trans   = trans
  )

  class(d) <- append('search_result', class(d))
  attributes(d)$swap  <- swap
  attributes(d)$k     <- k
  attributes(d)$r     <- r
  attributes(d)$trans <- trans

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
  swap=FALSE
) {
  d <- wrapper(
    FUN  = c_map,
    x    = synfilename,
    y    = gfffilename,
    swap = swap
  )

  class(d) <- append('map_result', class(d))
  attributes(d)$swap <- swap

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
  swap=FALSE
) {
  d <- wrapper(
    FUN  = c_count,
    x    = synfilename,
    y    = gfffilename,
    swap = swap
  )

  class(d) <- append('count_result', class(d))
  attributes(d)$swap <- swap

  d
}
