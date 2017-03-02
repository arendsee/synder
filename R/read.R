#' Read a synteny map
#'
#' @param file Synteny map file name
#' @return A dataframe
#' @export
read_synmap <- function(file) {
  d <- readr::read_tsv(file, col_types='ciiciidc', col_names=FALSE)

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
  d <- readr::read_tsv(file, col_types='ccciicccc', col_names=FALSE)
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
