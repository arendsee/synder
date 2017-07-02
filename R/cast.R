#' Synder cast functions
#'
#' cast data as Synder objects
#'
#' @param d input type
#' @name synder_cast
NULL


#' @rdname synder_cast
#' @export
as_synmap <- function(d) {
  d <- as.data.frame(d)

  # assert correct number of columns
  stopifnot(ncol(d) == 8)

  # Set column names
  names(d) <- names(SYNMAP_COLS)

  # Assert column types match expectations
  stopifnot(SYNMAP_COLS == lapply(d, class))

  # Assign class membership
  class(d) <- append('synmap', class(d))

  return(d)
}

#' @rdname synder_cast
#' @export
as_gff <- function(d) {
  d <- as.data.frame(d)

  # Assert the correct number of columns were read
  stopifnot(ncol(d) == 9)

  # Set column names
  names(d) <- names(GFF_COLS)

  # Assert column types match expectations
  for(i in 1:ncol(d)){
    if(class(d[[i]]) != GFF_COLS[i]){
      stop(paste("GFFError: in column", i, "expected type", GFF_COLS[i],
                 ", found", class(d[[i]])))
    }
  }

  # Assign class membership
  class(d) <- append('gff', class(d))

  return(d)
}

#' @rdname synder_cast
#' @export
as_hitmap <- function(d) {
  d <- as.data.frame(d)

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

#' @rdname synder_cast
#' @export
as_conlen <- function(d) {
  d <- as.data.frame(d)

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
