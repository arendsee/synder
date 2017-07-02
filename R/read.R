#' Synder read functions
#'
#' read synder classes from files
#'
#' @param file TAB-delimited file
#' @name synder_read
NULL


#' @rdname synder_read
#' @export
read_synmap <- function(file) {
  d <- readr::read_tsv(
    file,
    col_types = 'ciiciidc',
    col_names = FALSE,
    comment   = "#"
  )

  as_synmap(d)
}

#' @rdname synder_read
#' @export
read_gff <- function(file) {
  d <- readr::read_tsv(
    file,
    col_types = 'ccciincic',
    col_names = FALSE,
    na        = ".",
    comment   = "#"
  )

  as_gff(d)
}

#' @rdname synder_read
#' @export
read_hitmap <- function(file) {
  d <- readr::read_tsv(file, col_names=FALSE, comment="#")

  as_hitmap(d)
}

#' @rdname synder_read
#' @export
read_conlen <- function(file) {
  d <- readr::read_tsv(file, col_names=FALSE, comment="#")

  as_conlen(d)
}
