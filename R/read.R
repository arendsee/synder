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
  readr::read_tsv(
    file,
    col_types = 'ciiciidc',
    col_names = names(SYNMAP_COLS),
    comment   = "#"
  ) %>% as_synmap
}

#' @rdname synder_read
#' @export
read_gff <- function(file) {
  readr::read_tsv(
    file,
    col_types = 'ccciincic',
    col_names = names(GFF_COLS),
    na        = ".",
    comment   = "#"
  ) %>% as_gff
}

#' @rdname synder_read
#' @export
read_hitmap <- function(file) {
  readr::read_tsv(
    file,
    col_names = FALSE,
    comment   = "#"
  ) %>% as_hitmap
}

#' @rdname synder_read
#' @export
read_conlen <- function(file) {
  readr::read_tsv(
    file,
    col_names=names(CON_LENGTH, comment="#")
  ) %>% as_conlen
}
