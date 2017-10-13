#' Find queries that map to the ends of a scaffold
#'
#' @param x SearchResult object
#' @export
find_truncations <- function(x){
  stopifnot(class(x) == "SearchResult")
  stop("This function is a stub")
}

#' Find queries that are in syntenically scambled locations
#'
#' @param x SearchResult object
#' @export
find_obfuscations <- function(x){
  stopifnot(class(x) == "SearchResult")
  stop("This function is a stub")
}

#' Calculate a metric for the overall density of the synteny map
#'
#' @param x Synmap object
#' @param qgen query-side scaffold lengths
#' @param tgen target-side scaffold lengths
#' @export
syntenic_density <- function(x, qgen=NULL, tget=NULL){
  stopifnot(class(x) == "Synmap")
  qgen <- as_conlen(qgen)
  tgen <- as_conlen(tgen)
  stop("This function is a stub")
}

#' Calculate a metric for the overall scrambledness of the synteny map
#'
#' @param x Synmap object
#' @param qgen query-side scaffold lengths
#' @param tgen target-side scaffold lengths
#' @export
syntenic_scatter <- function(x, qgen=NULL, tget=NULL){
  stopifnot(class(x) == "SearchResult")
  qgen <- as_conlen(qgen)
  tgen <- as_conlen(tgen)
  stop("This function is a stub")
}

#' Calculate the local scrambledness in the region around each query
#'
#' @param x SearchResult object
#' @param synmap Synmap object
#' @export
neighborhood_madness <- function(x, synmap){
  stopifnot(class(x) == "SearchResult")
  stopifnot(class(synmap) == "Synmap")
  stop("This function is a stub")
}

#' Summarize the flags returned from a search
#'  
#' @param x SearchResult object
#' @export
flag_summary <- function(x){
  stopifnot(class(x) == "SearchResult")
  stop("This function is a stub")
}
