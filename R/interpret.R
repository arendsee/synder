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
#' Creates a single row for each query, summarizing the set of search intervals
#' with the following descriptors:
#'
#' \enumerate{
#'   \item inbetween - no search interval overlaps a syntenic block
#'   \item lo_bound - at least one search interval has a lower overlapping edge
#'   \item hi_bound - at least one search interval has an upper overlapping edge
#'   \item doubly_bound - at least one search interval overlaps on both edges
#'   \item unbound - at least one search interval is unbound on both edges
#'   \item beyond - at least on search interval is on a target scaffold edge
#' }
#'
#' @param x SearchResult object
#' @export
flag_summary <- function(x){
  stopifnot(class(x) == "SearchResult")
  as.data.frame(x) %>%
    dplyr::group_by(.data$attr) %>%
    dplyr::summarize(
      inbetween    = all(.data$inbetween),
      lo_bound     = any(.data$l_flag < 2),
      hi_bound     = any(.data$r_flag < 2),
      doubly_bound = any(.data$l_flag < 2 & .data$r_flag < 2),
      unbound      = all(.data$l_flag > 1 & .data$r_flag > 1),
      beyond       = any(.data$l_flag == 4 | .data$r_flag == 4)
    )
}
