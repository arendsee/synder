#' Synder print functions
#'
#' Prints the special attributes of the class and then prints a parent class.
#'
#' @param x table of a synder class
#' @param ... additional arguments that are currently ignored
#' @name synder_print
NULL

#' @rdname synder_print
#' @export
print.Synmap <- function(x, ...){
  NextMethod(x, ...)
}

#' @rdname synder_print
#' @export
print.GFF <- function(x, ...){
  NextMethod(x, ...)
}

#' @rdname synder_print
#' @export
print.DumpResult <- function(x, ...){
  NextMethod(x)
  cat(sprintf("swap=%s  trans=%s  offsets=%s\n",
    x@swap,
    x@trans,
    paste(x@offsets, collapse="")
  ))
}

#' @rdname synder_print
#' @export
print.SearchResult <- function(x, ...){
  NextMethod(x)
  cat(sprintf("swap=%s  trans=%s  k=%s  r=%s  offsets=%s\n",
    x@swap,
    x@trans,
    x@k,
    x@r,
    paste(x@offsets, collapse="")
  ))
}
