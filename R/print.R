#' Synder print functions
#'
#' Prints the special attributes of the class and then prints a parent class.
#'
#' @param x table of a synder class
#' @param ... additional arguments that are currently ignored
#' @name synder_print
NULL

simple_print_wrapper <- function(x, ...){
  cat(sprintf("# %s\n", class(x)[1]))
  class(x) <- class(x)[-1]
  print(x)
}

paramaterized_print_wrapper <- function(x, args){
  arg2str <- function(p) {
    sprintf("%s=%s", p, eval(parse(text=sprintf("attributes(x)$%s", p))))
  }
  argstr <- lapply(args, arg2str) %>% unlist %>% paste0(collapse=", ")
  sprintf("# %s -- %s\n", class(x)[1], argstr) %>% cat
  class(x) <- class(x)[-1]
  print(x)
}

#' @rdname synder_print
#' @export
print.synmap <- function(x, ...){
  simple_print_wrapper(x)
}

#' @rdname synder_print
#' @export
print.gff <- function(x, ...){
  simple_print_wrapper(x)
}

#' @rdname synder_print
#' @export
print.hitmap <- function(x, ...){
  simple_print_wrapper(x)
}

#' @rdname synder_print
#' @export
print.dump_result <- function(x, ...){
  paramaterized_print_wrapper(x, c('strans', 'swap'))
}

#' @rdname synder_print
#' @export
print.search_result <- function(x, ...){
  paramaterized_print_wrapper(x, c('swap', 'trans', 'k', 'r', 'offsets'))
}

#' @rdname synder_print
#' @export
print.filter_result <- function(x, ...){
  paramaterized_print_wrapper(x, c('swap', 'trans', 'k', 'r', 'offsets'))
}

#' @rdname synder_print
#' @export
print.map_result <- function(x, ...){
  paramaterized_print_wrapper(x, c('swap', 'offsets'))
}

#' @rdname synder_print
#' @export
print.count_result <- function(x, ...){
  paramaterized_print_wrapper(x, c('swap', 'offsets'))
}
