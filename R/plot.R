#' Synder plot functions
#'
#' plots synder classes
#'
#' @param x table of a synder class
#' @param y synmap object (for context)
#' @param ... additional arguments that are currently ignored 
#' @name synder_plot
NULL

# ----------------------
# Plot utility functions
# ----------------------

to_global <- function(x, prefix='q') {
  xseqid <- paste0(prefix, 'seqid')
  xstop  <- paste0(prefix, 'stop')
  xstart <- paste0(prefix, 'start')
  xv <- dplyr::group_by_(x, xseqid) %>%
    dplyr::mutate_(xmax=max(xstop)) %>%
    dplyr::select_(xseqid, 'xmax')  %>%
    unique                          %>%
    dplyr::arrange_('xmax')         %>%
    {
      xlv <- c(0, .$xmax[-nrow(.)])
      names(xlv) <- .[[xseqid]]
      xlv
    }
  x[xstart] = x[[xstart]] + xv[x[[xseqid]]]
  x[xstop]  = x[[xstop]]  + xv[x[[xseqid]]]
  x
}

zero_base <- function(x) {
  qmin <- min(x$qstart)
  tmin <- min(x$tstart)

  x$qstart <- x$qstart - qmin
  x$tstart <- x$tstart - tmin
  x$qstop  <- x$qstop  - qmin
  x$tstop  <- x$tstop  - tmin

  x
}

handle_inversions <- function(x) {
  tstart <- ifelse(x$strand == '+', x$tstart, x$tstop)
  tstop <- ifelse(x$strand == '+', x$tstop, x$tstart)
  x$tstart <- tstart
  x$tstop <- tstop
  x
}

# ----------------------
# Plot functions
# ----------------------

#' @rdname synder_plot
#' @method plot synmap
#' @export
plot.synmap <- function(x, ...){
  x <- to_global(x, prefix='q')
  x <- to_global(x, prefix='t')
  x <- zero_base(x)
  x <- handle_inversions(x)

  ggplot2::ggplot(x) +
    ggplot2::geom_segment(
      ggplot2::aes_string(
        x="qstart", xend="qstop",
        y="tstart", yend="tstop"
      )
    )
}

#' @rdname synder_plot
#' @method plot gff
#' @export
plot.gff <- function(x, ...){ }

#' @rdname synder_plot
#' @method plot hitmap
#' @export
plot.hitmap <- function(x, ...){ }

#' @rdname synder_plot
#' @method plot dump_result
#' @export
plot.dump_result <- function(x, ...){
  x <- to_global(x, prefix='q')
  x <- to_global(x, prefix='t')
  x <- zero_base(x)

  x$cset <- factor(as.integer(x$cset))
  csets <- dplyr::group_by_(x, "cset") %>%
    dplyr::summarize_(
      qseqid = "qseqid[1]",
      qstart = "min(qstart)",
      qstop  = "max(qstop)",
      tseqid = "tseqid[1]",
      tstart = "min(tstart)",
      tstop  = "max(tstop)",
      strand = "strand[1]"
    )

  x     <- handle_inversions(x)
  csets <- handle_inversions(csets)

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data=x,
      ggplot2::aes_string(
        x="qstart", xend="qstop",
        y="tstart", yend="tstop"
      )
    ) +
    ggplot2::geom_segment(
      data=csets,
      ggplot2::aes_string(
        x="qstart", xend="qstop",
        y="tstart", yend="tstop"
      ),
      alpha=0.3,
      size=3
    )
}

#' @rdname synder_plot
#' @method plot search_result
#' @export
plot.search_result <- function(x, y, ...){
  stopifnot('synmap' %in% class(y))
  x <- x[2:8]
  x$group <- 'search'
  y$group <- 'synmap'
  y$score <- NULL
  z <- rbind(x, y)

  z <- to_global(z, prefix='q')
  z <- to_global(z, prefix='t')
  z <- zero_base(z)

  y <- dplyr::filter_(z, "group == 'synmap'")
  x <- dplyr::filter_(z, "group == 'search'")

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data=y,
      ggplot2::aes_string(
        x="qstart", xend="qstop",
        y="tstart", yend="tstop"
      )
    ) +
    ggplot2::geom_rect(
      data=x,
      ggplot2::aes_string(
        xmin="qstart", xmax="qstop",
        ymin="tstart", ymax="tstop"
      ),
      color='red',
      alpha=0.2
    )
}

#' @rdname synder_plot
#' @method plot filter_result
#' @export
# TODO: complete
plot.filter_result <- function(x, y, ...){ }

#' @rdname synder_plot
#' @method plot map_result 
#' @export
# TODO: complete
plot.map_result <- function(x, y, ...){ }

#' @rdname synder_plot
#' @method plot count_result
#' @export
# TODO: complete
plot.count_result <- function(x, y, ...){ }
