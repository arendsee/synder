#' Synder plot functions
#'
#' plots synder classes
#'
#' @param x table of a synder class
#' @param y synmap object (for context)
#' @name synder_plot
NULL

# ----------------------
# Plot utility functions
# ----------------------

to_global <- function(x, prefix='q') {
  xcon <- paste0(prefix, 'con')
  xstop  <- paste0(prefix, 'stop')
  xstart <- paste0(prefix, 'start')
  xv <- dplyr::group_by_(x, xcon)   %>%
    dplyr::mutate_(xmax=max(xstop)) %>%
    dplyr::select_(xcon, 'xmax')    %>%
    unique                          %>%
    dplyr::arrange_('xmax')         %>%
    {
      xlv <- c(0, .$xmax[-nrow(.)])
      names(xlv) <- .[[xcon]]
      xlv
    }
  x[xstart] = x[[xstart]] + xv[x[[xcon]]]
  x[xstop]  = x[[xstop]]  + xv[x[[xcon]]]
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

# ----------------------
# Plot functions
# ----------------------

#' @rdname synder_plot
#' @export
plot.synmap <- function(x){
  x <- to_global(x, prefix='q')
  x <- to_global(x, prefix='t')
  x <- zero_base(x)
  ggplot2::ggplot(x) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x=qstart, xend=qstop,
        y=tstart, yend=tstop
      )
    )
}

#' @rdname synder_plot
#' @export
plot.gff <- function(x){ }

#' @rdname synder_plot
#' @export
plot.hitmap <- function(x){ }

#' @rdname synder_plot
#' @export
plot.dump_result <- function(x){
  x <- to_global(x, prefix='q')
  x <- to_global(x, prefix='t')
  x <- zero_base(x)

  x$cset <- factor(as.integer(x$cset))
  csets <- dplyr::group_by(x, cset) %>%
    dplyr::summarize(
      qcon   = qcon[1],
      qstart = min(qstart),
      qstop  = max(qstop),
      tcon   = tcon[1],
      tstart = min(tstart),
      tstop  = max(tstop),
      strand = strand[1]
    )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data=x,
      ggplot2::aes(
        x=qstart, xend=qstop,
        y=tstart, yend=tstop
      )
    ) +
    ggplot2::geom_segment(
      data=csets,
      ggplot2::aes(
        x=qstart, xend=qstop,
        y=tstart, yend=tstop
      ),
      alpha=0.3,
      size=3
    )
}

#' @rdname synder_plot
#' @export
plot.search_result <- function(x, y){
  stopifnot('synmap' %in% class(y))
  x <- x[2:8]
  x$group <- 'search'
  y$group <- 'synmap'
  y$score <- NULL
  z <- rbind(x, y)

  z <- to_global(z, prefix='q')
  z <- to_global(z, prefix='t')
  z <- zero_base(z)

  y <- subset(z, group == 'synmap')
  x <- subset(z, group == 'search')

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data=y,
      ggplot2::aes(
        x=qstart, xend=qstop,
        y=tstart, yend=tstop
      )
    ) +
    ggplot2::geom_rect(
      data=x,
      ggplot2::aes(
        xmin=qstart, xmax=qstop,
        ymin=tstart, ymax=tstop
      ),
      color='red',
      alpha=0.2
    )
}

#' @rdname synder_plot
#' @export
plot.filter_result <- function(x, y){ }

#' @rdname synder_plot
#' @export
plot.map_result <- function(x, y){ }

#' @rdname synder_plot
#' @export
plot.count_result <- function(x, y){ }
