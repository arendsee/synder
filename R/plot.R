# ----------------------
# Plot utility functions
# ----------------------

to_global <- function(x, prefix='q') {
  xseqid <- paste0(prefix, 'seqid')
  xstop  <- paste0(prefix, 'stop')
  xstart <- paste0(prefix, 'start')
  xv <- dplyr::group_by(x, .data[[xseqid]]) %>%
    dplyr::mutate(xmax=max(.data[[xstop]])) %>%
    dplyr::select_(xseqid, 'xmax')                 %>%
    base::unique()                                 %>%
    dplyr::arrange(.data$xmax)              %>%
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

#' Plot a synteny map
#'
#' @method plot Synmap
#' @export
#' @param x Synmap object
#' @param ... Additional arguments (not currently used)
plot.Synmap <- function(x, ...){

  x <- as.data.frame(x)
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

#' Print the result of dumping a synteny map
#'
#' @method plot DumpResult
#' @export
#' @param x DumpResult object
#' @param ... Additional arguments (not currently used)
plot.DumpResult <- function(x, ...){

  x <- as.data.frame(x)
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

#' Plot the result of a search operations
#'
#' This plot function takes both a SearchResult object and the Synmap object
#' that was used to produce it. The latter is required to plot the search
#' intervals in the synmap context.
#'
#' @method plot SearchResult
#' @export
#' @param x SearchResult object
#' @param y Synmap object
#' @param rank logical Should intervals be printed by rank?
#' @param ... Additional arguments (not currently used)
plot.SearchResult <- function(x, y, rank=FALSE, ...){

  stopifnot(is_synmap(y))

  x <- as.data.frame(x)
  y <- as.data.frame(y)

  x <- x[2:8]
  x$group <- 'search'
  y$group <- 'synmap'
  y$score <- NULL
  z <- rbind(x, y)

  if(rank){
    zq <- z[,1:3] %>%
      reshape2::melt() %>%
      dplyr::group_by(.data$qseqid) %>%
      dplyr::mutate(value = as.numeric(ordered(.data$value))) %>% {
        data.frame(
          qseqid = .$qseqid[.$variable == "qstart"],
          qstart = .$value[.$variable == "qstart"],
          qstop  = .$value[.$variable == "qstop"]
        )
      }

    zt <- z[,4:6] %>%
      reshape2::melt() %>%
      dplyr::group_by(.data$tseqid) %>%
      dplyr::mutate(value = as.numeric(ordered(.data$value))) %>% {
        data.frame(
          tseqid = .$tseqid[.$variable == "tstart"],
          tstart = .$value[.$variable == "tstart"],
          tstop  = .$value[.$variable == "tstop"]
        )
      }

    z <- cbind(zq, zt, z[,7:8])
  }

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
