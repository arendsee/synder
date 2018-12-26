.namedlist <- function(...){
  ns <- unlist(lapply(as.list(substitute(list(...)))[-1], deparse))
  xs <- list(...)
  names(xs) <- ns
  xs
}

.as_grange <- function(x, xid, xa, xb){
  GenomicRanges::GRanges(
    seqnames=x[[xid]],
    ranges=IRanges::IRanges(
      start = pmin(x[[xa]], x[[xb]]),
      end = pmax(x[[xa]], x[[xb]])
    )
  )
}

.remove_rownames <- function(d){
  rownames(d) <- NULL
  d
}

.as_map <- function(d, xkey=1, xval=2){
  stopifnot(ncol(d) >= max(xkey, xval))
  d <- d[, c(xkey, xval)]
  d <- dplyr::distinct(d)
  x <- d[[2]]
  names(x) <- d[[1]]
  x
}

# make column `n` the first column, shifting the other columns up
.as_first_column <- function(d, n){
  i <- which(names(d) == n) 
  d[, c(i, (1:ncol(d))[-i])]
}

#' Join two tables by overlapping rows
#'
#' @param x data.frame: table with interval info in each row (with arbitrary
#' additional columns).
#' @param xid character: name of the scaffold column in x
#' @param xa character: name of the start position column in x
#' @param xb character: name of the stop position column in x
#' @param yid character: name of the scaffold column in y
#' @param ya character: name of the start position column in y
#' @param yb character: name of the stop position column in y
#' @param add_id logical: whether to add a column of indices holding the column
#' orders of the original tables (named "xid" and "yid", respectively)
#' @return data.frame with columns corresponding to those of x and y and rows
#' being cases where intervals overlap between the two.
.overlaps <- function(x, xid, xa, xb, y, yid=xid, ya=xa, yb=xb, add_id=TRUE){
  # add unique ID to each row in each input table
  if(add_id){
    x$xid <- 1:nrow(x)
    y$yid <- 1:nrow(y)
  }
  # convert each input table to a GRanges object
  xrng <- .as_grange(x, xid, xa, xb)
  GenomicRanges::mcols(xrng) <- x[, which(!(names(x) %in% c(xid, xa, xb)))]
  yrng <- .as_grange(y, yid, ya, yb)
  GenomicRanges::mcols(yrng) <- y[, which(!(names(y) %in% c(yid, ya, yb)))]
  # find the interval overlaps
  hits <- GenomicRanges::findOverlaps(xrng, yrng, ignore.strand=TRUE) 
  # extract a table of overlapping rows from the GRanges object
  dodt <- function(z, ids, id, a, b){
    z <- as.data.frame(z)[ids, ]
    z$strand <- NULL
    z$width <- NULL
    names(z)[1:3] <- c(id, a, b)
    for(n in c("seqname", "start", "stop", "width", "strand")){
      if(any(paste0(n, ".1") %in% names(z))){
        names(z)[names(z) %in% paste0(n, ".1")] <- n
      }
    }
    z
  }
  x2 <- dodt(xrng, S4Vectors::queryHits(hits), xid, xa, xb)
  y2 <- dodt(yrng, S4Vectors::subjectHits(hits), yid, ya, yb)
  # since yid == xid for all overlaps, remove yid (avoid duplicated columns)
  y2[[yid]] <- NULL
  # bind the two tables together
  cbind(x2, y2)
}
