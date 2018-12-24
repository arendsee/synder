#' Map features between GFFs based on a synder SearchResult
#'
#' @param srcres SearchResult object from a synder::search call
#' @param tgff GFF object for target species
#' @param fgff GFF object for focal species
#' @export
featureMap <- function(srcres, tgff, fgff){

  hits <- GenomicRanges::findOverlaps(CNEr::second(srcres), tgff, ignore.strand=TRUE)
  s <- as.data.frame(srcres)[S4Vectors::queryHits(hits), ] %>%
    dplyr::rename(
      ffeature = attr,
      fchr = qseqid,
      fstart = qstart,
      fstop = qstop,
      si_start = tstart,
      si_stop = tstop,
      tchr = tseqid,
      orientation = strand,
      synder_score = score
    )
  # hit order
  s$order <- 1:nrow(s)

  tg <- as.data.frame(tgff)[S4Vectors::subjectHits(hits), ] %>%
    dplyr::select(
      tfeature = attr,
      ttype = type,
      tsource = source,
      tchr = seqid,
      tstart = start,
      tstop = stop,
      tstrand = strand
    )
  # hit order
  tg$order <- 1:nrow(tg)

  fg <- as.data.frame(fgff)[match(s$ffeature, fgff$attr), ] %>%
    dplyr::select(
      ffeature = attr,
      ftype = type,
      fsource = source,
      fchr = seqid,
      fstart = start,
      fstop = stop,
      fstrand = strand
    )

  # the search result queries are the focal features
  stopifnot(fg$ffeature == s$ffeature)
  stopifnot(fg$fchr == s$fchr)
  stopifnot(fg$fstart == s$fstart)
  stopifnot(fg$fstop == s$fstop)
  # the search intervals are match the target GFF chromosome/scaffold names
  stopifnot(tg$tchr == s$tchr)
  # all have the same number of rows
  stopifnot(nrow(tg) == nrow(fg))
  stopifnot(nrow(tg) == nrow(s))

  d1 <- merge(s, dplyr::distinct(fg[, -c(4,5,6)]), by="ffeature")
  stopifnot(nrow(d1) == nrow(s))

  d2 <- merge(d1, tg, by="order")
  stopifnot(nrow(d2) == nrow(s))
  stopifnot(d2$tchr.x == d2$tchr.y)

  d2$strands_agree <- with(d2, 
    # if the syntenic orientation is the same, then they will be on the same strand 
    (orientation == "+" & (fstrand == tstrand)) |
    # otherwise they will be on opposite strands
    (orientation == "-" & (fstrand != tstrand))
  )

  d3 <- dplyr::rename(d2, tchr=tchr.x) %>% dplyr::select(
    ## focal feature locations
    ffeature, ftype, fsource, fchr, fstart, fstop, fstrand,
    ## target feature info
    tfeature, ttype, tsource, tchr, tstart, tstop, tstrand,
    ## synder search interval info
    si_start, si_stop, orientation, synder_score, cset, l_flag, r_flag, inbetween,
    ## derived data
    strands_agree
  )
}
