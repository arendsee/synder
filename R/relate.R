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
      fseqid = attr,
      fchr = qseqid,
      fstart = qstart,
      fstop = qstop,
      si_start = tstart,
      si_stop = tstop,
      tchr = tseqid,
      orientation = strand,
      si_score = score
    )
  # hit order
  s$order <- 1:nrow(s)

  tg <- as.data.frame(tgff)[S4Vectors::subjectHits(hits), ] %>%
    dplyr::select(
      tseqid = attr,
      ttype = type,
      tsource = source,
      tchr = seqid,
      tstart = start,
      tstop = stop,
      tstrand = strand
    )
  # hit order
  tg$order <- 1:nrow(tg)

  fg <- as.data.frame(fgff)[match(s$fseqid, fgff$attr), ] %>%
    dplyr::select(
      fseqid = attr,
      ftype = type,
      fsource = source,
      fchr = seqid,
      fstart = start,
      fstop = stop,
      fstrand = strand
    )

  # the search result queries are the focal features
  stopifnot(fg$fseqid == s$fseqid)
  stopifnot(fg$fchr == s$fchr)
  stopifnot(fg$fstart == s$fstart)
  stopifnot(fg$fstop == s$fstop)
  # the search intervals are match the target GFF chromosome/scaffold names
  stopifnot(tg$tchr == s$tchr)
  # all have the same number of rows
  stopifnot(nrow(tg) == nrow(fg))
  stopifnot(nrow(tg) == nrow(s))

  d1 <- merge(s, dplyr::distinct(fg[, -c(4,5,6)]), by="fseqid")
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
    fseqid, ftype, fsource, fchr, fstart, fstop, fstrand,
    ## target feature info
    tseqid, ttype, tsource, tchr, tstart, tstop, tstrand,
    ## synder search interval info
    si_start, si_stop, orientation, si_score, cset, l_flag, r_flag, inbetween,
    ## derived data
    strands_agree
  )
}

#' make gene map
#'
#' @param tblastn data.frame
#' @param tgff data.frame
#' @export
make_tblastn_gene_map <- function(tblastn, tgff){
  .overlaps(
    x=tblastn, y=as.data.frame(tgff),
    xid="sseqid", xa="sstart", xb="send",
    yid="seqid", ya="start", yb="stop"
  ) %>%
    dplyr::select(
      ## the focal protein that was BLASTed
      fprot_id = qseqid,
      fprot_start = qstart,
      fprot_stop = qend,
      ## the target genomic intervals with translated coding similarity
      tchr = sseqid,
      tn_start = sstart,
      tn_stop = send,
      tn_frame = sframe,
      tn_evalue = evalue,
      tn_ppos = ppos,
      ## the target feature that overlaps the region of similarity
      tfeat_name = attr,
      tfeat_start = start,
      tfeat_stop = stop,
      tfeat_src = source,
      tfeat_type = type,
      tfeat_strand = strand
    ) %>% .remove_rownames
}

#' tBLASTn search interval map
#'
#' @param tblastn data.frame
#' @param synres data.frame
#' @export
make_tblastn_si_map <- function(tblastn, synres){
  x <- dplyr::select(tblastn,
    # the focal protein that was BLASTed
    fprot_id = qseqid,
    fprot_start = qstart,
    fprot_stop = qend,
    # the target genomic intervals with translated coding similarity
    tchr = sseqid,
    tn_start = sstart,
    tn_stop = send,
    tn_frame = sframe,
    tn_evalue = evalue,
    tn_ppos = ppos
  )
  y <- dplyr::select(as.data.frame(synres),
    # synder search interval info
    ## * focal query interval
    fseqid = attr,
    fchr = qseqid,
    si_fstart = qstart,
    si_fstop = qstop,
    ## * target search interval
    si_tchr = tseqid,
    si_tstart = tstart,
    si_tstop = tstop,
    orientation = strand,
    ## * synder stats
    si_score = score,
    cset = cset,
    l_flag = l_flag,
    r_flag = r_flag,
    inbetween = inbetween
  )
  .overlaps(x, y,
           xid="tchr", xa="tn_start", xb="tn_stop",
           yid="si_tchr", ya="si_tstart", yb="si_tstop", add_id=FALSE) %>%
    dplyr::select(-si_tchr) %>%
    dplyr::select(
      # focal protein
      fprot_id, fprot_start, fprot_stop,
      # region of coding similarity found by tBLASTn 
      tchr, tn_start, tn_stop, tn_frame,
      # tBLASTn stats
      tn_evalue, tn_ppos,
      # synder search interval
      fseqid, cset,
      # - focal region
      fchr, si_fstart, si_fstop,
      # - target region
      si_fstart, si_fstop, orientation,
      # - synder stats
      si_score, l_flag, r_flag, inbetween
    ) %>% .remove_rownames
}

make_blastp_si_map <- function(blastp, srcres, fgff, tgff){
  NULL
}

make_blastp_gene_map <- function(blastp, srcres, fgff, tgff){
  NULL
}
