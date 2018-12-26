.tidy_gff <- function(gff, prefix="", fromAttr=identity, rm_col=c("score","phase","attr")){
  gff <- as.data.frame(synder::as_gff(gff))
  gff <- dplyr::rename(gff, chr=seqid)
  gff$seqid <- fromAttr(gff$attr)
  gff <- .as_first_column(gff, "seqid")
  gff <- gff[, -which(names(gff) %in% rm_col)]
  names(gff) <- paste0(prefix, names(gff))
  gff
}

.tidy_searchResult <- function(srcres, prefix="si_", fromAttr=identity){
  s <- as.data.frame(srcres)
  s$fseqid = fromAttr(s$attr)
  s <- .as_first_column(s, "fseqid")
  s$attr <- NULL
  s <- dplyr::rename(s,
    fchr = qseqid,
    fstart = qstart,
    fstop = qstop,
    tchr = tseqid,
    tstart = tstart,
    tstop = tstop,
    orientation = strand
  )
  names(s)[c(3,4,6,7,9)] <- paste0(prefix, names(s)[c(3,4,6,7,9)])
  s
}

#' Map features between GFFs based on a synder SearchResult
#'
#' @param srcres tidy search result (from .tidy_searchResult)
#' @param tgff data.frame for GFF object for target species (from .tidy_gff)
#' @param fgff data.frame for GFF object for focal species (from .tidy_gff)
#' @export
featureMap <- function(srcres, tgff, fgff){

  s <- .overlaps(
    x=srcres, xid="tchr", xa="si_tstart", xb="si_tstop",
    y=tgff,   yid="tchr", ya="tstart",    yb="tstop"
  )
  s <- .as_first_column(s, "fseqid")

  s <- merge(s, fgff, by='fseqid')
  stopifnot(s$fchr.x == s$fchr.y)

  s$strands_agree <- with(s, 
    # if the syntenic orientation is the same, then they will be on the same strand 
    (orientation == "+" & (fstrand == tstrand)) |
    # otherwise they will be on opposite strands
    (orientation == "-" & (fstrand != tstrand))
  )

  dplyr::select(s,
    ## focal feature locations
    fseqid, ftype, fsource, fchr=fchr.x, fstart, fstop, fstrand,
    ## target feature info
    tseqid, ttype, tsource, tchr, tstart, tstop, tstrand,
    ## synder search interval info
    si_fstart, si_fstop, si_tstart, si_tstop, orientation,
    ## synder stats
    si_score, cset, l_flag, r_flag, inbetween,
    ## derived data
    strands_agree
  )
}

.read.blastfile <- function(filename, qtype="fseqid", ttype="tchr", prefix="tbn_"){
  out <- read.table(
    filename,
    header=TRUE,
    stringsAsFactors=FALSE
  )
  names(out)[which(names(out) == "sseqid")] <- ttype
  names(out)[which(names(out) == "qseqid")] <- qtype
  out <- .as_first_column(out, ttype)
  out <- .as_first_column(out, qtype)
  out$id <- 1:nrow(out)
  i <- which(!(names(out) %in% c(ttype, qtype)))
  names(out)[i] <- paste0(prefix, names(out)[i])
  out
}

#' make gene map
#'
#' @param tblastn_file filename for tblastn results
#' @param tgff a tidy GFF file (tidied with .tidy_gff)
#' @param fromAttr a function for extracting the sequence ID from the GFF attribute column
#' @export
make_tblastn_gene_map <- function(tblastn_file, tgff, fmap=NULL){
  tblastn <- .load_tblastn_file(tblastn_file, fmap=fmap)
  out <- .overlaps(
    x=tblastn, y=tgff,
    xid="tchr", xa="tn_sstart", xb="tn_send",
    yid="tchr", ya="tstart", yb="tstop"
  ) %>%
    dplyr::select(
      ## the focal protein that was BLASTed
      fseqid, fprotein,
      prot_fstart = tn_qstart,
      prot_fstop = tn_qend,
      ## the target genomic intervals with translated coding similarity
      tchr,
      tn_tstart = tn_sstart,
      tn_tstop = tn_send,
      tn_tframe = tn_sframe,
      tn_evalue,
      tn_ppos,
      ## the target feature that overlaps the region of similarity
      feat_tname = tseqid,
      feat_tstart = tstart,
      feat_tstop = tstop,
      feat_tsrc = tsource,
      feat_ttype = ttype,
      feat_tstrand = tstrand
    ) %>% .remove_rownames
  list(blastmap=out, blastraw=tblastn)
}

.load_tblastn_file <- function(tblastn_file, fmap=NULL){
  tblastn <- .read.blastfile(tblastn_file, qtype="fprotein", ttype="tchr", prefix="tn_")
  if(!is.null(fmap)){
    tblastn <- merge(dplyr::distinct(fmap), tblastn, by='fprotein')
  } else {
    tblastn$fseqid <- tblastn$fprotein
  }
  tblastn <- .as_first_column(tblastn, "fprotein")
  tblastn <- .as_first_column(tblastn, "fseqid")
  tblastn
}

#' tBLASTn search interval map
#'
#' @param tblastn_file filename of tblastn result
#' @param synres data.frame
#' @export
make_tblastn_si_map <- function(tblastn_file, srcres, fmap=NULL){
  tblastn <- .load_tblastn_file(tblastn_file, fmap=fmap)
  x <- dplyr::rename(tblastn,
    # the focal protein that was BLASTed
    prot_fstart = tn_qstart,
    prot_fstop = tn_qend,
    # the target genomic intervals with translated coding similarity
    tn_tstart = tn_sstart,
    tn_tstop = tn_send
  )
  srcres <- dplyr::rename(srcres, si_fseqid=fseqid)
  out <- .overlaps(x, srcres,
           xid="tchr", xa="tn_tstart",  xb="tn_tstop",
           yid="tchr", ya="si_tstart", yb="si_tstop", add_id=FALSE) %>%
    dplyr::select(
      # focal protein
      fseqid, fprotein, prot_fstart, prot_fstop,
      # region of coding similarity found by tBLASTn 
      tchr, tn_tstart, tn_tstop, tn_tframe=tn_sframe,
      # tBLASTn stats
      tn_evalue, tn_ppos,
      # synder search interval
      si_fseqid, cset,
      # - focal region
      fchr, si_fstart, si_fstop,
      # - target region
      si_fstart, si_fstop, orientation,
      # - synder stats
      si_score, l_flag, r_flag, inbetween
    ) %>% .remove_rownames
  list(blastmap=out, blastraw=tblastn)
}

#' Map intervals for protein blast
#'
#' Returns info on the focal and target intervals where 1) there is a
#' significant protein BLAST hit to the proteins encoded by the focal and
#' target features and 2) the hit overlaps a search interval on the target
#' side. The columns \code{fseqid} and \code{si_fseqid} hold the protein
#' encoding focal feature and the focal interval used in the synder search,
#' respectively. These two will not always be the same.
#'
#' @param blastp_file filename: with columns {qseqid, sseqid, qlen, slen, length, evalue, ppos}
#' @param fgff GFF: focal species (from \code{.tidy_gff})
#' @param tgff GFF: target species (from \code{.tidy_gff})
#' @param srcres SearchResult object (from \code{.tidy_searchResult})
#' @param fmap data.frame: a map from blastp qseqid (a protein ID) to the seqid
#' in the SearchResult (probably a gene ID for an interval)
#' @param tmap data.frame: like fmap but for the target proteins
#' @return list containing the blast results filtered through synder and the raw blast results
#' @export
make_blastp_map <- function(blastp_file, fgff, tgff, srcres, fmap=NULL, tmap=NULL){
  d1 <- .read.blastfile(blastp_file, qtype="fprotein", ttype="tprotein", prefix="bp_")

  # map the target protein ID to the feature ID used in the target GFF 
  if(is.null(tmap)){
    d1$tseqid <- d1$tprotein
  } else {
    stopifnot(d1$tprotein %in% tmap[[1]])
    d1$tseqid <- .as_map(tmap)[d1$tprotein]
  }
  d1 <- .as_first_column(d1, "tseqid")

  # map the focal protein ID to the feature ID used in the query GFF 
  if(is.null(fmap)){
    d1$fseqid <- d1$fprotein
  } else {
    stopifnot(d1$fprotein %in% fmap[[1]])
    d1$fseqid <- .as_map(fmap)[d1$fprotein]
  }
  d1 <- .as_first_column(d1, "fseqid")

  d2 <- d1 %>%
    merge(by="fseqid", fgff) %>%
    merge(by="tseqid", tgff)

  y <- dplyr::rename(srcres, si_fseqid = fseqid, si_fchr=fchr)
  out <- .overlaps(
    x=d2, xid="tchr", xa="tstart", xb="tstop",               
    y=y,  yid="tchr", ya="si_tstart", yb="si_tstop" 
  ) %>%
    dplyr::select(
      # key ids
      fseqid, si_fseqid, tseqid, fprotein, tprotein,
      # focal feature interval
      fchr, fsource, ftype, fstart, fstop, fstrand,
      # target feature interval
      tchr, tsource, ttype, tstart, tstop, tstrand,
      # blastp stats
      bp_qlen, bp_slen, bp_length, bp_evalue, bp_ppos,
      # search interval
      si_fchr, si_fstart, si_fstop,
      si_tstart, si_tstop, orientation,
      # search stats
      cset, si_score, l_flag, r_flag, inbetween
    )
  list(blastmap=out, blastraw=d1)
}
