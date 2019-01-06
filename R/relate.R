.tidy_gff <- function(gff, prefix="", fromAttr=identity, rm_col=c("score","phase","attr")){
  if(! is.data.frame(gff)){
    gff <- as.data.frame(as_gff(gff))
    gff <- dplyr::rename(gff, chr=seqid)
    gff$seqid <- fromAttr(gff$attr)
    gff <- .as_first_column(gff, "seqid")
    gff <- gff[, -which(names(gff) %in% rm_col)]
    names(gff) <- paste0(prefix, names(gff))
  }
  gff
}

.tidy_searchResult <- function(srcres, prefix="si_", fromAttr=identity){
  if(! is.data.frame(srcres)){
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
  } else {
    rcol <- c("fseqid", "fchr", "si_fstart", "si_fstop", "tchr", "si_tstart",
              "si_tstop", "orientation", "si_score", "cset", "l_flag",
              "r_flag", "inbetween")
    if(! all(names(srcres) == rcol)){
      stop("In search interval, missing expected columns: ", setdiff(rcol, names(srcres)))
    }
    srcres 
  }
}

#' Map features between GFFs based on a synder SearchResult
#'
#' @param srcres tidy search result (from .tidy_searchResult)
#' @param tgff data.frame for GFF object for target species (from .tidy_gff)
#' @param fgff data.frame for GFF object for focal species (from .tidy_gff)
#' @param tfromAttr function (optional) extract IDs attr column in target GFF
#' @param ffromAttr function (optional) extract IDs attr column in focal GFF
#' @export
featureMap <- function(srcres, tgff, fgff, tfromAttr=identity, ffromAttr=identity){
  srcres <- .tidy_searchResult(srcres)
  fgff <- .tidy_gff(fgff, prefix="f", fromAttr=ffromAttr)
  tgff <- .tidy_gff(tgff, prefix="t", fromAttr=tfromAttr)

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
#' @param tgff GFF file or object
#' @param fmap an optional map from target protein to gene name (fseqid)
#' @param fromAttr function mapping 'attr' column to gene name
#' @return data.frame
#' @export
make_tblastn_gene_map <- function(tblastn_file, tgff, fmap=NULL, fromAttr=identity){
  tblastn <- load_tblastn_file(tblastn_file, fmap=fmap)
  tgff <- .tidy_gff(tgff, prefix="t", fromAttr=fromAttr)
  .overlaps(
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
}

#' Load a tBLASTn file
#'
#' @param tblastn_file filename
#' @param fmap an optional map from target protein to gene name (fseqid)
#' @export
load_tblastn_file <- function(tblastn_file, fmap=NULL){
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
#' @param srcres data.frame
#' @param fmap an optional map from target protein to gene name (fseqid)
#' @export
make_tblastn_si_map <- function(tblastn_file, srcres, fmap=NULL){
  if(! is.data.frame(srcres)){
    srcres <- .tidy_searchResult(srcres)
  }
  tblastn <- load_tblastn_file(tblastn_file, fmap=fmap)
  x <- dplyr::rename(tblastn,
    # the focal protein that was BLASTed
    prot_fstart = tn_qstart,
    prot_fstop = tn_qend,
    # the target genomic intervals with translated coding similarity
    tn_tstart = tn_sstart,
    tn_tstop = tn_send
  )
  srcres <- dplyr::rename(srcres, si_fseqid=fseqid)
  .overlaps(x, srcres,
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
}

#' Load a protein BLAST file
#'
#' @param blastp_file filename
#' @param fmap data.frame: a map from blastp qseqid (a protein ID) to the seqid
#' in the SearchResult (probably a gene ID for an interval)
#' @param tmap data.frame: like fmap but for the target proteins
#' @return data.frame
#' @export
load_blastp_file <- function(blastp_file, fmap=NULL, tmap=NULL){
  d <- .read.blastfile(blastp_file, qtype="fprotein", ttype="tprotein", prefix="bp_")

  # map the target protein ID to the feature ID used in the target GFF 
  if(is.null(tmap)){
    d$tseqid <- d$tprotein
  } else {
    stopifnot(d$tprotein %in% tmap[[1]])
    d$tseqid <- .as_map(tmap)[d$tprotein]
  }
  d <- .as_first_column(d, "tseqid")

  # map the focal protein ID to the feature ID used in the query GFF 
  if(is.null(fmap)){
    d$fseqid <- d$fprotein
  } else {
    stopifnot(d$fprotein %in% fmap[[1]])
    d$fseqid <- .as_map(fmap)[d$fprotein]
  }
  d <- .as_first_column(d, "fseqid")

  d
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
#' @param fgff GFF: focal species
#' @param tgff GFF: target species
#' @param srcres SearchResult object (from \code{.tidy_searchResult})
#' @param fmap data.frame: a map from blastp qseqid (a protein ID) to the seqid
#' in the SearchResult (probably a gene ID for an interval)
#' @param tmap data.frame: like fmap but for the target proteins
#' @param tfromAttr function (optional) extract IDs attr column in target GFF
#' @param ffromAttr function (optional) extract IDs attr column in focal GFF
#' @return list containing the blast results filtered through synder and the raw blast results
#' @export
make_blastp_map <- function(blastp_file, fgff, tgff, srcres, fmap=NULL, tmap=NULL, ffromAttr=identity, tfromAttr=identity){
  fgff <- .tidy_gff(fgff, prefix="f", fromAttr=ffromAttr)
  tgff <- .tidy_gff(tgff, prefix="t", fromAttr=tfromAttr)
  x <- load_blastp_file(blastp_file, fmap, tmap) %>%
    merge(by="fseqid", fgff) %>%
    merge(by="tseqid", tgff)
  y <- dplyr::rename(srcres, si_fseqid = fseqid, si_fchr=fchr)
  .overlaps(
    x=x, xid="tchr", xa="tstart", xb="tstop",               
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
}
