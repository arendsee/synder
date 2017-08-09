#' Synder cast functions
#'
#' Cast data as synder objects
#'
#' \code{as_synmap}
#' \code{as_gff}
#' \code{as_conlen}
#'
#' @param d input type
#' @param ... additional arguments
#' @name synder_cast
NULL

.make_GRanges <- function(seqnames, start, stop, strand='+', ...){
  GenomicRanges::GRanges(
    seqnames = seqnames,
    strand = strand,
    ranges = IRanges::IRanges(start=start, end=stop),
    ...
  )
}

.base_GRangePairs_to_df <- function(x, ordering=NULL){

  a <- CNEr::first(x)
  b <- CNEr::second(x)

  if(any(GenomicRanges::strand(a) != '+')){
    stop(
      "Cannot have a negative strand in a GRangePairs' first GRange object. ",
      "In synder, strand is relative to the first genome, so only the second ",
      "genome will have negative sense."
    )
  }

  d <- data.frame(
    qseqid = as.character(GenomicRanges::seqnames(a)),
    qstart = GenomicRanges::start(a),
    qstop  = GenomicRanges::end(a),
    tseqid = as.character(GenomicRanges::seqnames(b)),
    tstart = GenomicRanges::start(b),
    tstop  = GenomicRanges::end(b),
    stringsAsFactors=FALSE
  )

  a_meta <- GenomicRanges::mcols(a) %>% as.data.frame
  b_meta <- GenomicRanges::mcols(b) %>% as.data.frame
  c_meta <- GenomicRanges::mcols(x) %>% as.data.frame

  # For some silly reason, Bioconductor uses '*' for unknown strand, even
  # though '.' is the GFF convention. Synder also considers '.' as missig data.
  # So here I convert back.
  if('strand' %in% names(c_meta)){
    c_meta$strand <- as.character(c_meta$strand)
    c_meta$strand <- ifelse(c_meta$strand == '*', '.', c_meta$strand)
  }

  if(ncol(a_meta) > 0)
    d <- cbind(d, a_meta)
  if(ncol(b_meta) > 0)
    d <- cbind(d, b_meta)
  if(ncol(c_meta) > 0)
    d <- cbind(d, c_meta)

  if(!is.null(ordering)){
    stopifnot(all(ordering %in% names(d)))
    d <- d[, ordering]
  }

  # Convert from DataFrame (some Bioc nonsense) to the normal data.frame
  as.data.frame(d)
}

as.data.frame.Synmap <- function(x, ...){
  .base_GRangePairs_to_df(x, ordering=names(SYNMAP_COLS))
}

as.data.frame.DumpResult <- function(x, ...){
  .base_GRangePairs_to_df(x, ordering=names(DUMP_COLS))
}

as.data.frame.SearchResult <- function(x, ...){
  .base_GRangePairs_to_df(x, ordering=names(SI_COLS))
}

as.data.frame.Seqinfo <- function(x, ...){
  data.frame(
    seqid = as.character(GenomeInfoDb::seqnames(x)),
    length = GenomeInfoDb::seqlengths(x),
    stringsAsFactors=FALSE
  )
}

as.data.frame.GFF <- function(
  x,
  source_tag = "source",
  type_tag   = "type",
  score_tag  = "score",
  phase_tag  = "phase",
  id_tag     = "attr",
  ...
){
  .maybe_meta <- function(x, field, default=NA, caster=identity){
    if(field %in% names(GenomicRanges::mcols(x))){
      caster(GenomicRanges::mcols(x)[[field]])
    } else {
      default
    }
  }

  strand <- as.character(GenomicRanges::strand(x))
  strand <- ifelse(strand == '*', '.', strand)

  data.frame(
    seqid  = as.character(GenomicRanges::seqnames(x)),
    source = .maybe_meta(x, source_tag, NA_character_, as.character),
    type   = .maybe_meta(x, type_tag, NA_character_, as.character),
    start  = GenomicRanges::start(x),
    stop   = GenomicRanges::end(x),
    score  = .maybe_meta(x, score_tag, NA_real_, as.numeric),
    strand = strand,
    phase  = .maybe_meta(x, phase_tag, NA_integer_, as.integer),
    attr   = .maybe_meta(x, id_tag, NA_character_, as.character),
    stringsAsFactors=FALSE
  )
}



#' @rdname synder_cast
#' @export
as_synmap <- function(x, ...){
  UseMethod('as_synmap', x)
}

as_synmap.Synmap <- function(x, ...) x

as_synmap.character <- function(x, ...){
  if(file.exists(x)){
    read_synmap(x, ...)
  } else {
    stop(sprintf("Cannot read synmap file '%s'", x))
  }
}

as_synmap.Axt <- function(x, seqinfo_a=NULL, seqinfo_b=NULL){
  a = CNEr::queryRanges(x)
  b = CNEr::targetRanges(x)

  if(seqinfo(a) == NULL)
    GenomeInfoDb::seqinfo(a) <- as_conlen(seqinfo_a)
  if(seqinfo(b) == NULL)
    GenomeInfoDb::seqinfo(b) <- as_conlen(seqinfo_b)

  Synmap(CNEr::GRangePairs(
    first  = a,
    second = b,
    score  = CNEr::score(x),
    strand = GenomicRanges::strand(a) # NOTE: strand stored relative to query
  ))
}

as_synmap.data.frame <- function(x, seqinfo_a=NULL, seqinfo_b=NULL) {
  Synmap(CNEr::GRangePairs(
    .make_GRanges(
      seqnames=x$qseqid,
      start=x$qstart,
      stop=x$qstop,
      seqinfo=as_conlen(seqinfo_a)
    ),
    .make_GRanges(
      seqnames=x$tseqid,
      start=x$tstart,
      stop=x$tstop,
      seqinfo=as_conlen(seqinfo_b)
    ),
    score=x$score,
    strand=x$strand
  ))
}

as_synmap.GRangePairs <- function(x, seqinfo_a=NULL, seqinfo_b=NULL){
  if(!is.null(seqinfo_a))
    first(x)$seqinfo <- as_conlen(seqinfo_a)
  if(!is.null(seqinfo_b))
    second(x)$seqinfo <- as_conlen(seqinfo_b)
  Synmap(x)
}


#' @rdname synder_cast
#' @export
as_gff <- function(x, ...){
  UseMethod('as_gff', x)
}

as_gff.GFF <- function(x, ...) x

as_gff.character <- function(x, ...){
  if(file.exists(x)){
    read_gff(x, ...)
  } else {
    stop(sprintf("Cannot read gff file '%s'", x))
  }
}

as_gff.GRanges <- function(x, seqinfo_=NULL) {
  if(!is.null(seqinfo_)){
    seqinfo(x) <- seqinfo_
  }
  GFF(x)
}

as_gff.data.frame <- function(x, seqinfo=NULL){
  GFF(.make_GRanges(
    seqnames = x$seqid,
    start    = x$start,
    stop     = x$stop,
    source   = ifelse(x$source == '.', NA_character_, x$source),
    type     = ifelse(x$type   == '.', NA_character_, x$type),
    score    = ifelse(x$score  == '.', NA_real_,      x$score),
    strand   = ifelse(x$strand == '.', '*',           x$strand),
    phase    = ifelse(x$phase  == '.', NA_integer_,   x$phase),
    attr     = x$attr,
    seqinfo  = as_conlen(seqinfo)
  ))
}


#' @rdname synder_cast
#' @export
as_conlen <- function(x, ...) {
  UseMethod('as_conlen', x)
}

as_conlen.NULL <- function(x, ...) NULL

as_conlen.Seqinfo <- function(x, ...) x

as_conlen.character <- function(x, ...){
  if(file.exists(x)){
    read_conlen(x, ...)
  } else {
    stop(sprintf("Cannot read conlen file '%s'", x))
  }
}

as_conlen.data.frame <- function(x, ...) {
  GenomeInfoDb::Seqinfo(
    seqnames   = x$seqid,
    seqlengths = x$length,
    ...
  )
}
