#' Synder cast functions
#'
#' cast data as Synder objects
#'
#' @param d input type
#' @param ... additional arguments
#' @name synder_cast
NULL

as_synder_data_frame <- function(x, ...){
  UseMethod('as_synder_data_frame', x)
}

as_synder_data_frame.GRangePairs <- function(x, ...){

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

  a_meta <- GenomicRanges::mcols(a)
  b_meta <- GenomicRanges::mcols(b)

  if(ncol(a_meta) > 0)
    d <- cbind(d, a_meta)
  if(ncol(b_meta) > 0)
    d <- cbind(d, b_meta)

  # For some silly reason, Bioconductor uses '*' for unknown strand, even
  # though '.' is the GFF convention. Synder also considers '.' as missig data.
  # So here I convert back.
  str <- as.character(BiocGenerics::strand(b))
  str <- ifelse(str == '*', '.', str)
  d$strand = str

  if(setequal(names(d), names(SI_COLS))){
    d <- d[, names(SI_COLS)]
  }

  # Convert from DataFrame (some Bioc nonsense) to the normal data.frame
  as.data.frame(d)
}

as_synder_data_frame.Seqinfo <- function(x, ...){
  data.frame(
    seqid = as.character(GenomeInfoDb::seqnames(x)),
    length = GenomeInfoDb::seqlengths(x),
    stringsAsFactors=FALSE
  )
}

as_synder_data_frame.GRanges <- function(
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

as_synmap <- function(x, ...){
  UseMethod('as_synmap', x)
}

as_synmap.character <- function(x, ...){
  if(file.exists(x)){
    read_synmap(x)
  } else {
    stop(sprintf("Cannot read synmap file '%s'", x))
  }
}

as_synmap.Axt <- function(x, seqinfo_a=NULL, seqinfo_b=NULL){
  a = CNEr::queryRanges(x)
  b = CNEr::targetRanges(x)

  if(seqinfo(a) == NULL)
    GenomeInfoDb::seqinfo(a) <- seqinfo_a
  if(seqinfo(b) == NULL)
    GenomeInfoDb::seqinfo(b) <- seqinfo_b

  GenomicRanges::mcols(b)$score <- CNEr::score(x)
  # NOTE: CNEr stores the relative strand in the query
  GenomicRanges::mcols(b)$strand <- GenomicRanges::strand(a)
  CNEr::GRangePairs(first=a, second=b)
}

as_synmap.data.frame <- function(x, seqinfo_a=NULL, seqinfo_b=NULL) {
  CNEr::GRangePairs(
    .make_GRanges(
      seqnames=x$qseqid,
      start=x$qstart,
      stop=x$qstop,
      seqinfo=seqinfo_a
    ),
    .make_GRanges(
      seqnames=x$tseqid,
      start=x$tstart,
      stop=x$tstop,
      seqinfo=seqinfo_b,
      score=x$score,
      strand=x$strand
    )
  )
}

as_synmap.GRangePairs <- function(x, seqinfo_a=NULL, seqinfo_b=NULL){
  if(!is.null(seqinfo_a))
    first(x)$seqinfo <- seqinfo_a
  if(!is.null(seqinfo_b))
    second(x)$seqinfo <- seqinfo_b
  x
}

as_gff <- function(x, ...){
  UseMethod('as_gff', x)
}

as_gff.character <- function(x, ...){
  if(file.exists(x)){
    read_gff(x)
  } else {
    stop(sprintf("Cannot read gff file '%s'", x))
  }
}

as_gff.GRanges <- function(x, seqinfo_=NULL) {
  if(is.null(seqinfo(x))){
    seqinfo(x) <- seqinfo_
  }
  x
}

as_gff.data.frame <- function(x, seqinfo=FALSE){
  .make_GRanges(
    seqnames = x$seqid,
    start    = x$start,
    stop     = x$stop,
    source   = ifelse(x$source == '.', NA_character_, x$source),
    type     = ifelse(x$type   == '.', NA_character_, x$type),
    score    = ifelse(x$score  == '.', NA_real_,      x$score),
    strand   = ifelse(x$strand == '.', '*',           x$strand),
    phase    = ifelse(x$phase  == '.', NA_integer_,   x$phase),
    attr     = x$attr
  )
}

as_hitmap <- function(x){
  UseMethod('as_hitmap', x) 
}

as_hitmap.data.frame <- function(x) {
  as_bioc_hitmap(x)
}

as_hitmap.GRanges <- function(x) {
  x
}


as_conlen <- function(x) {
  UseMethod('as_conlen', x)
}

as_conlen.character <- function(x, ...){
  if(file.exists(x)){
    read_conlen(x)
  } else {
    stop(sprintf("Cannot read conlen file '%s'", x))
  }
}

as_conlen.Seqinfo <- function(x) {
  x
}

as_conlen.data.frame <- function(x, isCircular=NA, genome=NA) {
  GenomeInfoDb::Seqinfo(
    seqnames   = x$seqid,
    seqlengths = x$length,
    isCircular = isCircular,
    genome     = genome
  )
}


.make_GRanges <- function(seqnames, start, stop, strand='+', ...){
  GenomicRanges::GRanges(
    seqnames = seqnames,
    strand = strand,
    ranges = IRanges::IRanges(start=start, end=stop),
    ...
  )
}
