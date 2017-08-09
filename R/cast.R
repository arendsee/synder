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

.as_bioc_strand <- function(x){
  x <- as.character(x)
  ifelse(x == '.', '*', x)
}

.as_gff_strand <- function(x){
  x <- as.character(x)
  ifelse(x == '*', '.', x)
}

.make_GRanges <- function(seqnames, start, stop, strand='+', ...){
  GenomicRanges::GRanges(
    seqnames = seqnames,
    strand = .as_bioc_strand(strand),
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
    c_meta$strand <- .as_gff_strand(c_meta$strand)
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


#' @rdname synder_cast
#' @method as.data.frame Synmap
#' @export
as.data.frame.Synmap <- function(x, ...){
  .base_GRangePairs_to_df(x, ordering=names(SYNMAP_COLS))
}

#' @rdname synder_cast
#' @method as.data.frame DumpResult
#' @export
as.data.frame.DumpResult <- function(x, ...){
  .base_GRangePairs_to_df(x, ordering=names(DUMP_COLS))
}

#' @rdname synder_cast
#' @method as.data.frame SearchResult
#' @export
as.data.frame.SearchResult <- function(x, ...){
  .base_GRangePairs_to_df(x, ordering=names(SI_COLS))
}

#' @rdname synder_cast
#' @method as.data.frame Seqinfo
#' @export
as.data.frame.Seqinfo <- function(x, ...){
  data.frame(
    seqid = as.character(GenomeInfoDb::seqnames(x)),
    length = GenomeInfoDb::seqlengths(x),
    stringsAsFactors=FALSE
  )
}

#' @rdname synder_cast
#' @method as.data.frame GFF
#' @export
as.data.frame.GFF <- function(x, ...){

  strand <- .as_gff_strand(GenomicRanges::strand(x))

  met <- GenomicRanges::mcols(x)

  data.frame(
    seqid  = as.character(GenomicRanges::seqnames(x)),
    source = met$source,
    type   = met$type,
    start  = GenomicRanges::start(x),
    stop   = GenomicRanges::end(x),
    score  = met$score,
    strand = strand,
    phase  = met$phase,
    attr   = met$attr,
    stringsAsFactors=FALSE
  )
}

#' Convert to a synder Synmap object
#'
#' @export
#' @param x Thing to be converted
#' @param seqinfo_a Seqinfo object for first GRanges entry
#' @param seqinfo_b Seqinfo object for second GRanges entry
#' @param ... Additional arguments
as_synmap <- function(x, ...){
  UseMethod('as_synmap', x)
}

#' @rdname as_synmap
#' @export
as_synmap.Synmap <- function(x, ...) x

#' @rdname as_synmap
#' @export
as_synmap.character <- function(x, ...){
  if(file.exists(x)){
    read_synmap(x, ...)
  } else {
    stop(sprintf("Cannot read synmap file '%s'", x))
  }
}

#' @rdname as_synmap
#' @export
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
    # NOTE: strand stored relative to query
    strand = .as_bioc_strand(GenomicRanges::strand(a))
  ))
}

#' @rdname as_synmap
#' @export
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
    score  = x$score,
    strand = .as_bioc_strand(x$strand)
  ))
}

#' @rdname as_synmap
#' @export
as_synmap.GRangePairs <- function(x, seqinfo_a=NULL, seqinfo_b=NULL){
  if(!is.null(seqinfo_a))
    first(x)$seqinfo <- as_conlen(seqinfo_a)
  if(!is.null(seqinfo_b))
    second(x)$seqinfo <- as_conlen(seqinfo_b)
  Synmap(x)
}


#' Convert to a synder GFF object
#'
#' as_gff.character will attempt to load a file, which is expected to be
#' TAB-delimited and headerless.
#'
#' @export
#' @param x Thing to be converted
#' @param ... Additional arguments
#' @param seqinfo_ Seqinfo object
#' @param source character vector for GFF column source meta-column
#' @param type character vector for GFF column type meta-column
#' @param score numeric vector for GFF column score meta-column
#' @param phase integer vector for GFF column phase meta-column
#' @param id character vector for GFF column attr, is currently used as the
#'        name of the query
as_gff <- function(x, ...){
  UseMethod('as_gff', x)
}

#' @rdname as_gff
#' @export
as_gff.GFF <- function(x, ...) x

#' @rdname as_gff
#' @export
as_gff.GRanges <- function(
  x,
  seqinfo_ = NULL,
  source   = NULL,
  type     = NULL,
  score    = NULL,
  phase    = NULL,
  id       = NULL
) {
  if(is.null(seqinfo_))
    seqinfo_ <- GenomicRanges::seqinfo(x)

  .maybe_meta <- function(x, field, default=NA, caster=identity){
    if(is.null(field)){
      rep(default, length(x))
    } else {
      caster(field)
    }
  }

  GFF(
    GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(x),
      ranges   = GenomicRanges::ranges(x),
      seqinfo  = seqinfo_,
      strand   = .as_bioc_strand(GenomicRanges::strand(x)),
      source   = .maybe_meta(x, source, NA_character_, as.character),
      type     = .maybe_meta(x, type,   NA_character_, as.character),
      score    = .maybe_meta(x, score,  NA_real_,      as.numeric),
      phase    = .maybe_meta(x, phase,  NA_integer_,   as.integer),
      attr     = .maybe_meta(x, id,     NA_character_, as.character)
    )
  )
}

#' @rdname as_gff
#' @export
as_gff.character <- function(x, ...){
  if(file.exists(x)){
    read_gff(x, ...)
  } else {
    stop(sprintf("Cannot read gff file '%s'", x))
  }
}

#' @rdname as_gff
#' @export
as_gff.data.frame <- function(x, seqinfo=NULL){
  GFF(.make_GRanges(
    seqnames = x$seqid,
    start    = x$start,
    stop     = x$stop,
    source   = as.character(ifelse(x$source == '.', NA_character_, x$source)),
    type     = as.character(ifelse(x$type   == '.', NA_character_, x$type)),
    score    = as.numeric(ifelse(x$score    == '.', NA_real_,      x$score)),
    strand   = .as_bioc_strand(x$strand),
    phase    = as.integer(ifelse(x$phase    == '.', NA_integer_,   x$phase)),
    attr     = as.character(x$attr),
    seqinfo  = as_conlen(seqinfo)
  ))
}


#' Convert to a contig length (Seqinfo) object
#'
#' @export
#' @param x Thing to be converted
#' @param ... Additional arguments
as_conlen <- function(x, ...) {
  UseMethod('as_conlen', x)
}

#' @rdname as_conlen
#' @export
as_conlen.NULL <- function(x, ...) NULL

#' @rdname as_conlen
#' @export
as_conlen.Seqinfo <- function(x, ...) x

#' @rdname as_conlen
#' @export
as_conlen.character <- function(x, ...){
  if(file.exists(x)){
    read_conlen(x, ...)
  } else {
    stop(sprintf("Cannot read conlen file '%s'", x))
  }
}

#' @rdname as_conlen
#' @export
as_conlen.data.frame <- function(x, ...) {
  GenomeInfoDb::Seqinfo(
    seqnames   = x$seqid,
    seqlengths = x$length,
    ...
  )
}
