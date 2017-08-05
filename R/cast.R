#' Synder cast functions
#'
#' cast data as Synder objects
#'
#' @param d input type
#' @param ... additional arguments
#' @name synder_cast
NULL

.pairs_to_synmap <- function(x, y, strand, score){
  if(! ('score' %in% names(GenomicRanges::mcols(y)))){
    stop('A "score" field is required in the second GRanges object ',
         'for synteny maps expressed as GRangesPairs objects.')
  }

  data.frame(
    qseqid = as.character(GenomicRanges::seqnames(x)),
    qstart = GenomicRanges::start(x),
    qstop  = GenomicRanges::end(x),
    tseqid = as.character(GenomicRanges::seqnames(y)),
    tstart = GenomicRanges::start(y),
    tstop  = GenomicRanges::end(y),
    score  = GenomicRanges::mcols(y)$score,
    strand = as.character(BiocGenerics::strand(y)),
    stringsAsFactors=FALSE
  )
}

#' @rdname synder_cast
#' @export
as_synmap <- function(d) {

  if('synmap' %in% class(d)) return(d)

  d <- if(("Axt" %in% class(d)) || ("GRangePairs" %in% class(d))){
    .pairs_to_synmap(x=CNEr::first(d), y=CNEr::last(d))
  } else if (is.character(d)){
    if(file.exists(d)){
      read_synmap(d)
    } else {
      stop(sprintf("Cannot find synteny map file '%s'", d))
    }
  } else {
    as.data.frame(d)
  }

  # assert correct number of columns
  stopifnot(ncol(d) == 8)

  # Set column names
  names(d) <- names(SYNMAP_COLS)

  # Assert column types match expectations
  stopifnot(SYNMAP_COLS == lapply(d, class))

  # Assign class membership
  class(d) <- append('synmap', class(d))

  return(d)
}

.maybe_meta <- function(x, field, default=NA, caster=identity){
  if(field %in% names(GenomicRanges::mcols(x))){
    caster(GenomicRanges::mcols(x)[[field]])
  } else {
    default
  }
}

.GRanges_to_GFF <- function(
  x,
  source_tag = "source",
  type_tag   = "type",
  score_tag  = "score",
  phase_tag  = "phase",
  id_tag     = "name"
){
  data.frame(
    seqid   = as.character(GenomicRanges::seqnames(x)),
    source  = .maybe_meta(x, source_tag, NA_character_, as.character),
    type    = .maybe_meta(x, type_tag, NA_character_, as.character),
    start   = GenomicRanges::start(x),
    stop    = GenomicRanges::end(x),
    score   = .maybe_meta(x, score_tag, NA_real_, as.numeric),
    strand  = as.character(GenomicRanges::strand(x)),
    phase   = .maybe_meta(x, phase_tag, NA_integer_, as.integer),
    attr    = .maybe_meta(x, id_tag, NA_character_, as.character),
    stringsAsFactors=FALSE
  )
}

#' @rdname synder_cast
#' @export
as_gff <- function(d, ...) {

  if('gff' %in% class(d)) return(d)

  d <- if("GRanges" %in% class(d)){
    .GRanges_to_GFF(d, ...)
  } else if (is.character(d)){
    if(file.exists(d)){
      read_gff(d)
    } else {
      stop(sprintf("Cannot find GFF file '%s'", d))
    }
  } else {
    as.data.frame(d)
  }

  # Assert the correct number of columns were read
  stopifnot(ncol(d) == 9)

  # Set column names
  names(d) <- names(GFF_COLS)

  # Assert column types match expectations
  for(i in 1:ncol(d)){
    if(class(d[[i]]) != GFF_COLS[i]){
      stop(paste("GFFError: in column", i, "expected type", GFF_COLS[i],
                 ", found", class(d[[i]])))
    }
  }

  # Assign class membership
  class(d) <- append('gff', class(d))

  return(d)
}

#' @rdname synder_cast
#' @export
as_hitmap <- function(d) {
  d <- as.data.frame(d)

  # Assert the correct number of columns were read
  stopifnot(ncol(d) >= 6)

  # Set column names
  names(d)[1:6] <- names(SYNMAP_COLS)[1:6]

  # Assert column types match expectations
  stopifnot(SYNMAP_COLS[1:6] == lapply(d[,1:6], class))

  # Assign class membership
  class(d) <- append('hitmap', class(d))

  return(d)
}

#' @rdname synder_cast
#' @export
as_conlen <- function(d) {

  if('conlen' %in% class(d)) return(d)

  d <- if("Seqinfo" %in% class(d)){
    data.frame(
      seqid = as.character(GenomeInfoDb::seqnames(d)),
      length = GenomeInfoDb::seqlengths(d),
      stringsAsFactors=FALSE
    )
  } else if (is.character(d)){
    if(file.exists(d)){
      read_conlen(d)
    } else {
      stop(sprintf("Cannot find genome lengths file '%s'", d))
    }
  } else {
    as.data.frame(d)
  }

  # Assert the correct number of columns were read
  stopifnot(ncol(d) == 2)

  # Set column names
  names(d) <- names(CON_LENGTH)

  # Assert column types match expectations
  stopifnot(CON_LENGTH == lapply(d, class))

  # Assign class membership
  class(d) <- append('conlen', class(d))

  return(d)
}
