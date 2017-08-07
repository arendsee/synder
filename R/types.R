SYNMAP_COLS <- c(
  "qseqid" = "character",
  "qstart" = "integer",
  "qstop"  = "integer",
  "tseqid" = "character",
  "tstart" = "integer",
  "tstop"  = "integer",
  "score"  = "numeric",
  "strand" = "character"
)

GFF_COLS <- c(
  "seqid"   = "character",
  "source"  = "character",
  "type"    = "character",
  "start"   = "integer",
  "stop"    = "integer",
  "score"   = "numeric",
  "strand"  = "character",
  "phase"   = "integer",
  "attr"    = "character"
)

CON_LENGTH <- c(
  "seqid"  = "character",
  "length" = "integer"
)

SI_COLS <- c(
  "attr"      = "character",
  "qseqid"    = "character",
  "qstart"    = "integer",
  "qstop"     = "integer",
  "tseqid"    = "character",
  "tstart"    = "integer",
  "tstop"     = "integer",
  "strand"    = "character",
  "score"     = "numeric",
  "cset"      = "integer",
  "l_flag"    = "integer",
  "r_flag"    = "integer",
  "inbetween" = "logical"
)

DUMP_COLS <- c(
  "qseqid" = "character",
  "qstart" = "integer",
  "qstop"  = "integer",
  "tseqid" = "character",
  "tstart" = "integer",
  "tstop"  = "integer",
  "score"  = "numeric",
  "strand" = "logical",
  "cset"   = "integer"
)


.synder_is <- function(x, types, bioc_base_type, ordered=TRUE){
  if(is.data.frame(x)){
    all(names(x) == names(types))
  } else if(bioc_base_type %in% class(x)) {
    if(ordered){
      all(names(as_synder_data_frame(x)) == names(types))
    } else {
      setequal(names(as_synder_data_frame(x)), names(types))
    }
  } else {
    FALSE
  }
}

is_synmap <- function(x, ...){
  .synder_is(x, SYNMAP_COLS, 'GRangePairs', ...)
}

is_gff <- function(x, ...){
  .synder_is(x, GFF_COLS, 'GRanges', ...)
}

is_search_result <- function(x, ...){
  .synder_is(x, SI_COLS, 'GRangePairs', ...)
}

is_dump <- function(x, ...){
  .synder_is(x, DUMP_COLS, 'GRangePairs', ...)
}

is_conlen <- function(x, ...){
  .synder_is(x, CON_LENGTH, 'Seqinfo', ...)
}
