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


.synder_is <- function(x, types, bioc_type){
  (is.data.frame(x) && all(names(x) == names(types))) || class(x) == bioc_type
}

is_synmap <- function(x, ...){
  .synder_is(x, SYNMAP_COLS, 'Synmap', ...)
}

is_gff <- function(x, ...){
  .synder_is(x, GFF_COLS, 'GFF', ...)
}

is_search_result <- function(x, ...){
  .synder_is(x, SI_COLS, 'SearchResult', ...)
}

is_dump <- function(x, ...){
  .synder_is(x, DUMP_COLS, 'DumpResult', ...)
}

is_conlen <- function(x, ...){
  .synder_is(x, CON_LENGTH, 'Seqinfo', ...)
}
