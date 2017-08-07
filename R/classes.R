Synmap <- setClass(
  'Synmap',
  contains = 'GRangePairs'
)

GFF <- setClass(
  'GFF',
  contains = 'GRanges'
)

DumpResult <- setClass(
  'DumpResult',
  contains = 'GRangePairs',
  slots = list(
    swap    = 'logical',
    trans   = 'character',
    offsets = 'integer'
  )
)
DumpResult <- function(x, offsets=c(1,1,1,1,1,1), ...){
  # TODO: add checking
  new("DumpResult", x, offsets=as.integer(offsets), ...)
}

setClass(
  'SearchResult',
  contains = 'GRangePairs',
  slots = list(
    swap    = 'logical',
    trans   = 'character',
    k       = 'integer',
    r       = 'numeric',
    offsets = 'integer'
  )
)
SearchResult <- function(x, offsets=c(1,1,1,1,1,1), k=0, ...){
  # TODO: add checking
  new("SearchResult", x, offsets=as.integer(offsets), k=as.integer(k), ...)
}
