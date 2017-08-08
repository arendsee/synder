#' @importClassesFrom S4Vectors Pairs Vector Annotated
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom CNEr GRangePairs
NULL

#' Hold a synteny map
#'
#' @exportClass Synmap
Synmap <- setClass(
  'Synmap',
  contains = 'GRangePairs'
)

#' Hold a GFF map
#'
#' @exportClass GFF
GFF <- setClass(
  'GFF',
  contains = 'GRanges'
)

#' Hold a GFF map
#'
#' @exportClass DumpResult
DumpResult <- setClass(
  'DumpResult',
  contains = 'GRangePairs',
  slots = list(
    swap    = 'logical',
    trans   = 'character',
    offsets = 'integer'
  ),
  prototype = list(
    swap    = FALSE,
    trans   = 'i',
    offsets = c(1L,1L,1L,1L,1L,1L)
  )
)

#' Hold a SearchResult
#'
#' @exportClass SearchResult
SearchResult <- setClass(
  'SearchResult',
  contains = 'GRangePairs',
  slots = list(
    swap    = 'logical',
    trans   = 'character',
    k       = 'integer',
    r       = 'numeric',
    offsets = 'integer'
  ),
  prototype = list(
    swap    = FALSE,
    trans   = 'i',
    k       = 0L,
    r       = 0,
    offsets = c(1L,1L,1L,1L,1L,1L)
  )
)
