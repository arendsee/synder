#' @importClassesFrom S4Vectors Pairs Vector Annotated
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom CNEr GRangePairs
NULL


#' Synder classes
#'
#' Synder uses the input classes \code{Synmap}, \code{GFF} and \code{Seqinfo}
#' (the last of which is from bioconductor); and the output classes
#' SearchResult and DumpResult.
#'
#' Each of these classes can be converted from a \code{data.frame} with
#' \code{as_*} and to a data.frame with \code{as.data.frame}.
#'
#' @section Synteny maps and the Synmap object:
#'
#' As a data.frame, a synteny map must have the following columns:
#'
#' \enumerate{
#'   \item qseqid - query contig id (e.g. Chr1)
#'   \item qstart - query interval start
#'   \item qstop  - query interval stop
#'   \item sseqid - target contig id
#'   \item sstart - target interval start
#'   \item sstop  - target interval stop
#'   \item score  - score of the syntenic match*
#'   \item strand - relative orientation
#' }
#'
#' \code{score} can be any numeric value, it will be transformed as specified
#' by the -x option
#'
#' It can be from a file with the, where it is required to be TAB-delimited
#' with no header.
#'
#' The target and query genome lengths files must be TAB-delimited with
#' columns: <name>, <length>
#'
#' The \code{as_synmap} function converts such a data.frame to a Synmap object.
#' A Synmap object is a GRangePairs object (from the bioconductor CNEr package)
#' with socre and strand as meta-columns.
#'
#' @section GFF files and the GFF object:
#'
#' GFF files are used to represent queries that will be mapped against the
#' synteny map. These can be loaded into GFF objects with read_gff or converted
#' from data.frames with as_gff. Either way, the resul is a GRange object with
#' \code{source}, \code{type}, \code{score}, \code{phase}, and \code{attr}
#' meta-columns.
#'
#' @section Contig length files and Seqinfo objects:
#'
#' Properly handling boundary cases when mapping requires knowing the length of
#' each contig. This information is provided by the biodonductor Seqinfo
#' objects. The most important information these objects hold is the name and
#' length of all contigs in the assembly. This data can be stored in a
#' two-column data.frame (or headerless, TAB-delimited file).
#'
#' @section SearchResult:
#'
#' Holds the results of a successful \code{search} run. The class inherits from GRangePairs.
#'
#' @section DumpResult:
#'
#' Holds the results of the \code{dump} command. The class inherits from GRangePairs.
#'
#' @name synder_classes
NULL

#' @rdname synder_classes
#' @exportClass Synmap
Synmap <- setClass(
  'Synmap',
  contains = 'GRangePairs'
)

#' @rdname synder_classes
#' @exportClass GFF
GFF <- setClass(
  'GFF',
  contains = 'GRanges'
)

#' @rdname synder_classes
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
    offsets = c(1L,1L)
  )
)

#' @rdname synder_classes
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
    offsets = c(1L,1L)
  )
)
