#' @useDynLib synder
#' @importFrom Rcpp sourceCpp
#' @importFrom magrittr "%>%"
utils::globalVariables(c("%>%", "."))
NULL

#' synder: Trace intervals between genomes
#'
#' A synteny map is a set of linked intervals between two genomes. The primary
#' function of Synder is to use a synteny map to trace an interval in one
#' genome to a narrow search space on another genome.
#' 
#' The main functions exported by Synder are 
#'
#' \itemize{
#'   \item search - map intervals in A to search intervals in B
#'   \item anon_search - simplified version of search
#'   \item filter - find links that agree with a synteny map
#'   \item map - find intervals in B that overlap intervals in A
#'   \item count - count links in A overlapping given intervals
#'   \item dump - dump synteny map with added contiguous set ids
#' }
#'
#' Synder also defines several classes. These all have specialized plot and
#' print functions.
#'
#' \itemize{
#'    \item synmap
#'    \item gff
#'    \item hitmap
#'    \item dump_result
#'    \item search_result
#'    \item filter_result
#'    \item map_result
#'    \item count_result
#' }
#'
#' @docType package
#' @name synder
NULL

#' Synder Commands
#'
#' @section Inputs:
#'
#' The synteny map must be TAB-delimited, with no header, and must have the
#' following fields:
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
#'   * score can be any numeric value, it will be
#'   transformed as specified by the -x option
#'
#' The target and query genome lengths files must be TAB-delimited with
#' columns: <name>, <length>
#'
#' @section Search Command:
#'
#' \emph{search} predicts search intervals.
#' 
#' This is the primary function of synder. A set of intervals in one genome (the
#' query) is mapped to a set of intervals in another genome (the target). The
#' input intervals may fall between query-side syntenic intervals, in which case,
#' the search interval also will fall inbetween target syntenic intervals.
#' 
#' The output is a table with the following fields:
#' \enumerate{
#'    \item  query interval name (e.g. AT1G20300)
#'    \item  query chromosome name
#'    \item  query start position
#'    \item  query stop position
#'    \item  target chromosome name
#'    \item  search interval start position on target chromsome
#'    \item  search interval stop position on target chromsome
#'    \item  search interval strand ('+' / '-')
#'    \item  score
#'    \item contiguous set id
#'    \item lower flag
#'      \itemize{
#'        \item 0 lower bound is inside a syntenic interval
#'        \item 1 lower bound is between intervals in a contiguous set
#'        \item 2 lower bound does not overlap the contiguous set
#'        \item 3 lower bound is beyond any syntenic interval (near end of scaffold)
#'      }
#'    \item upper flag - see lower flag
#'    \item inbetween flag - TRUE if the query interval overlaps no syntenic intervals 
#' }
#'
#' @section Filter command:
#' 
#' \emph{filter} removes links that disagree with the synteny map.
#'
#' The filter function takes two syntenic maps and finds the congruent links.
#' Given two syntenic maps, A and B, the query-side intervals in A are mapped to
#' target side search intervals using the syntenic map B. Then the target-side
#' intervals in A that overlap the predicted search intervals are printed.
#' 
#' The most obvious usage case takes as input the results of BLASTing a sequence
#' against a genome. Often such searches will have many hits in a swooping e-value
#' gradient. The highest scoring hit may not be the orthologous one (for example,
#' a weakly similar, but long hit may outscore the nearly identical, but
#' truncated, true ortholog). `synder filter` will find the hits that are
#' concordant with genomic context, reducing possibly thousands of hits to only
#' a few.
#'
#' @section Map command:
#' 
#' Find all query-side syntenic blocks that overlap the input interval. Then map
#' these blocks to the target-side and print the result. If an input interval
#' overlaps no query bock, the flanks are printed.
#' 
#' The output will have the columns:o 
#' \enumerate{
#'   \item input interval name (e.g. AT1G01010)
#'   \item target contig name (e.g. Chr1)
#'   \item target start position
#'   \item target stop position
#'   \item missing flag, 0 if input overlaps no block, 1 otherwise
#' }
#'
#' @section Count command:
#'
#' Like map except it counts the number of overlaps, rather than printing them.
#' The output is a TAB-delimited list of sequence names and counts.
#'
#' @section Dump command:
#'
#' Builds the internal synteny datastructure and prints the results.
#' 
#' For example, given the file
#' 
#' \preformatted{
#' que   100    200    tar   1100   1200   100   +
#' que   1100   1800   tar   1500   1900   100   +
#' que   1400   1700   tar   1600   2000   100   +
#' que   1200   1600   tar   1700   2100   100   +
#' que   1300   1900   tar   1800   2200   100   +
#' }
#' 
#' \preformatted{
#' R> synder::dump('que-tar.syn', trans='p')
#' que  100   200   tar  1100  1200  101.000000  +  0
#' que  1100  1900  tar  1500  2200  700.084964  +  0
#' }
#' 
#' In the above example, \code{dump} shows that blocks 2-5 are merged (since
#' they are doubly-overlapping) and shows the results of the score
#' transformations. It is this dumped synteny map that would have been used in
#' any filter or search operations.
#' 
#' Also note the addition of a 9th column. This column specifies that contiguous
#' set. In this case, all blocks, after merging, are in the same set.
#'
#' @param syn synteny map file name or object
#' @param gff GFF file of input intervals
#' @param tcl target genome lengths file or object
#' @param qcl query genome lengths file or object
#' @param swap reverse direction of synteny map (target -> query)
#' @param trans synteny map score transform (Synder requires additive scores)
#'   \itemize{
#'     \item \eqn{i -> f(S) = S}            (default, no transformation)
#'     \item \eqn{d -> f(S) = L * S}        (transform from score densities)
#'     \item \eqn{p -> f(S) = L * S / 100}  (transform from percent identity)
#'     \item \eqn{l -> f(S) = -log(S)}      (transform from e-values or p-values)
#'   }
#'   Where S is input score and L interval length
#' @param k Number of interrupting intervals allowed before breaking contiguous
#' set.
#' @param r Score decay rate.
#' @param offsets Start and stop offsets (0 or 1) for synteny map, GFF file,
#' and output.
#' @param hit hit map or object
#' @name synder_commands
NULL

# NOTE: Handling of offsets:
# --------------------------
# Start and stop base offsets (0 or 1) vary between between formats.  The main
# synteny build I use, Satsuma, is 0-based relative to start, and 1-based
# relative to stop. GFF files are required to be 0-based. bioconductor (and R
# in general) is 1-based. Internally, Synder is 0-based. There are 3 sets of
# start/stop offsets to consider: synteny map base, input GFF or hitmap base,
# and output base. I will let C-side synder handle input bases, and R-side
# Synder handle output base.

do_offsets <- function(d, offsets){
  d$qstart <- d$qstart + offsets[5]
  d$qstop  <- d$qstop  + offsets[6]
  d$tstart <- d$tstart + offsets[5]
  d$tstop  <- d$tstop  + offsets[6]
  d
}

# Changes data.frames to temporary files
# FIXME: This is a bit of a hack. The old C++ code took files. Here a convert
# perfectly good R data.frames, which have already been loaded, back to files.
# I should just pass the data frames to C-side, Rcpp can handle it.
df2file <- function(x) {
  if(!is.null(x) && 'data.frame' %in% class(x)){
    xfile <- tempfile()
    readr::write_tsv(x, path = xfile, col_names = FALSE)
    x <- xfile
    class(x) <- append(class(x), 'tmp')
  }
  x
}

check_parameters <- function(
  offsets = NULL,
  k       = NULL,
  r       = NULL,
  swap    = NULL,
  trans   = NULL,
  ...
){
  stopifnot(is.null(offsets) || all(offsets %in% c(1,0)))
  stopifnot(is.null(k)       || is.numeric(k))
  stopifnot(is.null(r)       || is.numeric(r))
  stopifnot(is.null(swap)    || is.logical(swap))
  stopifnot(is.null(trans)   || trans %in% c('i', 'd', 'p', 'l'))
}

wrapper <- function(FUN, x, y=NULL, ...) {

  x <- df2file(x)
  y <- df2file(y)

  check_parameters(...)

  if(is.null(x)){
    d <- NULL
    warning("x is NULL")
  } else if(file.exists(x) && is.null(y)){
    d <- FUN(x, ...) %>% tibble::as_data_frame()
  } else if (file.exists(x) && file.exists(y)) {
    d <- FUN(x, y, ...) %>% tibble::as_data_frame()
  } else {
    warning("Failed to open files")
    d <- NULL
  }

  # Remove temporary files
  if('tmp' %in% class(x)) file.remove(x)
  if('tmp' %in% class(y)) file.remove(y)

  d
}

#' Wrapper for search allowing simple interval inputs
#'
#' The function \code{search} allows powerful mapping to target search
#' intervals using a feature table (GFF). This approach requires building a GFF
#' table, which is a needless hassle if you want to search your own intervals
#' (rather than known features). \code{anon_search} addresses this problem by
#' allowing searches given only start and stop positions and contig name.
#'
#' @param syn   synteny map file or object
#' @param a,b   start and stop locations
#' @param seqid the name of the reference query contig
#' @param ...   additional arguments sent to synder::search
#' @export
anon_search <- function(syn, a, b, seqid, ...){
  stopifnot(length(a) == c(length(b)))
  stopifnot(b >= a)
  N <- length(a)
  gff <- tibble::data_frame(
    seqid   = seqid,
    source  = NA_character_,
    type    = NA_character_,
    start   = as.integer(a),
    stop    = as.integer(b),
    score   = NA_real_,
    strand  = NA_character_,
    phase   = NA_integer_,
    attr    = paste0('seq_', 1:N)
  )
  search(syn, gff, ...)
}

#' @rdname synder_commands
#' @export
search <- function(
  syn,
  gff,
  tcl     = "",
  qcl     = "",
  swap    = FALSE,
  trans   = 'i',
  k       = 0,
  r       = 0,
  offsets = c(1,1,1,1,1,1),
  bioc    = FALSE
) {

  # If syn is a GRangePairs, try to infer the contig lengths from  the seqinfo
  # for each GRange object.
  if('GRangePairs' %in% class(syn)){
    a <- CNEr::first(syn)
    b <- CNEr::last(syn)
    if(!is.null(GenomeInfoDb::seqlengths(a)) &&
       !is.null(GenomeInfoDb::seqnames(a)))
      tcl <- GenomeInfoDb::seqinfo(a)
    if(!is.null(GenomeInfoDb::seqlengths(b)) &&
       !is.null(GenomeInfoDb::seqnames(b)))
      qcl <- GenomeInfoDb::seqinfo(b)
  }

  if(!(is.character(tcl) && tcl == "")) tcl <- as_conlen(tcl) 
  if(!(is.character(qcl) && qcl == "")) qcl <- as_conlen(qcl) 

  d <- wrapper(
    FUN     = c_search,
    x       = as_synmap(syn),
    y       = as_gff(gff),
    tcl     = df2file(tcl),
    qcl     = df2file(qcl),
    swap    = swap,
    k       = k,
    r       = r,
    trans   = trans,
    offsets = offsets[1:4]
  )

  d$attr   <- as.character(d$attr)
  d$qseqid <- as.character(d$qseqid)
  d$tseqid <- as.character(d$tseqid)
  d$strand <- as.character(d$strand)

  class(d) <- append('search_result', class(d))
  attributes(d)$swap  <- swap
  attributes(d)$k     <- k
  attributes(d)$r     <- r
  attributes(d)$trans <- trans

  d <- do_offsets(d, offsets)

  if(bioc)
    d <- as_bioc(d, seqinfo_a=as_bioc(qcl), seqinfo_b=as_bioc(tcl))

  d
}


#' @rdname synder_commands
#' @export
filter <- function(
  syn,
  hit,
  swap    = FALSE,
  trans   = 'i',
  k       = 0,
  r       = 0,
  offsets = c(1,1,1,1,1,1)
) {
  d <- wrapper(
    FUN     = c_filter,
    x       = as_synmap(syn),
    y       = hit,
    swap    = swap,
    k       = k,
    r       = r,
    trans   = trans,
    offsets = offsets[1:4]
  )
  if(is.null(d)) return(NULL)

  d <- sub(d, pattern="\n", replacement="") %>%
    strsplit(split="\t")                    %>%
    do.call(what=rbind)                     %>%
    tibble::as_data_frame()
  names(d)[1:6] <- names(SYNMAP_COLS)[1:6]
  d$qstart <- as.numeric(d$qstart)
  d$qstop  <- as.numeric(d$qstop)
  d$tstart <- as.numeric(d$tstart)
  d$tstop  <- as.numeric(d$tstop)

  class(d) <- append('filter_result', class(d))
  attributes(d)$swap  <- swap
  attributes(d)$k     <- k
  attributes(d)$r     <- r
  attributes(d)$trans <- trans

  d <- do_offsets(d, offsets)

  d
}


#' @rdname synder_commands
#' @export
map <- function(
  syn,
  gff,
  swap    = FALSE,
  offsets = c(1,1,1,1,1,1)
) {
  d <- wrapper(
    FUN     = c_map,
    x       = as_synmap(syn),
    y       = as_gff(gff),
    swap    = swap,
    offsets = offsets[1:4]
  )

  class(d) <- append('map_result', class(d))
  attributes(d)$swap <- swap

  d <- do_offsets(d, offsets)

  d
}


#' @rdname synder_commands
#' @export
count <- function(
  syn,
  gff,
  swap    = FALSE,
  offsets = c(1,1,1,1,1,1)
) {
  d <- wrapper(
    FUN     = c_count,
    x       = as_synmap(syn),
    y       = as_gff(gff),
    swap    = swap,
    offsets = offsets[1:4]
  )

  d$attr = as.character(d$attr)

  class(d) <- append('count_result', class(d))
  attributes(d)$swap <- swap

  d
}


#' @rdname synder_commands
#' @export
dump <- function(
  syn,
  swap    = FALSE,
  trans   = 'i',
  offsets = c(1,1,1,1,1,1)
) {
  d <- wrapper(
    FUN     = c_dump,
    x       = as_synmap(syn),
    swap    = swap,
    trans   = trans,
    offsets = offsets[1:4]
  )
  d$qseqid <- as.character(d$qseqid)
  d$tseqid <- as.character(d$tseqid)

  # Assign class and attributes
  class(d) <- append('dump_result', class(d))
  attributes(d)$swap  = swap
  attributes(d)$trans = trans

  d <- do_offsets(d, offsets)

  d
}
