#' @useDynLib synder
#' @importFrom methods new
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
#' The primary function exported by Synder are
#'
#' \itemize{
#'   \item search - map intervals in A to search intervals in B
#'   \item anon_search - simplified version of search
#'   \item dump - dump synteny map with added contiguous set ids
#' }
#'
#' @docType package
#' @name synder
NULL

#' Synder Commands
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
  if(!is.character(x) && !is.null(x)){ 
    xfile <- tempfile()
    readr::write_tsv(as.data.frame(x), path = xfile, col_names = FALSE)
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
  k       = 0L,
  r       = 0,
  offsets = c(1L,1L,1L,1L,1L,1L)
) {

  # If syn is a GRangePairs, try to infer the contig lengths from  the seqinfo
  # for each GRange object.
  if('GRangePairs' %in% class(syn)){
    a <- CNEr::first(syn)
    b <- CNEr::last(syn)
    if(!all(is.na(GenomeInfoDb::seqlengths(a))))
      tcl <- GenomeInfoDb::seqinfo(a)
    if(!all(is.na(GenomeInfoDb::seqlengths(b))))
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

  d <- do_offsets(d, offsets)

  if(is.character(qcl) && qcl == "") qcl <- NULL
  if(is.character(tcl) && tcl == "") tcl <- NULL

  SearchResult(
    CNEr::GRangePairs(
      first  = .make_GRanges(
        as.character(d$qseqid),
        d$qstart,
        d$qstop,
        seqinfo=qcl
      ),
      second = .make_GRanges(
        as.character(d$tseqid),
        d$tstart,
        d$tstop,
        seqinfo=tcl
      ),
      attr      = as.character(d$attr),
      strand    = d$strand,
      score     = d$score,
      cset      = d$cset,
      l_flag    = d$l_flag,
      r_flag    = d$r_flag,
      inbetween = d$inbetween
    ),
    swap    = swap,
    trans   = trans,
    k       = k,
    r       = r,
    offsets = offsets
  )

}

#' @rdname synder_commands
#' @export
dump <- function(
  syn,
  swap    = FALSE,
  trans   = 'i',
  offsets = c(1L,1L,1L,1L,1L,1L)
) {

  syn <- as_synmap(syn)

  d <- wrapper(
    FUN     = c_dump,
    x       = syn,
    swap    = swap,
    trans   = trans,
    offsets = offsets[1:4]
  )

  d <- do_offsets(d, offsets)

  DumpResult(
    CNEr::GRangePairs(
      first  = .make_GRanges(
        as.character(d$qseqid),
        d$qstart,
        d$qstop,
        seqinfo=GenomeInfoDb::seqinfo(CNEr::first(syn))
      ),
      second = .make_GRanges(
        as.character(d$tseqid),
        d$tstart,
        d$tstop,
        seqinfo=GenomeInfoDb::seqinfo(CNEr::second(syn))
      ),
      strand = d$strand,
      score  = d$score,
      cset   = d$cset
    ),
    swap    = swap,
    trans   = trans,
    offsets = offsets[1:4]
  )
}

