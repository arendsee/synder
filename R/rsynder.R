#' @useDynLib synder
#' @importFrom Rcpp sourceCpp
NULL

#' Synder Commands
#'
#' There are several commands
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
#'    \item between flag - 0 if query overlaps a syntenic interval, 1 otherwise
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
#' @section Terms:
#'
#' \itemize{
#'   \item query genome - the genome referenced by the input intervals
#' 
#'   \item target genome - the genome to which input intervals are mapped
#' 
#'   \item synteny map - a set of query/target interval pairs inferred to be syntenic
#' 
#'   \item block - a single pair of intervals from the synteny map
#' 
#'   \item contig - a chromosome, scaffold, contig (`synder` doesn't distinguish
#'   between them)
#' 
#'   \item interval adjacency - two intervals are adjacent if they are on the
#'   same contig and no other interval is fully contained between them
#' 
#'   \item block adjacency - two blocks are adjacent if 1) the intervals are
#'   adjacent on both the query and target sides and 2) both have the same sign
#' 
#'   \item query context - all blocks that overlap or are adjacent to the query
#'   interval
#' 
#'   \item contiguous set - a set where block *i* is adjacent to block *i+1*
#' 
#'   \item search intervals - a set of intervals on the target genome where the
#'   ortholog of the query interval is expected to be (the ortholog search space)
#' }
#'
#' @section Algorithm:
#'
#' 
#' Execution order for the `synder search` command
#' \enumerate{
#'   \item load synteny map
#'   \item transform scores
#'   \item merge doubly-overlapping blocks
#'   \item determine contiguous sets
#'   \item find query-side, contextually-relevant blocks for each input sequence
#'   \item map to overlapping contiguous sets
#'   \item calculate score for input sequence relative to each contiguous set
#'   \item calculate search interval relative to each contiguous set
#' }
#' 
#' \strong{1. Load Synteny Map}
#' 
#' The input synteny map must have the following columns:
#' 
#' \enumerate{
#'   \item{qseqid} query contig id (e.g. Chr1)
#'   \item{qstart} query interval start
#'   \item{qstop}  query interval stop
#'   \item{sseqid} target contig id
#'   \item{sstart} target interval start
#'   \item{sstop}  target interval stop
#'   \item{score}  score of the syntenic match
#'   \item{strand} relative orientation
#' }
#' 
#' \strong{2. Transform scores}
#' 
#' Internally, scores for each block are assumed to be additive. However, there
#' are many possible forms of input scores. Low numbers may represent good matches
#' (e.g. e-value) or high numbers (e.g. bit scores). Scores may be additive
#' (bitscores) or averaged (percent similarity).
#' 
#' For this reason, the user must specify a transform for the score column (column
#' 7 of synteny map). This transform can be one of the following:
#' 
#' \tabular{ll}{
#'   `S' = S`           \tab default, no transformation) \cr
#'   `S' = L * S`       \tab transform from score densities \cr
#'   `S' = L * S / 100` \tab transform from percent identity \cr
#'   `S' = -log(S)`     \tab transform from e-values or p-values
#' }
#' 
#' Where S is input score, S' is the transformed score, and L interval length
#' 
#' \strong{3. Merge doubly-overlapping blocks}
#' 
#' If two blocks overlap on both the query and target side, they should be merged.
#' In a perfect synteny map, this would never occur. But in practice, it is
#' common, and convolutes the formation of contiguous sets.
#' 
#' For example the following syntenic map
#' 
#' \preformatted{
#' que 100 200 tar 100 200 100 +
#' que 300 400 tar 300 400 100 +
#' que 310 390 tar 310 390 100 +
#' que 500 600 tar 500 600 100 +
#' }
#' 
#' would be separated into two contiguous sets with the bounds (100, 400) and
#' (310, 600).
#' 
#' Now if a input intervals at (que, 250, 450) is searched, two search intervals
#' with the same bounds will obtain: (200, 500) and (200, 500). Each will be
#' flagged as unbound (e.g. an interval edge is between contiguous sets).
#' 
#' However, if the second and third blocks are merged, we obtain a single search
#' interval flagged as bound on both ends.
#' 
#' In practice, many repetitive regions of the genome can be massively
#' doubly-overlapping, resulting in many redundant search intervals. This throws
#' of statistics and wastes time searching extra space.
#' 
#' To avoid such problems, \code{synder} merges any blocks that overlap on both the
#' query and target sides.
#' 
#' The bounds of the merged blocks are simply the union.
#' 
#' But the scores also need to be merged. Originally, I used the equation:
#' 
#' \preformatted{
#'  da (la - lo) + db (lb - lo) + lo (da + db) / 2            E1
#' }
#' 
#' Where
#' \itemize{
#'   \item *da* and *db* are the score densities of blocks *a* and *b*
#'   \item *la*, *lb* and *lo* are the lengths of *a*, *b*, and their overlap
#' }
#' 
#' E1 is problematic for two reasons. First, it doesn't really make sense
#' to average the overlap scores. Second, iterative pairwise averaging is
#' assymetric, giving higher weight to the blocks merged later.
#' 
#' I replaced this approach with taking the max of the overlap scores:
#' 
#' \preformatted{
#'  da (la - lo) + db (lb - lo) + lo * max(da, db)            E2
#' }
#' 
#' This fixes the first problem.
#' 
#' The second problem is a bit trickier, since it would require tracking
#' sub-intervals. It would be a pain to implement and probably would have
#' little effect on any real dataset. So I am content with an imperfect
#' solution for now.
#' 
#' 
#' \strong{4. Determine contiguous sets}
#' 
#' Build an interval adjacency matrix for the query context intervals and for the
#' target context intervals. AND them together to get a block adjacency matrix.
#' From this matrix, extract paths of adjacent blocks. Synder will merge any
#' blocks that overlap on both genomes, this ensures there is a unique path
#' through the adjacency matrix.
#' 
#' \strong{5. Find contextual blocks}
#' 
#' Mapping query intervals to target intervals would be trivial for a
#' perfectly linear map, e.g.
#' 
#' \preformatted{
#'                                       [-----------]
#' =====      ============     ==========             ============
#'   |              |               |                       |  
#'   |              |               |                       |  
#' =====      ============     ==========   <--->     ============
#' 
#' Where <---> is the query interval and [----] the search interval
#' }
#' 
#' Synder reduces all blocks in the genome into into contiguous sets, which are
#' non-overlapping sets of intervals where all blocks are adjacent.
#' 
#' \itemize{
#' 
#'   \item interval adjacency - two intervals are adjacent if they are on the same
#'   contig and no other interval is fully contained between them
#' 
#'   \item block adjacency - two blocks are adjacent if 1) the intervals are
#'   adjacent on both the query and target sides and 2) both have the same sign
#' 
#'   \item query context - all blocks that overlap or are adjacent to the query
#'   interval
#' 
#'   \item contiguous set - a set where block *i* is adjacent to block \eqn{i+1}
#' 
#'   \item search intervals - a set of intervals on the target genome where the
#'   ortholog of the query interval is expected to be (the ortholog search space)
#' 
#' }
#' 
#' 
#' \strong{6. Reduce to overlapping sets}
#' 
#' Once all overlapping and flanking query-side blocks are identified.
#' 
#' \strong{7. Calculate scores relative to contiguous sets}
#' 
#' Especially with the high \code{k}, it is important to be able to rank search
#' intervals. Queries that heavily overlap elements of contiguous sets are more
#' reliable than ones that fall inbetween. Likewise, queries in dense contiguous
#' sets, should rank higher than those in sparse ones (all else being equal).
#' 
#' Input scores for syntenic blocks are additive (assuming the user entered the
#' correct transformation).
#' 
#' \preformatted{
#' s = 0
#' for i in [a..b]
#'     if i in any block in c AND i in q
#'         s += d / L
#'     else if i in any block in c
#'         s += d * e^(-r * (dist(q, i)))
#' where
#'     c := the contiguous set
#'     q := query interval
#'     a := contig lower bound
#'     b := contig upper bound
#'     d := query syntenic score
#'     L := query length
#'     dist := function calculating distance
#' }
#' 
#' So all blocks in the contiguous set contribute to the total score. The scores
#' of blocks that do not overlap the query decay exponentially with distance from
#' the query bound.
#' 
#' The score decay rate is controlled by the parameter \code{r}. A value of 0.001,
#' the current default, indicates weight will fall to 0.5 by 1000 bases from the
#' query. k=0 would give equal weight to all elements in the contiguous set, i.e.
#' be more affected by context. A high value, e.g. \eqn{r=100}, would completely
#' ignore context, basing score only on overlapping elements. 
#' 
#' 
#' \strong{8. Calculate search interval relative to contiguous set}
#' 
#' Exactly one search interval is created for each previously selected contiguous
#' set.
#' 
#' The chromosome length file tells synder how long each chromosome is. Then if
#' some input query is closer to the end than anything in the synteny map, the
#' search interval bound can be set to the end of the chromosome, rather than
#' infinity.
#' 
#' \enumerate{
#'   \item Find input interval, *i*, context in the query genome. This context
#'   consists of all query intervals *Sq* that either overlap or flank the input.
#'   \item From *Sq* determine which contiguous sets ,*Sc*, the query overlaps or
#'   borders.
#'   \item For each contiguous set, *c*, in *Sc*, classify each bound of *i* as:
#'   \itemize{
#'     \item anchored - if overlaps a member of *c*
#'     \item bound - if is inbetween two members of *c*
#'     \item unbound - if is more extreme than any member of the set, but is
#'     inbetween two contiguous sets
#'     \item extreme - No entry is found
#'   }
#' }
#' 
#' The snapping rules are detailed below:
#' 
#' \preformatted{
#' KEY:
#' |x--  --y| - bounds of the contiguous block; start==x, stop==y
#' a========b - a syntenic block with start == a and stop == b
#'   <---q    - the query interval, with stop == q (start doesn't matter)
#' a==b--c==d - query bounding blocks in the contiguous set
#' [===]      - a non-bounding block in the same contiguous set
#'  ...  F=== - nearest non-adjacent block ***ON THE TARGET***, F=start
#'      ^     - search interval bound
#' 
#' Possible snap positions (relative to query)
#'   |x...[===]-----a=======b-----c=======d-----[===]...y|  ...  F===
#'                  ^       ^     ^       ^     ^                ^
#' 
#' We always snap to one bound of the relevant block. If the bound falls between
#' the blocks (i1, j1) and (i2, j2), it would be reasonable to set the search
#' interval to (j1 + 1, i2 - 1), however, this results in negative lenghts when
#' the blocks are adjacent. So instead we set such intervals to to (j1, i2).
#' 
#' =====================================================================
#' Query bound is precedes both contiguous set bounds
#' 
#'          |x-----y|
#'       ...a=======b
#'          ^
#'   <---q
#' 
#' -------------------------------------------------------------------------------
#' Query bound falls within a contiguous set element
#' 
#'   |x...a=======b...y|
#'                ^
#'        <--q
#' 
#' -------------------------------------------------------------------------------
#' Query bound falls between elements of the contiguous set
#' 
#'   |x...a=======b-------c=======d...y|
#'                        ^
#'                  <--q
#' 
#' -------------------------------------------------------------------------------
#' Query is beyond the bounds of the contiguous set
#' 
#'   Target side           A=======B      F===   1. map b-\>B   
#'                                 |------^      2. map B->F
#'                                 |             3. set F as SI bound
#'   Query side     |x...--a=======b| ...
#'                                 <---q
#' 
#' -------------------------------------------------------------------------------
#' Query is beyond anything in the synteny map
#' 
#'   Target side           A=======B      THE\_END   
#'                                 |------^      
#'                                 |             1. map b->B   
#'   Query side     |x...--a=======b| ...        2. B is contig extremum
#'                                 <---q         3. set SI bound to contig length
#' 
#' -------------------------------------------------------------------------------
#' }
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
#'
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
#' @param syn synteny map file or object
#' @param a,b start and stop locations
#' @param conid the name of the reference query contig
#' @export
anon_search <- function(syn, a, b, conid, ...){
  stopifnot(length(a) == c(length(b)))
  stopifnot(b >= a)
  N <- length(a)
  gff <- tibble::data_frame(
    conid   = conid,
    source  = '.',
    type    = '.',
    start   = as.integer(a),
    stop    = as.integer(b),
    score   = '.',
    strand  = '.',
    phase   = '.',
    seqname = paste0('seq_', 1:N)
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
  offsets = c(0,0,0,0,0,0)
) {
  d <- wrapper(
    FUN     = c_search,
    x       = syn,
    y       = gff,
    tcl     = tcl,
    qcl     = qcl,
    swap    = swap,
    k       = k,
    r       = r,
    trans   = trans,
    offsets = offsets[1:4]
  )

  d$seqname <- as.character(d$seqname)
  d$qcon    <- as.character(d$qcon)
  d$tcon    <- as.character(d$tcon)
  d$strand  <- as.character(d$strand)

  class(d) <- append('search_result', class(d))
  attributes(d)$swap  <- swap
  attributes(d)$k     <- k
  attributes(d)$r     <- r
  attributes(d)$trans <- trans

  d <- do_offsets(d, offsets)

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
  offsets = c(0,0,0,0,0,0)
) {
  d <- wrapper(
    FUN     = c_filter,
    x       = syn,
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
  offsets = c(0,0,0,0,0,0)
) {
  d <- wrapper(
    FUN     = c_map,
    x       = syn,
    y       = gff,
    swap    = swap,
    offsets = offsets[1:4]
  )

  d$seqname <- as.character(d$seqname)
  d$qcon    <- as.character(d$qcon)
  d$tcon    <- as.character(d$tcon)
  d$strand  <- as.character(d$strand)

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
  offsets = c(0,0,0,0,0,0)
) {
  d <- wrapper(
    FUN     = c_count,
    x       = syn,
    y       = gff,
    swap    = swap,
    offsets = offsets[1:4]
  )

  d$seqname = as.character(d$seqname)

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
  offsets = c(0,0,0,0,0,0)
) {
  d <- wrapper(
    FUN     = c_dump,
    x       = syn,
    swap    = swap,
    trans   = trans,
    offsets = offsets[1:4]
  )
  d$qcon <- as.character(d$qcon)
  d$tcon <- as.character(d$tcon)

  # Assign class and attributes
  class(d) <- append('dump_result', class(d))
  attributes(d)$swap  = swap
  attributes(d)$trans = trans

  d <- do_offsets(d, offsets)

  d
}
