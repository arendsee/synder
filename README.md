# Manual stub

# Definitions

 * query interval - an interval on the query genome

 * block - a syntenic pair of intervals linking a query and target genome

 * synteny map - a set of blocks for a pair of genomes

 * interval adjacency - two intervals are adjacent if they are on the same
   chromosome (or scaffold) and no other interval is fully contained between
   them

 * block adjacency - two blocks are adjacent if the intervals are adjacent on
   both the query and target sides

 * query context - all blocks that overlap or are adjacent to the query interval

 * contiguous set - a set of blocks where block *i* is adjacent to block *i+1*

 * search intervals - a set of intervals on the target genome where the
   ortholog of the query interval is expected to be (the ortholog search space)

# Getting contiguous sets from query context

Build an interval adjacency matrix for the query context intervals and for the
target context intervals. AND them together to get a block adjacency matrix.
From this matrix, extract paths of adjacent blocks. For biological synteny
maps, these paths *should* be unique.

# Getting search intervals from contiguous sets

No more than one search interval will be contained within a given contiguous set. 

## Cases

The cases can be described pretty easily with a picture, see figure 1. Note
that any query intervals that occur nearer to the end of a chromosome than any
of the syntenic blocks, will always be considered "syntenically scrambled".

 ![Contiguous set to search interval. Cases E and F are considered syntenically scrambled so no search interval is obtained.](figures/contiguous-set-to-search-interval.pdf)

# TODO list for version 1

 - [ ] Implement the query context to contiguous function
 - [ ] Implement the contiguous set to search interval function
 - [ ] Clean up the CLI user interface
 - [ ] Specify the synteny class (A-F) of the output search intervals 

# TODO list for version 2

 - [ ] Distinguish between tandem duplicates. Currently query-side duplicates
   will map to the same search spaces. But, at least for tandem duplicates, we
   may be able assign unique search spaces to each. But is this worth doing?
 - [ ] Identify target side insertions. This may not be necessary since the
   Cadmium pipeline will find these.
 - [ ] Account for uncertainty in syntenic blocks. This may be hard.
 - [ ] Possibly try to resolve cases E and F (perhaps as an optional feature)
