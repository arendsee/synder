**This program is under development, 

# Synder

    Map query intervals to target search spaces using a synteny map

# Installation

``` bash
make
make install
```

This will install the program into /usr/local. To install elsewhere, run (for
example)

``` bash
make
make install PREFIX=$HOME
```

To uninstall Synder, run

``` bash
make uninstall
```

Or

``` bash
make uninstall PREFIX=<PATH>
```

# Getting help

To get a usage statement and descriptions of commands, type

`synder -h`

# Contiguous Set Use

./synder -i GFF -s DB -c contig

output is like -c map with additon of flag
SEQNAME TARGETNAME SEARCH_INTERVAL_START SEARCH_INTERVAL_STOP FLAG

Flag is to keep track of edge dependability:

 0 -> the search interval is bound on both ends (is reliable)
 1 -> the start edge is unbounded
 2 -> the stop edge is unbound
 3 -> both edges are unbound, but there are internal overlaps
 4 -> query is to the left of a contig, no overlap
 5 -> query is to the right of a contig, no overlap

# Definitions

 * query interval - an interval on the query genome

 * block - a syntenic pair of intervals linking a query and target genome

 * synteny map - a set of blocks for a pair of genomes

 * interval adjacency - two intervals are adjacent if they are on the same
   chromosome (or scaffold) and no other interval is fully contained between
   them

 * block adjacency - two blocks are adjacent if 1) the intervals are adjacent on
   both the query and target sides and 2) both have the same sign.

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

Note: The image above is out of date. We now guess a search interval for cases
E and F.

# TODO list for version 1

 - [x] Implement the query context to contiguous function
 - [x] Implement the contiguous set to search interval function
 - [x] Clean up the CLI user interface
 - [x] Specify the synteny class (A-F) of the output search intervals 

# TODO list for version 2

 - [ ] Distinguish between tandem duplicates. Currently query-side duplicates
   will map to the same search spaces. But, at least for tandem duplicates, we
   may be able assign unique search spaces to each. But is this worth doing?
 - [ ] Identify target side insertions. This may not be necessary since the
   Fagin pipeline will find these.
 - [ ] Account for uncertainty in syntenic blocks. This may be hard.
 - [x] Possibly try to resolve cases E and F (perhaps as an optional feature)
