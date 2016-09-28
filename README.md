[![Build Status](https://travis-ci.org/arendsee/synder.svg?branch=master)](https://travis-ci.org/arendsee/synder)

**This program is under development**

# Synder

    Map query intervals to target search spaces using a synteny map

# Installation

``` bash
make
make test
make install
```

This will install the program into /usr/local. To install elsewhere, run (for
example)

``` bash
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
From this matrix, extract paths of adjacent blocks. Synder will merge any
blocks that overlap on both genomes, this ensures there is a unque path through
the adjacency matrix.

# Getting search intervals from contiguous sets

No more than one search interval will be contained within a given contiguous set. 

## Cases

The cases can be described pretty easily with a picture, see figure 1. Note
that any query intervals that occur nearer to the end of a chromosome than any
of the syntenic blocks, will always be considered "syntenically scrambled".

 ![Contiguous set to search interval. Cases E and F are considered syntenically scrambled so no search interval is obtained.](figures/contiguous-set-to-search-interval.pdf)

Note: The image above is out of date. We now guess a search interval for cases
E and F.

# Output

## Search subcommand

One table with the following fields:
 1.  query interval name (e.g. AT1G20300)
 2.  query chromosome name
 3.  query start position
 4.  query stop position
 5.  target chromosome name
 6.  search interval start position on target chromsome
 7.  search interval stop position on target chromsome
 8.  search interval strand ('+' / '-')
 9.  score
 10. contiguous set id
 11. lower flag
     - 0 lower bound is inside a syntenic interval
     - 1 lower bound is between intervals in a contiguous set
     - 2 lower bound does not overlap the contiguous set
     - 3 lower bound is beyond any syntenic interval (near end of scaffold)
 12. upper flag - see lower flag
 13. between flag - 0 if query overlaps a syntenic interval, 1 otherwise

# TODO list

 - [x] Add strand awareness to contiguity rules (so all contiguous sets are elements on the same strand)
 - [x] Determine direction of SI based on strand
 - [x] Snap search space boundaries for to nearest block on target side
 - [x] write tests for 3rd gen cases 1-3
 - [x] debug contiguous set builder for cases 1-3
 - [x] implement contiguous\_set structures
 - [x] write tests for case 4
 - [x] debug contiguous set builder for case 4
 - [x] write tests for dedicated double-overlapper test
 - [x] implement double-overlapper merging
 - [x] add score transforms to positive additive
 - [x] write search interval score to output
 - [x] write contiguous set id to output
 - [ ] write test code for scores
 - [x] test against fagin
 - [x] refactor to c++
 - [x] clean up IO
 - [x]  - replace getopt
 - [x]  - incorporate subcommands
 - [x]  - allow reading of GFF files with string sequence names
 - [ ]  - improve input file type checking, fail on misformatted files
 - [ ]  - extract name from GFF 9th column, i.e `s/.*ID=([^;]+).*/\1/`.
 - [ ]  - if we get an argument that is not in the subcommands list, should die
 - [x] directly parse synteny files, no database Bash script
 - [x] implement filter
 - [ ] write tests for filter
 - [x] reimplement dump blocks
 - [x] add quiet mode to runtests.sh
 - [ ] implement assembly checking
 - [ ] update README documentation
 - [ ] update Doxygen documentation
 - [ ] make Github wiki
 - [ ] make Github pages site

# Theoretical stuff

I do not know a good way to do these, nor am I certain of their value.

 - [ ] implement target side scoring
 - [ ] implement score thresholding
 - [ ] implement contiguous set scoring
