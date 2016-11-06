[![Travis-CI Build Status](https://travis-ci.org/arendsee/synder.svg?branch=master)](https://travis-ci.org/arendsee/synder)

# Synder

    Map query intervals to target search spaces using a synteny map

# Installation

In an R session
``` bash
library(devtools)
install_github('arendsee/synder')
```

# Documentation

See `?search` for full documentation
 
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
 - [x] write test code for merge scores
 - [ ] write test code for search interval scores
 - [x] test against fagin
 - [x] refactor to c++
 - [x] clean up IO
 - [x]  - replace getopt
 - [x]  - incorporate subcommands
 - [x]  - allow reading of GFF files with string sequence names
 - [ ]  - improve input file type checking, fail on misformatted files
 - [ ]  - extract name from GFF 9th column, i.e `s/.*ID=([^;]+).*/\1/`.
 - [x]  - if we get an argument that is not in the subcommands list, should die
 - [x] directly parse synteny files, no database Bash script
 - [x] implement filter
 - [ ] write tests for filter
 - [x] remove hard-coding of score decay parameter
 - [x] reimplement dump blocks
 - [x] add quiet mode to runtests.sh
 - [ ] add score to filter
 - [ ] implement assembly checking
 - [ ] make Github wiki
 - [ ] make Github pages site
 - [x] reimplement as R package
 - [x] determine when to use factors in outputs
 - [x] fix screwy search output
 - [x] add offset handling
 - [x] add all tests from test\_data
 - [x] debug last broken test cases
 - [ ] add tests for all read\_\* functions
 - [ ] add test for read\_synmap ~= dump
 - [ ] add tests to cover UI
 - [x] refactor R code into more modular files
 - [x] write simple vignette
 - [x] complete documentation for each function
 - [x] update README (with links to vignettes and manual)
 - [ ] write standalone command line script
 - [ ] add more interesting dataset
 - [ ] write generic diagnose and diagnose.\* functions
 - [x] write plot.\* functions
 - [x] write print.\* functions
 - [x] write plot\_context function
 - [ ] write fna to scaflen function
 - [ ] write n-strings function

# Theoretical stuff

I do not know a good way to do these, nor am I certain of their value.

 - [ ] implement contiguous set scoring
 - [ ] write synteny map density functions
