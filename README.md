[![unstable](http://badges.github.io/stability-badges/dist/unstable.svg)](http://github.com/badges/stability-badges)
[![Travis-CI Build Status](https://travis-ci.org/arendsee/synder.svg?branch=master)](https://travis-ci.org/arendsee/synder)
[![Coverage Status](https://img.shields.io/codecov/c/github/arendsee/synder/master.svg)](https://codecov.io/github/arendsee/synder?branch=master)
[![DOI](https://zenodo.org/badge/49688107.svg)](https://zenodo.org/badge/latestdoi/49688107)

# Synder

    Map intervals on one genome to search spaces on another genome using
    a synteny map

## Citation

If you use this tool, please cite:

    Zebulun Arendsee, Andrew Wilkey, Urminder Singh, Jing Li, Manhoi Hur, and Eve Wurtele. "synder: inferring genomic orthologs from synteny maps." BioRxiv (2019).

## Funding

This work is funded by the National Science Foundation grant:

[NSF-IOS 1546858](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1546858)
Orphan Genes: An Untapped Genetic Reservoir of Novel Traits

## Installation

In an R session
``` R
library(devtools)
install_github('arendsee/synder')
```

To build the vignettes in an R shell run

```
library(devtools)
devtools::build_vignettes()
```

If you use RStudio, then there is probably some button for this (GUIs are too
volatile for me to say anything terribly helpful).

## Troubleshooting

If you get an error during install saying: 

```R
** preparing package for lazy loading                                                              
Error : object 'extract' not found whilst loading namespace 'R.utils'
ERROR: lazy loading failed for package ‘synder’                                
* removing ‘/home/<username>/R/x86_64-pc-linux-gnu-library/3.4/synder’           
Error: Command failed (1) 
```

Try removing the packages `R.oo` and `R.utils` and then, in a vanilla setting
without `magrittr` loaded, reinstall them. Then try installing `synder`.

If you get weird errors in the data, for example `data(toy)` and then see a lot
of hash marks in the GFF files. Then start a new session, load `GenomicRanges`
and `Biostrings`. Then try reloading synder.

## Documentation

First read the "intro" vignette

```
library(synder)
vignette("intro", package="synder")
```

Other vignettes can be listed with `vignette(package="synder")`.

Information on specific command is available through normal channels:

```
library(synder)
?synder::search
?synder::anon_search
```
