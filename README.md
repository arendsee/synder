[![Travis-CI Build Status](https://travis-ci.org/arendsee/synder.svg?branch=master)](https://travis-ci.org/arendsee/synder)
[![Coverage Status](https://img.shields.io/codecov/c/github/arendsee/synder/master.svg)](https://codecov.io/github/arendsee/synder?branch=master)


# Synder

    Map intervals on one genome to search spaces on another genome using
    a synteny map

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
