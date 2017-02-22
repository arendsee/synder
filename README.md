[![Travis-CI Build Status](https://travis-ci.org/arendsee/synder.svg?branch=master)](https://travis-ci.org/arendsee/synder)
[![Coverage Status](https://img.shields.io/codecov/c/github/arendsee/synder/master.svg)](https://codecov.io/github/arendsee/synder?branch=master)


# Synder

    Map intervals on one genome to search spaces on another genome using
    a synteny map

# Installation

In an R session
``` R
library(devtools)
install_github('arendsee/synder', build_vignettes=TRUE)
```

# Documentation

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
