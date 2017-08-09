## ---- echo=TRUE----------------------------------------------------------
# imports required for this vignette
library(synder)
library(knitr)
library(magrittr)
library(ggplot2)
library(dplyr)
library(readr)

## ------------------------------------------------------------------------
data(intro_1)

## ---- echo=TRUE, results='asis'------------------------------------------
knitr::kable(intro_1$synmap)

