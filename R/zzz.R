# This file holds infrastructure for corrdinating how vignettes and such are built.

# This snippet is borrowed from stackoverflow:32286610, which in turn was
# borrowed from some unknown source, presumably of Mid-Hadleyn origin.
.onLoad <- function(libname, pkgname) {
  vig_list <- tools:vignetteEngine(package = 'knitr')
  vweave   <- vig_list[['knitr::knitr']][c('weave')][[1]]
  vtangle  <- vig_list[['knitr::knitr']][c('tangle')][[1]]
  tools::vignetteEngine(
      pkgname
    , weave    = vweave
    , tangle   = vtangle
    , pattern  = "[.]Rnd$"
    , package  = pkgname
  )
}
