context("basic")

data(toy)

test_that(
  "casts all have the correct type",
  {
    expect_equal(class(as_synmap(toy$synmap))[[1]], 'Synmap')
    expect_equal(class(as.data.frame.Synmap(as_synmap(toy$synmap)))[[1]], 'data.frame')
  }
)

test_that(
  "as_* and as.data.frame are symmetric",
  {
    expect_equal(
      GenomicRanges::ranges(CNEr::first(toy$synmap)),
      GenomicRanges::ranges(CNEr::first(as_synmap(as.data.frame.Synmap(toy$synmap))))
    )
    expect_equal(
      GenomicRanges::mcols(toy$synmap),
      GenomicRanges::mcols(as_synmap(as.data.frame.Synmap(toy$synmap)))
    )
    expect_equal(
      GenomicRanges::ranges(toy$qgff),
      GenomicRanges::ranges(as_gff(as.data.frame.GFF(toy$qgff)))
    )
    expect_equal(
      GenomicRanges::mcols(toy$qgff),
      GenomicRanges::mcols(as_gff(as.data.frame.GFF(toy$qgff)))
    )
  }
)

test_that(
  "seqinfo is preserved across search and dump operations",
  {
    expect_true( !any( search(toy$synmap, toy$qgff) %>% CNEr::first() %>% is.na ) )
    expect_true( !any( dump(toy$synmap)             %>% CNEr::first() %>% is.na ) )
  }
)

test_that(
  "seqinfo is swapped on swap",
  {
    expect_equal(
      dump(toy$synmap, swap=TRUE) %>% CNEr::first() %>% GenomeInfoDb::seqinfo(),
      toy$synmap %>% CNEr::second() %>% GenomeInfoDb::seqinfo()
    )
    expect_equal(
      search(toy$synmap, toy$tgff, swap=TRUE) %>% CNEr::first() %>% GenomeInfoDb::seqinfo(),
      toy$synmap %>% CNEr::second() %>% GenomeInfoDb::seqinfo()
    )
  }
)
