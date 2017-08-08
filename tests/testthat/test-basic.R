context("basic")

data(toy)

test_that(
  "casts all have the correct type",
  {
    expect_equal(class(as_synmap(toy$synmap))[[1]], 'Synmap')
    expect_equal(class(as.data.frame(as_synmap(toy$synmap)))[[1]], 'data.frame')
  }
)

test_that(
  "as_* and as.data.frame are symmetric",
  {
    expect_equal(
      toy$synmap,
      as.data.frame(as_synmap(toy$synmap))
    )
    expect_equal(
      toy$qgff[, c(1,4,5,7,9)],
      as.data.frame(as_gff(toy$qgff))[, c(1,4,5,7,9)]
    )
  }
)
