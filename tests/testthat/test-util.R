context("utilities")

test_that(
  ".namedlist",
  {
    x=1;y=2
    expect_equal(.namedlist(x,y), list(x=1,y=2))
  }
)

# assert rownames are removed
# assert seqids are not factored
test_that(".overlaps", {
  x <- data.frame(
    foo = c("b", "a", "a"),
    start = c(500, 10, 1000),
    stop = c(510, 20, 1010),
    stringsAsFactors=FALSE
  )
  y <- data.frame(
    bar = c("a", "a", "c"),
    a = c(1, 1005, 500),
    b = c(15, 1015, 510),
    stringsAsFactors=FALSE
  )
  xy <- data.frame(
    foo = c("a", "a"),
    start = c(10, 1000),
    stop = c(20, 1010),
    a = c(1, 1005),
    b = c(15, 1015),
    stringsAsFactors=FALSE
  )
  xy_ids <- data.frame(
    foo = c("a", "a"),
    start = c(10, 1000),
    stop = c(20, 1010),
    xid = 2:3,
    a = c(1, 1005),
    b = c(15, 1015),
    yid = 1:2,
    stringsAsFactors=FALSE
  )
  expect_equal(.overlaps(x=x, xid="foo", xa="start", xb="stop",
                         y=y, yid="bar", ya="a",     yb="b",
                         add_id=FALSE), xy)
  expect_equal(.overlaps(x=x, xid="foo", xa="start", xb="stop",
                         y=y, yid="bar", ya="a",     yb="b",
                         add_id=TRUE), xy_ids)
})
