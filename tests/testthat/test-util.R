context("utilities")

test_that(
  ".namedlist",
  {
    x=1;y=2
    expect_equal(.namedlist(x,y), list(x=1,y=2))
  }
)

test_that("overlaps", {
  x <- data.frame(
    foo = c("a", "a", "b"),
    start = c(10, 1000, 500),
    stop = c(20, 1010, 510)
  )
  y <- data.frame(
    bar = c("a", "a", "c"),
    a = c(1, 1005, 500),
    b = c(15, 1015, 510)
  )
  xy <- data.frame(
    foo = c("a", "a"),
    start = c(10, 1000),
    stop = c(20, 1010),
    bar = c("a", "a"),
    a = c(1, 1005),
    b = c(15, 1015)
  )
  expect_equal(.overlaps(x=x, xid="foo", xa="start", xb="stop", y=y, yid="bar", ya="a", yb="b"), xy)
})
