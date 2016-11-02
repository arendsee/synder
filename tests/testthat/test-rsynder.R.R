context("rsynder.R")

o1  <- synder::dump('overlap/map-1.syn',  trans='p')
o2  <- synder::dump('overlap/map-2.syn',  trans='p')
o3  <- synder::dump('overlap/map-3.syn',  trans='p')
o4  <- synder::dump('overlap/map-4.syn',  trans='p')
o5  <- synder::dump('overlap/map-5.syn',  trans='p')
o6  <- synder::dump('overlap/map-6.syn',  trans='p')
o7  <- synder::dump('overlap/map-7.syn',  trans='p')
o8  <- synder::dump('overlap/map-8.syn',  trans='p')
o9  <- synder::dump('overlap/map-9.syn',  trans='p')
o10 <- synder::dump('overlap/map-10.syn', trans='p')
o11 <- synder::dump('overlap/map-11.syn', trans='p')

test_that("Assert overlapping regions are correctly merged", {
  expect_equal(o1[[2]], c(100,400,600))
  expect_equal(o1[[5]], c(1100,1500,1700))
  expect_equal(o1[[7]], c(101,101,101))

  expect_equal(o2[[2]], c(100,400,600))
  expect_equal(o2[[5]], c(1100,1500,1700))
  expect_equal(o2[[7]], c(101,101,101))

  expect_equal(o3[[2]], c(100,400,600))
  expect_equal(o3[[5]], c(1100,1500,1700))
  expect_equal(o3[[7]], c(101,101,101))

  expect_equal(o4[[2]], c(100,400,600))
  expect_equal(o4[[5]], c(1100,1500,1700))
  expect_equal(o4[[7]], c(101,101,101))

  expect_equal(o5[[2]], c(100,400,600))
  expect_equal(o5[[5]], c(1100,1500,1700))
  expect_equal(round(o5[[7]]), c(101,141,101))

  expect_equal(o6[[2]], c(100,400,600))
  expect_equal(o6[[5]], c(1100,1500,1700))
  expect_equal(round(o6[[7]]), c(101,141,101))

  expect_equal(o7[[2]], c(100,400))
  expect_equal(o7[[5]], c(1100,1500))
  expect_equal(round(o7[[7]]), c(101,91))

  expect_equal(o8[[2]], c(100,1100))
  expect_equal(o8[[5]], c(1100,1500))
  expect_equal(round(o8[[7]]), c(101,700))

  expect_equal(o9[[2]], c(100,400,600,600))
  expect_equal(o9[[5]], c(1100,1500,600,600))
  expect_equal(o9[[7]], c(101,101,101,101))

  expect_equal(o10[[2]], c(100,100,300))
  expect_equal(o10[[5]], c(300,100,400))
  expect_equal(round(o10[[7]]), c(101,101,201))

  expect_equal(o11[[2]], c(100))
  expect_equal(o11[[5]], c(100))
  expect_equal(round(o11[[7]]), c(251))
})
