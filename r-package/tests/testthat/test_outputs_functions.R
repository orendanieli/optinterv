library(optinterv)

set.seed(223)

#create data
var <- rnorm(1000)
w <- sample(c(1,2,3), 1000, replace = T)

test_that("the mean of all bins is the same as the overall mean", {
  expect_equal(mean(wtd_bin(var, 30, w)), weighted.mean(var, w))
})

test_that("the mean of one quantile is the same as the overall mean", {
  expect_equal(wtd_bin(var, 1, w), weighted.mean(var, w))
})
