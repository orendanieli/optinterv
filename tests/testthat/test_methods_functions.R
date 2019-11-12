library(optinterv)

set.seed(123)

#create data
p        <- 10
j        <- 3
n        <- 500
vars     <- matrix(2 * runif(n * p), n, p)
controls <- 3 * matrix(runif(n * j), n, j)
outcome  <- rnorm(n) + apply(controls, 1, sum) + vars[,1]

test_that("controls weren't changed using the non-parametric method", {
  #calculate means
  wgt1 <- non_parm(outcome, vars, controls)
  c1   <- apply(controls ,2 ,function(v) weighted.mean(v, wgt1))
  c0   <- apply(controls, 2, mean)
  expect_true(sum(abs(c1 - c0)) < .000001)
})

test_that("non-parametric method returns reasonable weights", {
  quant90 <- quantile(vars[,1], probs = 0.9)
  wgt1 <- non_parm(outcome, vars, controls)
  high_wgt <- wgt1[vars[,1] > quant90]
  expect_true(median(high_wgt) > 1)
})

test_that("nearest-neighbors method returns reasonable weights", {
  quant90 <- quantile(vars[,1], probs = 0.9)
  wgt1 <- nn(outcome, vars, controls)
  high_wgt <- wgt1[vars[,1] > quant90]
  expect_true(median(high_wgt) > 1)
})

test_that("(continuous) controls weren't changed using the nearest-neighbors method", {
  #calculate means
  wgt1 <- nn(outcome, vars, controls)
  c1   <- apply(controls ,2 ,function(v) weighted.mean(v, wgt1))
  c0   <- apply(controls, 2, mean)
  expect_true(sum(abs(c1 - c0)) < .000001)
})

test_that("(catgorical) controls weren't changed using the nearest-neighbors method", {
  #change controls
  controls <- matrix(rbinom(2 * n, size = 1, prob = 0.5), n, 2)
  outcome  <- 5 + rnorm(n) + 2 * controls[,1] - controls[,2] + vars[,1]
  wgt1 <- nn(outcome, vars, controls)
  c1   <- apply(controls ,2 ,function(v) weighted.mean(v, wgt1))
  c0   <- apply(controls, 2, mean)
  expect_true(sum(abs(c1 - c0)) < .01)
})


test_that("controls (with one small group) weren't changed using the nearest-neighbors method", {
  #change controls
  controls <- matrix(rbinom(2 * n, size = 1, prob = 0.95), n, 2)
  outcome  <- 5 + rnorm(n) + 2 * controls[,1] - controls[,2] + vars[,1]
  wgt1 <- nn(outcome, vars, controls)
  c1   <- apply(controls ,2 ,function(v) weighted.mean(v, wgt1))
  c0   <- apply(controls, 2, mean)
  expect_true(sum(abs(c1 - c0)) < .01)
})
