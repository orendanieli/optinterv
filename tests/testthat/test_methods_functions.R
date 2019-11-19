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


test_that("nearest-neighbors returns results which getting closer to the non-parametric
          results as sigma getting bigger", {
  wgt1_np <- non_parm(outcome, vars, controls)
  wgt1_nn_1 <- nn(outcome, vars, controls)
  wgt1_nn_1000 <- nn(outcome, vars, controls, sigma = 1000)
  expect_true(sum(abs(wgt1_nn_1 - wgt1_np)) > sum(abs(wgt1_nn_1000 - wgt1_np)))
})

test_that("(continuous) controls weren't changed using the nearest-neighbors method", {
  wgt1 <- nn(outcome, vars, controls)
  #calculate means
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
  expect_true(sum(abs(c1 - c0)) < .000001)
})


test_that("controls (with one small group) weren't changed using the nearest-neighbors method", {
  #change controls
  controls <- matrix(rbinom(2 * n, size = 1, prob = 0.95), n, 2)
  outcome  <- 5 + rnorm(n) + 2 * controls[,1] - controls[,2] + vars[,1]
  wgt1 <- nn(outcome, vars, controls)
  c1   <- apply(controls ,2 ,function(v) weighted.mean(v, wgt1))
  c0   <- apply(controls, 2, mean)
  expect_true(sum(abs(c1 - c0)) < .000001)
})

test_that("nearest-neighbors method returns reasonable weights", {
  control <- as.matrix(2 * runif(n))
  x <- as.matrix(3 * runif(n))
  outcome  <- 5 + rnorm(n) + control[,1] + x[,1]
  char_mat <- cbind(x, control)
  vcov <- cov.wt(char_mat)$cov
  #find mehalanobis distance for all pairs
  dist <- distances::distances(char_mat, normalize = vcov)
  dist <- as.matrix(dist)
  #find representative observation
  rep_obs <- which.min(as.vector(apply(dist, 2, median)))
  p <- sum(outcome < outcome[rep_obs]) / n
  wgt1 <- nn(outcome, x, control)
  expect_equal(p > 0.5, wgt1[rep_obs] > 1)
})
