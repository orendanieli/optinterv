library(optinterv)

set.seed(123)

test_that("check if controls were changed", {
  #create data
  p        <- 10
  j        <- 3
  n        <- 500
  vars     <- matrix(2 * runif(n * p), n, p)
  controls <- 3 * matrix(runif(n * j), n, j)
  outcome  <- rnorm(n) + apply(X, 1, sum) + X[,1]
  #calculate means
  wgt1 <- non_parm(outcome, vars, controls)
  c1   <- apply(controls ,2 ,function(v) weighted.mean(v,wgt1,na.rm = T))
  c0   <- apply(controls, 2, function(v) mean(v,na.rm = T))
  expect_true(sum(abs(c1 - c0)) < .000001)
})

