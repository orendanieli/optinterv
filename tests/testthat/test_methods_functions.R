library(optinterv)

set.seed(123)

#create data
p        <- 10
j        <- 3
n        <- 500
vars     <- matrix(2 * runif(n * p), n, p)
controls <- 3 * matrix(runif(n * j), n, j)
outcome  <- rnorm(n) + apply(controls, 1, sum) + vars[,1]

test_that("controls were changed using the non-parametric method", {
  #calculate means
  wgt1 <- non_parm(outcome, vars, controls)
  c1   <- apply(controls ,2 ,function(v) weighted.mean(v, wgt1))
  c0   <- apply(controls, 2, mean)
  expect_true(sum(abs(c1 - c0)) < .000001)
})


test_that("controls were changed using the nearest-neighbors method", {
  #calculate means
  wgt1 <- nn(outcome, vars, controls)
  c1   <- apply(controls ,2 ,function(v) weighted.mean(v, wgt1))
  c0   <- apply(controls, 2, mean)
  expect_true(sum(abs(c1 - c0)) < .000001)
})

#write test for "reasonable" results

