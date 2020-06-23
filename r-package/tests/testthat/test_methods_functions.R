library(optinterv)

set.seed(123)

#create data
p        <- 10
j        <- 3
n        <- 500
vars     <- matrix(2 * runif(n * p), n, p)
controls <- 3 * matrix(runif(n * j), n, j)
outcome  <- rnorm(n) + apply(controls, 1, sum) + vars[,1]

test_that("log weights are linear in log original weights, wage,
          and controls (non-parametric method)", {
  #add weights:
  w <- 1 + rgeom(n, 0.5)
  dat <- as.data.frame(controls)
  dat$lwgt1 <- log(non_parm(outcome, vars, controls, w))
  dat$lwgt <-log(w)
  dat$ly <- log(outcome)
  reg <- lm(lwgt1 ~ ., dat)
  reg_sum <- tryCatch(summary(reg),
                      warning = function(w){suppressWarnings(summary(reg))})
  expect_true(reg_sum$r.squared >.9999)
})

test_that("log weights are linear in log original weights, wage,
          and controls (nearest-neighbors method)", {
  #add weights:
  w <- 1 + rgeom(n, 0.5)
  #calculate mehalanobis distance for all pairs
  char_mat <- cbind(vars, controls)
  vcov <- cov.wt(char_mat, wt = w)$cov
  dist <- distances::distances(char_mat, normalize = vcov)
  dist <- as.matrix(dist)
  #find representative observation
  wgt_mat <- nn_wgt(outcome, vars, controls, w, test = T)
  rep_obs <- which.max(as.vector(apply(wgt_mat, 2, median)))
  #check linearity
  dat <- data.frame(row.names = 1:n)
  dat$lwgt1 <- log(wgt_mat[,rep_obs])
  dat$lwgt <- log(w)
  dat$ly <- log(outcome)
  dat$dst2 <- dist[,rep_obs]^2
  #run a regression only on large values. lower values could be bad for rounding errors
  reg <- lm(lwgt1 ~ ., dat, wgt_mat[,rep_obs] > 10^-12)
  reg_sum <- tryCatch(summary(reg),
                      warning = function(w){suppressWarnings(summary(reg))})
  expect_true(reg_sum$r.squared >.9999)
})

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
  #calculate mehalanobis distance for all pairs
  dist <- distances::distances(char_mat, normalize = vcov)
  dist <- as.matrix(dist)
  #find representative observation
  rep_obs <- which.min(as.vector(apply(dist, 2, median)))
  p <- sum(outcome < outcome[rep_obs]) / n
  wgt1 <- nn(outcome, x, control, lambda = 1, sigma = 100)
  expect_equal(p > 0.5, wgt1[rep_obs] > 1)
})










