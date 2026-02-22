library(testthat)

# Step function: value 0 on [0,1), 0.2 on [1,2), ..., 0.8 on [4,5). Area = 0 + 0.2 + 0.4 + 0.6 + 0.8 = 2.0
test_that("FindAUC matches integrate on step function.", {
  times <- c(0, 1, 2, 3, 4, 5)
  probs <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
  auc <- FindAUC(times = times, probs = probs, tau = 5)
  expect_equal(auc, 2.0, tolerance = 0.01)
})

test_that("FindAOC equals tau minus FindAUC.", {
  times <- c(0, 1, 2, 3)
  probs <- c(0, 0.2, 0.5, 0.8)
  tau <- 3
  expect_equal(FindAOC(times, probs, tau), tau - FindAUC(times, probs, tau))
})

test_that("FindRate returns CIC value at tau.", {
  times <- c(0, 1, 2, 3)
  probs <- c(0, 0.2, 0.5, 0.8)
  expect_equal(FindRate(times, probs, tau = 2), 0.5)
  expect_equal(FindRate(times, probs, tau = 0.5), 0)
})

test_that("FindQuantile returns NA when q > max(probs).", {
  times <- c(0, 1, 2, 3)
  probs <- c(0, 0.2, 0.5, 0.8)
  expect_true(is.na(FindQuantile(times, probs, q = 0.9)))
})

test_that("FindQuantile returns value within time range.", {
  times <- c(0, 1, 2, 3, 4)
  probs <- c(0, 0.2, 0.4, 0.6, 0.8)
  q <- FindQuantile(times, probs, q = 0.5)
  expect_true(q >= min(times) && q <= max(times))
  expect_false(is.na(q))
})

test_that("FindStat returns numeric for each sum_stat.", {
  data <- data.frame(time = c(1, 2, 3, 4, 5), status = c(1, 1, 1, 1, 1))
  expect_type(FindStat(data$time, data$status, "AUC", param = 5), "double")
  expect_type(FindStat(data$time, data$status, "AOC", param = 5), "double")
  expect_type(FindStat(data$time, data$status, "Rate", param = 3), "double")
  expect_type(FindStat(data$time, data$status, "Quantile", param = 0.5), "double")
})

test_that("SumStats without strata returns contrasts vector.", {
  d <- GenTwoSampleData(n1 = 40, n0 = 40, tau = 5)
  d$strata <- 1
  withr::local_seed(10)
  s <- SumStats(data = d, sum_stat = "AUC", param = 5, return_strata = FALSE)
  expect_equal(names(s), c("diff", "ratio"))
  expect_type(s, "double")
})

test_that("SumStats with return_strata returns list with curves and contrasts.", {
  d <- GenTwoSampleData(n1 = 30, n0 = 30, tau = 5)
  d$strata <- 1
  withr::local_seed(11)
  s <- SumStats(data = d, sum_stat = "AUC", param = 5, return_strata = TRUE)
  expect_true(is.list(s))
  expect_true(all(c("contrasts", "curves", "marg_stats", "stratum_stats") %in% names(s)))
  expect_equal(names(s$contrasts), c("diff", "ratio"))
  expect_s3_class(s$curves, "data.frame")
  expect_s3_class(s$marg_stats, "data.frame")
})

test_that("SumStats with strata runs and returns stratum_stats.", {
  d <- GenTwoSampleData(n1 = 40, n0 = 40, tau = 5)
  d$strata <- rep(1:2, length.out = nrow(d))
  s <- SumStats(data = d, sum_stat = "Rate", param = 3, return_strata = TRUE)
  expect_true("strata" %in% names(s$stratum_stats))
})
