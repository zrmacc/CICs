library(testthat)

test_that("CICurve returns stepfun and evaluates.", {
  d <- data.frame(time = c(1, 2, 3), status = c(1, 1, 0))
  g <- CICurve(d)
  expect_true(is.function(g))
  expect_equal(g(0), 0)
  expect_true(g(2.5) >= 0 && g(2.5) <= 1)
})

test_that("CICurve respects custom status_name and time_name.", {
  d <- data.frame(t = c(1, 2, 3), s = c(1, 1, 0))
  g <- CICurve(d, status_name = "s", time_name = "t")
  expect_equal(g(0), 0)
  expect_true(is.numeric(g(2)))
})

test_that("SECurve returns stepfun.", {
  d <- data.frame(time = c(1, 2, 3), status = c(1, 1, 0))
  g <- SECurve(d)
  expect_true(is.function(g))
  expect_true(is.numeric(g(1)))
})

test_that("NARCurve returns stepfun with n at time 0.", {
  d <- data.frame(time = c(1, 2, 3), status = c(1, 1, 0))
  g <- NARCurve(d)
  expect_true(is.function(g))
  expect_equal(g(0), 3)
})

test_that("OneSampleCICDF returns data.frame with time and prob.", {
  d <- data.frame(time = c(1, 2, 3, 4), status = c(1, 1, 1, 0))
  df <- OneSampleCICDF(d, eval_points = 10, tau = 4)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("time", "prob") %in% names(df)))
  expect_equal(nrow(df), 10)
  expect_true(all(df$prob >= 0 & df$prob <= 1))
})

test_that("TwoSampleCICDF returns both arms.", {
  d <- GenTwoSampleData(n1 = 20, n0 = 20, tau = 5)
  withr::local_seed(40)
  df <- TwoSampleCICDF(d, eval_points = 5, tau = 5)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("time", "prob", "arm") %in% names(df)))
  expect_equal(length(unique(df$arm)), 2L)
  arms_num <- as.numeric(as.character(unique(df$arm)))
  expect_true(all(arms_num %in% c(0, 1)))
})

test_that("OneSampleNARDF returns time and nar.", {
  d <- data.frame(time = c(1, 2, 3), status = c(1, 1, 0))
  df <- OneSampleNARDF(d, x_breaks = c(0, 1, 2, 3))
  expect_true(all(c("time", "nar") %in% names(df)))
  expect_equal(df$time, c(0, 1, 2, 3))
})

test_that("TwoSampleNARDF returns nar_ctrl and nar_trt.", {
  d <- GenTwoSampleData(n1 = 15, n0 = 15, tau = 5)
  df <- TwoSampleNARDF(d, x_breaks = c(0, 2, 4))
  expect_true(all(c("time", "nar_ctrl", "nar_trt") %in% names(df)))
})

test_that("PlotOneSampleCIC returns ggplot.", {
  d <- data.frame(time = c(1, 2, 3), status = c(1, 1, 0))
  p <- PlotOneSampleCIC(d)
  expect_s3_class(p, "ggplot")
})

test_that("PlotCICs returns ggplot.", {
  d <- GenTwoSampleData(n1 = 10, n0 = 10, tau = 5)
  p <- PlotCICs(d)
  expect_s3_class(p, "ggplot")
})

test_that("PlotAUCIC returns ggplot.", {
  d <- GenTwoSampleData(n1 = 10, n0 = 10, tau = 5)
  p <- PlotAUCIC(d, which_arm = 0, tau = 5)
  expect_s3_class(p, "ggplot")
})

test_that("PlotOneSampleNARs returns ggplot.", {
  d <- data.frame(time = c(1, 2, 3), status = c(1, 1, 0))
  p <- PlotOneSampleNARs(d, x_breaks = c(0, 1, 2, 3))
  expect_s3_class(p, "ggplot")
})

test_that("PlotNARs returns ggplot.", {
  d <- GenTwoSampleData(n1 = 10, n0 = 10, tau = 5)
  p <- PlotNARs(d, x_breaks = c(0, 2, 4))
  expect_s3_class(p, "ggplot")
})
