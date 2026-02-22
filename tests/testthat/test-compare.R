library(testthat)

test_that("CompareCICs returns compCICs object with expected slots.", {
  d <- GenTwoSampleData(n1 = 40, n0 = 40, tau = 5)
  withr::local_seed(30)
  out <- CompareCICs(
    time = d$time, status = d$status, arm = d$arm,
    sum_stat = "AUC", param = 5, reps = 100
  )
  expect_s4_class(out, "compCICs")
  expect_equal(slotNames(out), c("cis", "cics", "pvals", "reps", "marg_stats", "stratum_stats"))
  expect_s3_class(out@cis, "data.frame")
  expect_s3_class(out@cics, "data.frame")
  expect_s3_class(out@pvals, "data.frame")
  expect_type(out@reps, "list")
  expect_true(all(c("boot", "perm") %in% names(out@reps)))
})

test_that("CompareCICs p-values are in [0, 1].", {
  d <- GenTwoSampleData(n1 = 30, n0 = 30, tau = 5)
  withr::local_seed(31)
  out <- CompareCICs(time = d$time, status = d$status, arm = d$arm, reps = 80)
  expect_true(all(out@pvals$perm_p >= 0 & out@pvals$perm_p <= 1))
  expect_true(out@pvals$boot_p[1] >= 0 && out@pvals$boot_p[1] <= 1)
})

test_that("CompareCICs works with strata.", {
  d <- GenTwoSampleData(n1 = 20, n0 = 20, tau = 5)
  d$strata <- rep(1:2, length.out = nrow(d))
  withr::local_seed(32)
  out <- CompareCICs(
    time = d$time, status = d$status, arm = d$arm, strata = d$strata,
    sum_stat = "Quantile", param = 0.5, reps = 60
  )
  expect_s4_class(out, "compCICs")
  expect_true("strata" %in% names(out@stratum_stats))
})

test_that("CompareCICs works for sum_stat AOC, AUC, Quantile, Rate.", {
  d <- GenTwoSampleData(n1 = 25, n0 = 25, tau = 5)
  withr::local_seed(33)
  for (stat in c("AOC", "AUC", "Quantile", "Rate")) {
    param <- if (stat == "Quantile") 0.5 else 5
    out <- CompareCICs(
      time = d$time, status = d$status, arm = d$arm,
      sum_stat = stat, param = param, reps = 50
    )
    expect_s4_class(out, "compCICs")
  }
})

test_that("CompareCICs rejects invalid sum_stat.", {
  d <- GenTwoSampleData(n1 = 10, n0 = 10, tau = 5)
  expect_error(
    CompareCICs(time = d$time, status = d$status, arm = d$arm, sum_stat = "invalid"),
    "Sum stat not available"
  )
})

test_that("print and show work for compCICs.", {
  d <- GenTwoSampleData(n1 = 15, n0 = 15, tau = 5)
  withr::local_seed(34)
  out <- CompareCICs(time = d$time, status = d$status, arm = d$arm, reps = 30)
  expect_s4_class(out, "compCICs")
  expect_output(print(out), "Marginal Stats")
  expect_output(show(out), "Marginal Stats")
})
