library(testthat)

test_that("GenData returns data.frame with time and status.", {
  withr::local_seed(1)
  d <- GenData(n = 50, event_rate = 0.3, death_rate = 0.2)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 50)
  expect_true(all(c("time", "status") %in% names(d)))
  expect_true(all(d$status %in% c(0, 1, 2)))
  expect_true(all(d$time > 0))
})

test_that("GenData with censoring produces status 0.", {
  withr::local_seed(2)
  d <- GenData(n = 200, event_rate = 0.2, death_rate = 0.2, censor_rate = 0.5, tau = 5)
  expect_true(any(d$status == 0))
  expect_true(all(d$time <= 5))
})

test_that("GenData with tau truncates time and sets status to 0 when truncated.", {
  withr::local_seed(3)
  d <- GenData(n = 100, event_rate = 0.5, death_rate = 0.5, tau = 1)
  expect_true(all(d$time <= 1))
  expect_true(any(d$status == 0))
})

test_that("GenTwoSampleData returns two arms with correct structure.", {
  withr::local_seed(4)
  d <- GenTwoSampleData(n1 = 30, n0 = 40, tau = 10)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 70)
  expect_true(all(c("time", "status", "arm") %in% names(d)))
  expect_equal(sum(d$arm == 1), 30)
  expect_equal(sum(d$arm == 0), 40)
  expect_true(all(d$arm %in% c(0, 1)))
  expect_true(all(d$status %in% c(0, 1, 2)))
})

test_that("GenTwoSampleData with different rates produces different CICs.", {
  withr::local_seed(5)
  d <- GenTwoSampleData(
    n1 = 200, n0 = 200,
    event_rate1 = 0.6, event_rate0 = 0.2,
    death_rate1 = 0.2, death_rate0 = 0.2,
    censor_rate1 = 0.1, censor_rate0 = 0.1,
    tau = 5
  )
  cic1 <- CalcCIC(status = d$status[d$arm == 1], time = d$time[d$arm == 1])
  cic0 <- CalcCIC(status = d$status[d$arm == 0], time = d$time[d$arm == 0])
  expect_gt(max(cic1$cic_event), max(cic0$cic_event))
})
