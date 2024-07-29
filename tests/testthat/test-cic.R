library(testthat)

test_that("Check tabulation of censorings, events, and deaths.", {
  
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(0, 1, 2, 0, 1)
  )
  obs <- TabulateEvents(time = data$time, status = data$status)
  expect_equal(obs$censor, c(0, 1, 0, 0, 1, 0))
  expect_equal(obs$event, c(0, 0, 1, 0, 0, 1))
  expect_equal(obs$death, c(0, 0, 0, 1, 0, 0))
  
})


test_that("Check calculation of CIC.", {
  
  # Case: no censoring or death.
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(1, 1, 1, 1, 1)
  )
  
  obs <- CalcCIC(time = data$time, status = data$status)
  expect_equal(obs$cic_event, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
  
  # Case: censoring but not death.
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(0, 1, 1, 1, 1)
  )
  obs <- CalcCIC(time = data$time, status = data$status)
  expect_equal(obs$cic_event, c(0.0, 0.0, 0.25, 0.50, 0.75, 1.0))
  
  # Case: death but not censoring.
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(2, 1, 1, 1, 1)
  )
  obs <- CalcCIC(time = data$time, status = data$status)
  expect_equal(obs$cic_event, c(0.0, 0.0, 0.2, 0.4, 0.6, 0.8))
  
  # Case: censoring and death.
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(0, 2, 1, 1, 1)
  )
  obs <- CalcCIC(time = data$time, status = data$status)
  expect_equal(obs$cic_event, c(0.0, 0.0, 0.0, 0.25, 0.50, 0.75))
  
})