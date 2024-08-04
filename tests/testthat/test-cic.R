library(testthat)

test_that("Check tabulation of censorings, events, and deaths.", {
  
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(0, 1, 2, 0, 1)
  )
  obs <- TabulateEvents(status = data$status, time = data$time)
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
  
  obs <- CalcCIC(status = data$status, time = data$time)
  expect_equal(obs$cic_event, c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))
  
  # Case: censoring but not death.
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(0, 1, 1, 1, 1)
  )
  obs <- CalcCIC(status = data$status, time = data$time)
  expect_equal(obs$cic_event, c(0.0, 0.0, 0.25, 0.50, 0.75, 1.0))
  
  # Case: death but not censoring.
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(2, 1, 1, 1, 1)
  )
  obs <- CalcCIC(status = data$status, time = data$time)
  expect_equal(obs$cic_event, c(0.0, 0.0, 0.2, 0.4, 0.6, 0.8))
  
  # Case: censoring and death.
  data <- data.frame(
    time = c(1, 2, 3, 4, 5),
    status = c(0, 2, 1, 1, 1)
  )
  obs <- CalcCIC(status = data$status, time = data$time)
  expect_equal(obs$cic_event, c(0.0, 0.0, 0.0, 0.25, 0.50, 0.75))
  
})


test_that("Check influence function calculation.", {
  
  withr::local_seed(101)
  n <- 1000
  
  # No censoring.
  data <- GenData(
    n = n,
    event_rate = 0.25,
    death_rate = 0.50,
    tau = 10
  )
  cic <- CalcCIC(status = data$status, time = data$time)
  
  tau <- 1.0
  eval_time <- max(cic$time[cic$time < tau])
  ref_var <- cic$var_cic_event[cic$time == eval_time]
  
  inf <- InfluenceCIC(status = data$status, time = data$time, trunc_time = tau)
  inf_var <- mean(inf^2)
  
  expect_equal(n * ref_var, inf_var, tolerance = 0.005)
  
  # With censoring.
  data <- GenData(
    n = n,
    event_rate = 0.25,
    death_rate = 0.50,
    censor_rate = 0.50,
    tau = 10
  )
  cic <- CalcCIC(status = data$status, time = data$time)
  
  tau <- 1.0
  eval_time <- max(cic$time[cic$time < tau])
  ref_var <- cic$var_cic_event[cic$time == eval_time]
  
  inf <- InfluenceCIC(status = data$status, time = data$time, trunc_time = tau)
  inf_var <- mean(inf^2)
  
  expect_equal(n * ref_var, inf_var, tolerance = 0.005)  
  
})

