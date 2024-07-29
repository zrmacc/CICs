#' Generate Competing Risks Data
#' 
#' Assumes exponential cause-specific hazards and independent exponential
#' censoring.
#' 
#' @param n Sample size.
#' @param event_rate Event rate.
#' @param death_rate Death rate.
#' @param censor_rate Censoring rate.
#' @param tau Optional truncation time.
#' @export 
GenData <- function(
  n,
  event_rate,
  death_rate,
  censor_rate = 0,
  tau = NULL
) {
  
  # Overall hazard.
  overall_hazard <- event_rate + death_rate
  
  # Arrivals.
  event_times <- stats::rexp(n = n, rate = overall_hazard)
  
  # Event type.
  pi <- c(event_rate, death_rate) / overall_hazard
  event_status <- stats::rmultinom(n = n, size = 1, prob = pi)
  event_status <- apply(event_status, 2, which.max)
  
  # Censoring.
  if (censor_rate > 0) {
    censor <- stats::rexp(n = n, rate = censor_rate)
    time <- pmin(event_times, censor)
    status <- (event_times <= censor) * event_status
  }  else {
    time <- event_times
    status <- event_status
  }
  
  # Truncation.
  if (!is.null(tau)) {
    tmp_time <- time  
    time <- pmin(time, tau)
    status[tmp_time > tau] <- 0
  }
  
  # Observed data.
  out <- data.frame(
    time = time,
    status = status
  )
  return(out)
  
}


# -----------------------------------------------------------------------------

#' Generate Two-Sample Data
#' 
#' Assumes exponential cause-specific hazards and independent exponential
#' censoring.
#' 
#' @param n1 Sample size in arm 1.
#' @param n0 Sample size in arm 0.
#' @param censor_rate1 Censoring rate in arm 1.
#' @param censor_rate0 Censoring rate in arm 0.
#' @param event_rate1 Event rate in arm 1.
#' @param event_rate0 Event rate in arm 0.
#' @param death_rate1 Death rate in arm 1.
#' @param death_rate0 Death rate in arm 0.
#' @param tau Optional truncation time needed if censoring rate is set to 0.
#' @return Data.frame with observation 'time', event 'status, and 
#'   treatment 'arm'.
#' @export
GenTwoSampleData <- function(
  n1,
  n0,
  censor_rate1 = 0.25,
  censor_rate0 = 0.25,
  event_rate1 = 0.50,
  event_rate0 = 0.50,
  death_rate1 = 0.25,
  death_rate0 = 0.25,
  tau = NULL
) {
  
  # Generate data for arm 1.
  data1 <- GenData(n1, censor_rate1, event_rate1, death_rate1, tau)
  data1$arm <- 1
  
  # Generate data for arm 0.
  data0 <- GenData(n0, censor_rate0, event_rate0, death_rate0, tau)
  data0$arm <- 0
  
  # Output.
  data <- rbind(data1, data0)
  return(data)
}