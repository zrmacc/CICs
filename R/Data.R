#' Generate Competing Risks Data
#' 
#' Assumes exponential cause-specific hazards and independent exponential
#' censoring.
#' 
#' @param n Sample size.
#' @param censor_rate Censoring rate.
#' @param event_rate Event rate.
#' @param death_rate Death rate.
#' @param tau Optional truncation time needed if censoring rate is set to 0.
#' @return Data.frame containing the observation 'time' and event 'status'.
#' @importFrom stats rexp rmultinom
#' @export 

GenData <- function(
  n,
  censor_rate,
  event_rate,
  death_rate,
  tau = NULL
) {
  
  # Censoring times.
  if (censor_rate > 0) {
    censor <- rexp(n = n, rate = censor_rate)
  } else {
    if (is.null(tau)) {stop("Tau is required if censoring is absent.")}
    censor <- rep(tau, n)
  }

  
  # Overall hazard.
  overall_hazard <- event_rate + death_rate
  
  # Arrivals.
  event_times <- rexp(n = n, rate = overall_hazard)
  
  # Event type.
  pi <- c(event_rate, death_rate) / overall_hazard
  status <- rmultinom(n = n, size = 1, prob = pi)
  status <- apply(status, 2, which.max)
  
  # Observed data.
  out <- data.frame(
    time = pmin(event_times, censor),
    status = (event_times <= censor) * status
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