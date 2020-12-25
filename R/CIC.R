# Purpose: Calculate the cumulative incidence curve.

# -----------------------------------------------------------------------------

#' Tabulate Events
#' 
#' Tabulate the number at risk and the number of events at each unique
#' observation time.
#' 
#' @param time Event time.
#' @param status Status, coded as 1 for the recurrent event, 0 otherwise.
#' @return Data.frame with the censorings, deaths, and events occurring
#'   at each distinct time point. 

TabulateEvents <- function(time, status){
  
  # Initial number at risk.
  init_nar <- length(time)
  
  # Split time by status.
  split_time <- split(x = time, f = status)
  
  # Censoring times.
  censor_times <- sort(split_time[["0"]])
  
  # Event times.
  event_times  <- sort(split_time[["1"]])
  
  # Death times.
  death_times  <- sort(split_time[["2"]])
  
  # Unique observation times.
  out <- data.frame(time = sort(unique(time)))
  n_time <- nrow(out)
  
  # Number of censorings.
  out$censor <- 0
  if (length(censor_times) > 0) {
    out$censor[out$time %in% censor_times] <- table(censor_times)
  }
  
  # Number of deaths.
  out$death <- 0
  if (length(death_times) > 0) {
    out$death[out$time %in% death_times] <- table(death_times)
  }
  
  # Number events.
  out$event <- 0
  if (length(event_times) > 0) {
    out$event[out$time %in% event_times] <- table(event_times)
  }
  
  # Number at risk.
  cum_removed <- cumsum(out$censor + out$event + out$death)
  cum_removed <- c(0, cum_removed[1:(n_time - 1)])
  out$nar <- init_nar - cum_removed
  out$prop_at_risk <- out$nar / init_nar
  
  # Output.
  return(out)
}


# -----------------------------------------------------------------------------

#' Calculate Cumulative Incidence Curve.
#' 
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @export 
#' @return Tabulate cumulative incidence curve. 

CalcCIC <- function(time, status) {
  
  # Tabulate events and numbers at risk.
  n <- length(time)
  out <- TabulateEvents(time, status)
  n_time <- nrow(out)
  out$death_rate <- out$death / out$nar
  out$event_rate <- out$event / out$nar
  out$haz <- (out$death + out$event) / out$nar
  
  # Survival at the beginning of the interval. 
  out$surv_init <- c(1, cumprod(1 - out$haz)[1:(n_time-1)])
  
  # Cumulative incidence curves.
  out$cic_event <- cumsum(out$surv_init * out$event_rate)
  out$cic_death <- cumsum(out$surv_init * out$death_rate)
  
  # Variance calculation.
  # See equation (3) of <https://pubmed.ncbi.nlm.nih.gov/9160487/>.
  ## dN_1 terms:
  var1 <- (out$cic_event)^2 * cumsum(out$event / out$nar^2) +
    cumsum((1 - out$cic_death)^2 * out$event / out$nar^2) - 
    2 * out$cic_event * cumsum((1 - out$cic_death) * out$event / out$nar^2)
  
  ## dN_2 terms
  var2 <- (out$cic_event)^2 * cumsum(out$death / out$nar^2) +
    cumsum((out$cic_event)^2 * out$death / out$nar^2) -
    2 * out$cic_event * cumsum(out$cic_event * out$death / out$nar^2)
  
  # Output.
  out$var_cic_event <- n * (var1 + var2)
  out$se_cic_event <- sqrt(var1 + var2)
  return(out)
}
