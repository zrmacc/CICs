# Purpose: Functions to perform bootstrap inference on the difference in event rates
# between two cumulative incidence curves.

#' Find Event Rate at Time tau.
#'
#' @param times CIC times.
#' @param probs CIC probabilities.
#' @param tau Event time.
#' @importFrom stats stepfun

EventRate <- function(times, probs, tau) {
  g <- stepfun(x = times, y = c(0, probs), right = TRUE)
  return(g(tau))
}


#' Calculate Difference and Ratio in Event Rates. 
#' 
#' @param data Data.frame containing `time`, `status`, `arm`.
#' @param tau Truncation time.
#' @param return_rates Return the event rates?
#' @importFrom cmprsk cuminc
#' @return If `return_rates = TRUE`, list containing:
#' \itemize{
#'   \item `rates` for each arm. 
#'   \item `stats`, including the difference and ratio of areas.
#' }
#'  If `return_rates = FALSE`, only `stats` is returned. 

EventRate.Stats <- function(data, tau, return_rates = FALSE) {
  
  # Fit competing risks model.
  fit <- cuminc(ftime = data$time, fstatus = data$status, group = data$arm)
  
  # Cumulative incidence curves.
  r0 <- EventRate(fit$`0 1`$time, fit$`0 1`$est, tau = tau)
  r1 <- EventRate(fit$`1 1`$time, fit$`1 1`$est, tau = tau)
  rates <- c('r0' = r0, 'r1' = r1)
  
  # Difference and ratio
  diff <- r1 - r0
  ratio <- r1 / r0
  stats <- c('diff' = diff, 'ratio' = ratio)
  
  # Output
  if(return_rates){
    out <- list('rates' = rates, 'stats' = stats)
  } else {
    out <- stats
  }
  return(out)
}


#' Bootstrap Inference for the Event Rate Difference.
#'
#' Bootstrap inference for the difference and ratio of event rates
#' at time tau. 
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param tau Event time.
#' @param reps Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom cmprsk cuminc
#' @import stats
#' @export
#' @return Data.frame containing these columns:
#' \describe{
#'   \item{Time}{Truncation time.}
#'   \item{Arm0}{Event rate for arm 0.}
#'   \item{Arm1}{Event rate for arm 1.}
#'   \item{Contrast}{'A1-A0' for the difference and 'A1/A0' for the ratio.}
#'   \item{Estimate}{The estimated difference and ratio of AUCs.}
#'   \item{L}{Confidence lower bound.}
#'   \item{U}{Confidence upper bound.}
#'   \item{P}{P-value.}
#' }

CompareEventRates <- function(
  time, 
  status, 
  arm, 
  tau, 
  reps = 2000, 
  alpha = 0.05
) {
  
  # Data.
  data <- data.frame(time, status, arm)
  n <- nrow(data)
  
  # Split data.
  data1 <- data[data$arm == 1, ]
  data0 <- data[data$arm == 0, ]
  
  # Observed difference.
  obs <- EventRate.Stats(data = data, tau = tau, return_rates = TRUE)
  obs_rates <- obs$rates
  obs_stats <- obs$stats
  
  # Bootstrap function.
  aux <- function(b) {
    
    # Bootstrap data sets.
    boot1 <- BootData(data1)
    boot0 <- BootData(data0)
    boot <- rbind(boot1, boot0)
    
    # Bootstrap statistics.
    boot_stats <- EventRate.Stats(data = boot, tau = tau)
    
    # Permuted data.
    boot$arm <- boot$arm[sample(n, n, replace = FALSE)]
    perm_stats <- EventRate.Stats(data = boot, tau = tau)
    
    # Output
    out <- boot_stats
    out[3] <- 1 * abs(perm_stats[1]) >= abs(obs_stats[1])
    out[4] <- 1 * abs(log(perm_stats[2])) >= abs(log(obs_stats[2]))
    return(out)
  }
  
  # Bootstrapping. 
  sim <- lapply(seq(1:reps), aux)
  sim <- do.call(rbind, sim)
  
  # Confidence interval
  alpha2 <- alpha / 2
  lower <- apply(
    X = sim[, 1:2], 
    MARGIN = 2, 
    FUN = function(x) {as.numeric(quantile(x, alpha2, na.rm = TRUE))}
  )
  upper <- apply(
    X = sim[, 1:2], 
    MARGIN = 2, 
    FUN = function(x){as.numeric(quantile(x, 1 - alpha2, na.rm = TRUE))}
  )
  
  # P-values
  pval <- apply(
    X = rbind(sim[, 3:4], c(1, 1)),
    MARGIN = 2,
    FUN = function(x){as.numeric(mean(x, na.rm = TRUE))}
  )
  
  # Output
  out <- data.frame(
    'Time' = c(tau, tau), 
    'Arm0' = rep(x = obs_rates[1], times = 2),
    'Arm1' = rep(x = obs_rates[2], times = 2)
  )
  out$Contrast <- c('A1-A0', 'A1/A0')
  out$Estimate <- obs_stats
  out$L <- lower
  out$U <- upper
  out$P <- pval
  return(out)
}

