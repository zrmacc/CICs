# Purpose: Functions to perform bootstrap inference on the difference in AUCs
# between two cumulative incidence curves.

#' Find Area Under the Cumulative Incidence Curve.
#'
#' @param times CIC times.
#' @param probs CIC probabilities.
#' @param tau Truncation time.

FindAUC <- function(times, probs, tau) {
  g <- stepfun(
    x = times,
    y = c(0, probs),
    right = TRUE
  )
  area <- integrate(f = g, lower = 0, upper = tau, subdivisions = 2e3)
  return(area$value)
}


#' Find Event Rate at Time tau.
#'
#' @param times CIC times.
#' @param probs CIC probabilities.
#' @param tau Event time.
#' @importFrom stats stepfun

FindRate <- function(times, probs, tau) {
  g <- stepfun(x = times, y = c(0, probs), right = TRUE)
  return(g(tau))
}


#' Find Quantile of Cumulative Incidence Curve.
#'
#' @param times CIC times.
#' @param probs CIC probabilities.
#' @param q Quantile.

FindQuantile <- function(times, probs, q) {
  if (q > max(probs)) {
    out <- NA
  } else {
    probs <- round(probs, digits = 16)
    idx1 <- (probs >= q)
    idx2 <- (probs > q)
    out <- mean(c(min(times[idx1]), min(times[idx2])))
  }
  return(out)
}


#' Find Summary Statistic
#' 
#' Calculate a summary statistic from a cumulative incidence curve. 
#'
#' @param times CIC times.
#' @param probs CIC probabilities.
#' @param sum_stat Summary statistic, from among 'AUC', 'Quantile', 'Rate'.
#' @param param Either the truncation time, if `sum_stat` is 'AUC' or 'Rate', or 
#'   the quantile probability, if `sum_stat` is 'Quantile'.
#' @return Numeric summary statistic.
#' @export 

FindStat <- function(
  times,
  probs,
  sum_stat,
  param
) {
  out <- NULL
  if (sum_stat == 'AUC') {
    out <- FindAUC(times = times, probs = probs, tau = param)
  } else if (sum_stat == 'Quantile') {
    out <- FindQuantile(times = times, probs = probs, q = param)
  } else if(sum_stat == 'Rate') {
    out <- FindRate(times = times, probs = probs, tau = param)
  }
  return(out)
}

# -----------------------------------------------------------------------------

#' Calculate Difference and Ratio in Summary Stats
#' 
#' @param data0 Data.frame containing `time` and `status` for arm 0.
#' @param data1 Data.frame containing `time` and `status` for arm 1.
#' @param sum_stat Summary statistic, from among 'AUC', 'Quantile', 'Rate'.
#' @param param Either the truncation time, if `sum_stat` is 'AUC' or 'Rate', or 
#'   the quantile probability, if `sum_stat` is 'Quantile'.
#' @param return_per_arm Return per-Arm stats and CICs?
#' @importFrom cmprsk cuminc
#' @return If `return_per_arm == TRUE`, list containing:
#' \itemize{
#'   \item `cic`, tabulated cumulative incidence curves.
#'   \item `per_arm`, the per-arm summary statistics. 
#'   \item `stats`, including the difference and ratio of summary statistics. 
#' }
#'  If `return_per_arm == FALSE`, only `stats` is returned. 

SumStats <- function(
  data0, 
  data1, 
  sum_stat = 'AUC',
  param,
  return_per_arm = FALSE
) {
  
  # Fit cumulative incidence curves (CICs). 
  fit0 <- cuminc(ftime = data0$time, fstatus = data0$status)
  fit1 <- cuminc(ftime = data1$time, fstatus = data1$status)
  
  # Tabulate CICs.
  tab0 <- TabulateCIC(fit0)
  tab1 <- TabulateCIC(fit1)
  
  # CICs.
  tab0$Arm <- 0
  tab1$Arm <- 1
  curves <- rbind(tab0, tab1)
  
  # Summary statistics. 
  stat0 <- FindStat(
    times = tab0$Time[tab0$Status == 1],
    probs = tab0$CIC[tab0$Status == 1],
    sum_stat = sum_stat,
    param = param
  )
  stat1 <- FindStat(
    times = tab1$Time[tab1$Status == 1],
    probs = tab1$CIC[tab1$Status == 1],
    sum_stat = sum_stat,
    param = param
  )
  
  # Per arm summary statistics. 
  per_arm_stats <- data.frame(
    'Arm' = c(0, 1),
    'N' = c(nrow(data0), nrow(data1)),
    'Stat' = rep(sum_stat, times = 2),
    'Observed' = c(stat0, stat1)
  )
  
  # Difference and ratio
  diff <- stat1 - stat0
  ratio <- stat1 / stat0
  stats <- c('diff' = diff, 'ratio' = ratio)
  
  # Output
  if(return_per_arm){
    out <- list('cic' = curves, 'per_arm' = per_arm_stats, 'stats' = stats)
  } else {
    out <- stats
  }
  return(out)
}