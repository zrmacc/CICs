# Purpose: Functions to perform bootstrap inference on the difference in AUCs
# between two cumulative incidence curves.

#' Find Area Under the Cumulative Incidence Curve.
#'
#' @param times CIC times.
#' @param probs CIC probabilities.
#' @param tau Truncation time.
#' @return Numeric.
FindAUC <- function(times, probs, tau) {
  g <- stats::stepfun(
    x = times,
    y = c(0, probs)
  )
  area <- stats::integrate(f = g, lower = 0, upper = tau, subdivisions = 2e3)
  return(area$value)
}


#' Find Area Over the Cumulative Incidence Curve.
#'
#' @param times CIC times.
#' @param probs CIC probabilities.
#' @param tau Truncation time.
FindAOC <- function(times, probs, tau) {
  area <- tau - FindAUC(times, probs, tau)
  return(area)
}


#' Find Event Rate at Time tau.
#'
#' @param times CIC times.
#' @param probs CIC probabilities.
#' @param tau Event time.
FindRate <- function(times, probs, tau) {
  g <- stats::stepfun(x = times, y = c(0, probs))
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


# -----------------------------------------------------------------------------

#' Find Summary Statistic
#' 
#' Calculate a summary statistic from a cumulative incidence curve. 
#'
#' @param time Observation time.
#' @param status Status indicator.
#' @param sum_stat Summary statistic, from among 'AUC', 'Quantile', 'Rate'.
#' @param param Either the truncation time, if `sum_stat` is 'AUC' or 'Rate', or 
#'   the quantile probability, if `sum_stat` is 'Quantile'.
#' @return Numeric summary statistic.
#' @export 
FindStat <- function(
  time,
  status,
  sum_stat,
  param
) {
  
  # Construct cumulative incidence curve.
  tab <- CalcCIC(status = status, time = time)
  times <- tab$time
  probs <- tab$cic_event
  
  out <- NULL
  if (sum_stat == 'AOC') {
    out <- FindAOC(times = times, probs = probs, tau = param)
  } else if (sum_stat == 'AUC') {
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
#' @param data Data.frame containing 'time', 'status', 'arm', 'strata'.
#' @param sum_stat Summary statistic, from among 'AOC', 'AUC', 'Quantile', 'Rate'.
#' @param param Truncation time, if `sum_stat` is 'AOC', 'AUC', 'Rate';
#'   the quantile probability, if `sum_stat` is 'Quantile'.
#' @param return_strata Return per_stratum stats and CICs?
#' @return If `return_per_arm`, list containing:
#' \itemize{
#'   \item {`contrasts`: difference and ratio of summary statistics.} 
#'   \item {`curves`: list of tabulated cumulative incidence curves.}
#'   \item {`marg`: marginal summary statistics.}
#'   \item {`weights`: per-stratum summary statistics.}
#' }
SumStats <- function(
  data,
  sum_stat = 'AUC',
  param,
  return_strata = FALSE
) {
  
  # Stratum sizes.
  strata <- NULL
  stratum_sizes <- data %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise("n" = dplyr::n(), .groups = "drop") 
  stratum_sizes$weight <- stratum_sizes$n / sum(stratum_sizes$n)
  
  # Stratum stats
  arm <- time <- status <- NULL
  stratum_stats <- data %>%
    dplyr::group_by(arm, strata) %>%
    dplyr::summarise(
      "n" = dplyr::n(),
      "stat" = sum_stat,
      "est" = FindStat(time, status, sum_stat = sum_stat, param = param),
      .groups = "drop"
    ) %>% 
    dplyr::inner_join(
      stratum_sizes[, c("strata", "weight")], 
      by = "strata"
    )
  
  # Marginal stats.
  est <- n <- weight <- NULL
  marg_stats <- stratum_stats %>%
    dplyr::group_by(arm) %>%
    dplyr::summarise(
      "stat" = sum_stat,
      "n" = sum(n),
      "est" = sum(weight * est[!is.na(est)]) / sum(weight[!is.na(est)]),
      .groups = "drop"
    ) 
  
  # Contrasts.
  stat1 <- marg_stats$est[marg_stats$arm == 1]
  stat0 <- marg_stats$est[marg_stats$arm == 0]
  contrasts <- c(
    "diff" = stat1 - stat0, 
    "ratio" = stat1 / stat0 
  )
  
  # Output.
  if (return_strata) {
    
    # Cumulative incidence curves.
    curves <- data %>%
      dplyr::group_by(strata, arm) %>%
      dplyr::summarise(
        CalcCIC(status = status, time = time),
        .groups = "drop"
      ) 
    
    # Per-stratum summary statistics.
    stat <- est0 <- est1 <- NULL
    stratum_stats <- stratum_stats %>% 
      tidyr::pivot_wider(
        id_cols = c(strata, weight, stat),
        names_from = arm,
        values_from = c(n, est),
        names_sep = ""
      ) %>%
      dplyr::mutate(
        diff = est1 - est0,
        ratio = est1 / est0
      )
    
    out <- list(
      contrasts = contrasts,
      curves = curves,
      marg_stats = marg_stats,
      stratum_stats = stratum_stats
    )
  } else {
    out <- contrasts
  }
  return(out)
}
