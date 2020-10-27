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


# -----------------------------------------------------------------------------

#' Find Summary Statistic
#' 
#' Calculate a summary statistic from a cumulative incidence curve. 
#'
#' @param data Data.frame containing `time` and `status`.
#' @param sum_stat Summary statistic, from among 'AUC', 'Quantile', 'Rate'.
#' @param param Either the truncation time, if `sum_stat` is 'AUC' or 'Rate', or 
#'   the quantile probability, if `sum_stat` is 'Quantile'.
#' @return Numeric summary statistic.
#' @export 

FindStat <- function(
  data,
  sum_stat,
  param
) {
  
  # Construct cumulative incidence curve.
  tab <- TabulateCIC(time = data$time, status = data$status)
  times <- tab$Time[tab$Status == 1]
  probs <- tab$CIC[tab$Status == 1]
  
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
#' @param data0 Data.frame containing `time` and `status` for arm 0.
#' @param data1 Data.frame containing `time` and `status` for arm 1.
#' @param sum_stat Summary statistic, from among 'AOC', 'AUC', 'Quantile', 'Rate'.
#' @param param Truncation time, if `sum_stat` is 'AOC', 'AUC', 'Rate';
#'   the quantile probability, if `sum_stat` is 'Quantile'.
#' @param return_strata Return per_stratum stats and CICs?
#' @return If `return_per_arm`, list containing:
#' \itemize{
#'   \item `contrasts`, difference and ratio of summary statistics. 
#'   \item `curves`, list of tabulated cumulative incidence curves.
#'   \item `marg`, marginal summary statistics.
#'   \item `weights`, per-stratum summary statistics.
#' }

SumStats <- function(
  data0, 
  data1, 
  sum_stat = 'AUC',
  param,
  return_strata = FALSE
) {
  
  # Partition by strata.
  data1_strata <- split(data1, data1$strata, drop = TRUE)
  data0_strata <- split(data0, data0$strata, drop = TRUE)
  
  # Stratum sizes.
  data1_sizes <- sapply(data1_strata, nrow)
  data0_sizes <- sapply(data0_strata, nrow)
  sizes <- data1_sizes + data0_sizes
  weights <- sizes / sum(sizes)
  
  # Stats.
  stratum_stat <- function(x) {FindStat(data = x, sum_stat = sum_stat, param = param)}
  data1_stats <- sapply(data1_strata, stratum_stat)
  data0_stats <- sapply(data0_strata, stratum_stat)
  
  # Final stats.
  stat1 <- sum(data1_stats * weights)
  stat0 <- sum(data0_stats * weights)
  
  # Contrasts.
  diff <- stat1 - stat0
  ratio <- stat1 / stat0
  contrasts <- c('diff' = diff, 'ratio' = ratio)
  
  # Output.
  if (return_strata) {
    
    # Cumulative incidence curves.
    stratum_curve <- function(x) {TabulateCIC(time = x$time, status = x$status)}
    strata_names <- names(data1_sizes)
    curves1 <- lapply(data1_strata, stratum_curve)
    curves0 <- lapply(data0_strata, stratum_curve)
    curves <- lapply(strata_names, function(name) {
      curve1 <- curves1[[name]]
      curve1$Arm <- 1
      curve0 <- curves0[[name]]
      curve0$Arm <- 0
      out <- rbind(curve1, curve0)
      out$Stratum <- name
      return(out)
    })
    curves <- do.call(rbind, curves)
    
    # Marginal summary statistics.
    marg_df <- data.frame(
      'Arm' = c(0, 1),
      'N' = c(sum(data0_sizes), sum(data1_sizes)),
      'Stat' = sum_stat,
      'Est' = c(stat0, stat1)
    )
    
    # Per-stratum summary statistics.
    weights_df <- data.frame(
      'Stratum' = names(data1_sizes),
      'Weight' = weights,
      'Stat' = sum_stat,
      'N0' = data0_sizes,
      'Stat0' = data0_stats,
      'N1' = data1_sizes,
      'Stat1' = data1_stats
    )
    weights_df$Diff <- weights_df$Stat1 - weights_df$Stat0
    weights_df$Ratio <- weights_df$Stat1 / weights_df$Stat0
    
    out <- list(
      'contrasts' = contrasts,
      'curves' = curves,
      'marg' = marg_df,
      'weights' = weights_df
    )
  } else {
    out <- contrasts
  }
  return(out)
}
