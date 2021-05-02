# Purpose: Functions for bootstrap inference.

#' Construct Confidence Interval
#'
#' @param x Bootstrap replicates.
#' @param alpha Alpha level.
#' @return Numeric vector containing the 'alpha' level, the 'lower' and 'upper'
#'   confidence limits.
#'   
#' @importFrom stats quantile

BootCI <- function (x, alpha = 0.05) {
  
  lower <- quantile(x = x, probs = alpha / 2, na.rm = TRUE)
  upper <- quantile(x = x, probs = 1 - alpha / 2, na.rm = TRUE)
  
  # Output.
  out <- c(
    "lower" = as.numeric(lower),
    "upper" = as.numeric(upper)
  )
  return(out)
}


# -----------------------------------------------------------------------------

#' Bootstrap Inference
#'
#' @param data Data.frame containing: time, status, arm, strata.
#' @param obs_stats Observed contrasts.
#' @param sum_stat Summary statistic.
#' @param param Parameter.
#' @param alpha Type I error.
#' @param reps Simulations replicates.
#' @return List containing the simulation replicates, and the bootstrap CI
#'   and p-value.

BootSim <- function(
  data,
  obs_stats,
  sum_stat,
  param,
  alpha,
  reps
) {
  
  # Partition data.
  data1 <- data[data$arm == 1, ]
  data0 <- data[data$arm == 0, ]
  
  # Bootstrap function.
  Loop <- function(b) {
    
    # Bootstrap data sets.
    boot0 <- StratBoot(data0)
    boot1 <- StratBoot(data1)
    boot <- rbind(boot0, boot1)
    
    # Bootstrap statistics.
    boot_stats <- SumStats(
      data = boot,
      sum_stat = sum_stat,
      param = param,
      return_strata = FALSE
    )
    names(boot_stats) <- paste0("boot_", names(boot_stats))
    
    # Bootstrap p-value indicators.
    # Indicator is 1 if the sign of the difference in areas is opposite that observed.
    is_diff_sign <- sign(boot_stats[1]) != sign(obs_stats[1])
    
    # Results
    out <- c(
      boot_stats,
      is_diff_sign
    )
    return(out)
  }
  
  sim <- lapply(seq_len(reps), Loop)
  sim <- data.frame(do.call(rbind, sim))
  colnames(sim) <- c("boot_diff", "boot_ratio", "is_diff_sign")
  
  # -------------------------------------------------------
  
  # Confidence intervals.
  boot_cis <- lapply(sim[, 1:2], function(x) {BootCI(x, alpha = alpha)})
  boot_cis <- do.call(rbind, boot_cis)
  
  cis <- data.frame(
    contrast = c('A1-A0', 'A1/A0'),
    est = as.numeric(obs_stats),
    boot_cis
  )
  
  # P-value
  pval <- CalcP(sim$is_diff_sign)
  
  # Output.
  out <- list(
    pval = pval,
    cis = cis,
    sim = sim
  )
  return(out)
}

