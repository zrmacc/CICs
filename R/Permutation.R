# Purpose: Functions for permutation inference. 

# -----------------------------------------------------------------------------
# P-value calculation.
# -----------------------------------------------------------------------------

#' Calculate P
#' 
#' Calculates p-value from 1-sided rejection indicator.
#' @param p Vector of 1/0 rejection indicators.
#' @return Numeric p-value.

CalcP <- function(p) {
  out <- 2 * mean(c(1, p), na.rm = TRUE)
  out <- min(out, 1)
  return(out)
}


#' Permutation Inference.
#'
#' @param data Data.frame containing: time, status, arm, strata.
#' @param obs_stats Observed contrasts.
#' @param sum_stat Summary statistic.
#' @param param Parameter.
#' @param alpha Type I error.
#' @param reps Simulations replicates.
#' @return List 

PermSim <- function(
  data,
  obs_stats,
  sum_stat,
  param,
  alpha,
  reps
) {
  
  # Permutation function.
  Loop <- function(b) {
    
    # Permute data.
    perm <- StratPerm(data)
    
    # Permutation statistics
    perm_stats <- SumStats(
      data = perm,
      sum_stat = sum_stat,
      param = param,
      return_strata = FALSE
    )
    names(perm_stats) <- paste0("perm_", names(perm_stats))
    
    # Permutation indicators.
    perm_diff <- perm_stats[1]
    obs_diff <- obs_stats[1]
    perm_ratio <- perm_stats[2]
    obs_ratio <- obs_stats[2]
    
    is_diff_sign <- (sign(perm_diff) != sign(obs_diff))
    is_diff_more_extreme <- abs(perm_diff) >= abs(obs_diff)
    is_ratio_more_extreme <- abs(log(perm_ratio)) >= abs(log(obs_ratio))
    perm_diff_1sided <- is_diff_sign * is_diff_more_extreme
    perm_ratio_1sided <- is_diff_sign * is_ratio_more_extreme
    
    # Results
    out <- c(
      perm_stats,
      "perm_diff_1sided" = perm_diff_1sided,
      "perm_ratio_1sided" = perm_ratio_1sided
    )
    return(out)
  }
  
  sim <- lapply(seq(1:reps), Loop)
  sim <- data.frame(do.call(rbind, sim))
  colnames(sim) <- c(
    "perm_diff",
    "perm_ratio",
    "perm_diff_1sided",
    "perm_ratio_1sided"
  )
  
  # -------------------------------------------------------
  
  # P-values.
  perm_pvals <- apply(sim[, 3:4], 2, CalcP)
  pvals <- data.frame(
    contrast = c('A1-A0', 'A1/A0'),
    est = as.numeric(obs_stats),
    perm_p = as.numeric(perm_pvals)
  )
  out <- list(
    pvals = pvals,
    sim = sim
  )
  return(out)
}