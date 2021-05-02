#' Inference on CICs. 
#'
#' Inference for the difference and ratio in summary statistics comparing
#' two CICs. 
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param strata Optional stratification factor. 
#' @param sum_stat Summary statistic, from among 'AOC', 'AUC', 'Quantile', 'Rate'.
#' @param param Truncation time, if `sum_stat` is 'AOC', 'AUC', or 'Rate';
#'   quantile probability, if `sum_stat` is 'Quantile'.
#' @param reps Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom methods new
#' @export
#' @return Object of class \code{compCICs} containing:
#' \itemize{
#'   \item `@cis`: Observed difference and ratio in areas with confidence intervals. 
#'   \item `@cics`: Cumulative incidence curves. 
#'   \item `@pvals`: Bootstrap and permutation p-values.
#'   \item `@reps`: Bootstrap and permutation realizations of the test statistics. 
#'   \item `@marg_stats`: Marginal test statistics. 
#'   \item `@stratum_stats`: Stratum test statistics.
#' }

CompareCICs <- function(
  time, 
  status, 
  arm, 
  strata = NULL,
  sum_stat = 'AUC',
  param = NULL,
  reps = 2000, 
  alpha = 0.05
) {
  
  # Check the choice of summary statistic.
  sum_stat_choices <- c('AOC', 'AUC', 'Quantile', 'Rate')
  if (!sum_stat %in% sum_stat_choices) {
    cat('Select sum_stat from among:\n')
    cat(sum_stat_choices)
    stop('Sum stat not available.')
  }
  
  # Create single stratum if no strata are provided. 
  if (is.null(strata)){
    strata <- rep(1, length(time))
  }
  
  # Create data.frame.
  data <- data.frame(time, status, arm, strata)
  
  # Default summary statistic parameter.
  if (is.null(param)) {
    if (sum_stat %in% c('AOC', 'AUC', 'Rate')) {
      param = min(
        max(data$time[data$arm == 1]), 
        max(data$time[data$arm == 0])
      )
    } else if (sum_stat == 'Quantile') {
      param = 0.5
    }
  }
  
  # -------------------------------------------------------
  
  # Observed difference.
  obs <- SumStats(
    data = data,
    sum_stat = sum_stat,
    param = param,
    return_strata = TRUE
  )
  obs_stats <- obs$contrasts
  
  # Bootstrap.
  boot <- BootSim(
    data = data,
    obs_stats = obs_stats,
    sum_stat = sum_stat,
    param = param,
    alpha = alpha,
    reps = reps
  )
  
  # Permutation.
  perm <- PermSim(
    data = data,
    obs_stats = obs_stats,
    sum_stat = sum_stat,
    param = param,
    alpha = alpha,
    reps = reps
  )
  
  pvals <- perm$pvals
  pvals$boot_p <- boot$pval
  
  # Output
  out <- new(
    Class = 'compCICs',
    cis = boot$cis,
    cics = obs$curves,
    reps = list("boot" = boot$sim, "perm" = perm$sim),
    pvals = pvals,
    marg_stats = obs$marg_stats,
    stratum_stats = obs$stratum_stats
  )
  return(out)
}
