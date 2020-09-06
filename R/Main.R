#' Inference on CICs. 
#'
#' Inference for the difference and ratio in summary statistics comparing
#' two CICs. 
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param sum_stat Summary statistic, from among 'AUC', 'Quantile', 'Rate'.
#' @param param Either the truncation time, if `sum_stat` is 'AUC' or 'Rate', or 
#'   the quantile probability, if `sum_stat` is 'Quantile'.
#' @param reps Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom methods new
#' @import stats
#' @export
#' @return Object of class \code{compCICs} containing:
#' \itemize{
#'   \item `@CIs`: Observed difference and ratio in areas with confidence intervals. 
#'   \item `@Curves`: Cumulative incidence curves. 
#'   \item `@Pvals`: Bootstrap and permutation p-values.
#'   \item `@Reps`: Bootstrap and permutation realizations of the test statistics. 
#'   \item `@Stats`: Observed test statistics. 
#' }

CompareCICs <- function(
  time, 
  status, 
  arm, 
  sum_stat = 'AUC',
  param = NULL,
  reps = 2000, 
  alpha = 0.05
) {
  
  # Check sum stat.
  sum_stat_choices <- c('AUC', 'Quantile', 'Rate')
  if (!sum_stat %in% sum_stat_choices) {
    cat('Select sum_stat from:\n')
    cat(sum_stat_choices)
    stop('Sum stat not available.')
  }
  
  # Data.
  data <- data.frame(time, status, arm)
  n <- nrow(data)
  
  # Split data.
  data1 <- data[data$arm == 1, ]
  data0 <- data[data$arm == 0, ]
  
  # Default sum stat param.
  if (is.null(param)) {
    if (sum_stat %in% c('AUC', 'Rate')) {
      param = min(max(data1$time), max(data0$time))
    } else if (sum_stat == 'Quantile') {
      param = 0.5
    }
  }
  
  # Observed difference.
  obs <- SumStats(
    data0 = data0, 
    data1 = data1, 
    sum_stat = sum_stat,
    param = param,
    return_per_arm = TRUE
  )
  obs_stats <- obs$stats
  
  # Bootstrap function.
  aux <- function(b) {
    
    # Bootstrap data sets.
    boot1 <- BootData(data1)
    boot0 <- BootData(data0)
    
    # Bootstrap statistics.
    boot_stats <- SumStats(
      data0 = boot0, 
      data1 = boot1, 
      sum_stat = sum_stat,
      param = param
    )
    names(boot_stats) <- paste0('boot_', names(boot_stats))
    
    # Bootstrap indicator.
    boot_ind <- 1 * (sign(boot_stats[1]) != sign(obs_stats[1]))
    
    # Permuted data.
    perm_data <- PermData(data)
    perm0 <- perm_data[perm_data$arm == 0, ]
    perm1 <- perm_data[perm_data$arm == 1, ]
    perm_stats <- SumStats(
      data0 = perm0, 
      data1 = perm1, 
      sum_stat = sum_stat,
      param = param
    )
    names(perm_stats) <- paste0('perm_', names(perm_stats))
    
    # Permutation indicator.
    perm_sign_diff <- 1 * (sign(perm_stats[1]) != sign(obs_stats[1]))
    perm_ind_diff <- perm_sign_diff * (abs(perm_stats[1]) >= abs(obs_stats[1]))
    perm_ind_ratio <- perm_sign_diff * (abs(log(perm_stats[2])) >= abs(log(obs_stats[2])))
    
    # Output
    out <- boot_stats
    out[3] <- 1 * abs(perm_stats[1]) >= abs(obs_stats[1])
    out[4] <- 1 * abs(log(perm_stats[2])) >= abs(log(obs_stats[2]))
    
    out <- c(
      boot_stats,
      perm_stats,
      boot_ind,
      boot_ind,
      perm_ind_diff,
      perm_ind_ratio
    )
    
    return(out)
  }
  
  # Bootstrapping. 
  sim <- lapply(seq(1:reps), aux)
  sim <- do.call(rbind, sim)
  
  # -------------------------------------------------------
  
  # Confidence intervals.
  cis0 <- data.frame(
    'Contrast' = c('A1-A0', 'A1/A0'),
    'Observed' = as.numeric(obs_stats)
  )
  
  cis1 <- lapply(1:2, function(i) {
    ConstructCI(x = sim[, i], alpha = alpha)
  })
  cis1 <- do.call(rbind, cis1)
  
  cis <- cbind(cis0, cis1)
  
  # -------------------------------------------------------
  
  # P-values
  pval <- data.frame(
    'Contrast' = rep(c('A1-A0', 'A1/A0'), times = 2),
    'Method' = rep(c('Boot', 'Perm'), each = 2)
  )
  
  # P-values
  pval$P <- apply(
    X = rbind(sim[, 5:8], rep(1, 4)),
    MARGIN = 2,
    FUN = mean
  )
  pval$P <- sapply(pval$P, function (p) {min(2 * p, 1)})
  
  # -------------------------------------------------------
  
  # Output
  out <- new(
    Class = 'compCICs',
    CIs = cis,
    Curves = obs$cic,
    Reps = sim[, 1:4],
    Pvals = pval,
    Stats = obs$per_arm
  )
  return(out)
}
