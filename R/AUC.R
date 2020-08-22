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


#' Calculate Difference and Ratio in AUCs.
#' 
#' @param data Data.frame containing `time`, `status`, `arm`.
#' @param tau Truncation time.
#' @param return_areas Return the AUCs?
#' @importFrom cmprsk cuminc
#' @return If `return_areas = TRUE`, list containing:
#' \itemize{
#'   \item `areas` under the cumulative incidence curve for each arm.
#'   \item `stats`, including the difference and ratio of areas.
#' }
#'  If `return_areas = FALSE`, only `stats` is returned. 

AUCIC.Stats <- function(data, tau, return_areas = FALSE) {
  
  # Fit competing risks model.
  fit <- cuminc(ftime = data$time, fstatus = data$status, group = data$arm)
  
  # Cumulative incidence curves.
  a0 <- FindAUC(times = fit$`0 1`$time, probs = fit$`0 1`$est, tau = tau)
  a1 <- FindAUC(times = fit$`1 1`$time, probs = fit$`1 1`$est, tau = tau)
  areas <- c('a0' = a0, 'a1' = a1)
  
  # Difference and ratio
  diff <- a1 - a0
  ratio <- a1 / a0
  stats <- c('diff' = diff, 'ratio' = ratio)
  
  # Output
  if(return_areas){
    out <- list('areas' = areas, 'stats' = stats)
  } else {
    out <- stats
  }
  return(out)
}


#' Bootstrap Inference for the Difference in AUCs.
#'
#' Bootstrap inference for the difference and ratio in areas under the
#' cumulative incidence curve up to time tau. 
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param tau Truncation time.
#' @param reps Bootstrap replicates.
#' @param alpha Alpha level.
#' @import stats
#' @export
#' @return Data.frame containing these columns:
#' \describe{
#'   \item{Time}{Truncation time.}
#'   \item{Arm0}{AUC for arm 0.}
#'   \item{Arm1}{AUC for arm 1.}
#'   \item{Contrast}{'A1-A0' for the difference and 'A1/A0' for the ratio.}
#'   \item{Estimate}{The estimated difference and ratio of AUCs.}
#'   \item{L}{Confidence lower bound.}
#'   \item{U}{Confidence upper bound.}
#'   \item{P}{P-value.}
#' }

CompareAUCs <- function(
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
  obs <- AUCIC.Stats(data = data, tau = tau, return_areas = TRUE)
  obs_areas <- obs$areas
  obs_stats <- obs$stats
  
  # Bootstrap function.
  aux <- function(b) {
    
    # Bootstrap data sets.
    boot1 <- BootData(data1)
    boot0 <- BootData(data0)
    boot <- rbind(boot1, boot0)
    
    # Bootstrap statistics.
    boot_stats <- AUCIC.Stats(data = boot, tau = tau)
    
    # Permuted data.
    boot$arm <- boot$arm[sample(n, n, replace = FALSE)]
    perm_stats <- AUCIC.Stats(data = boot, tau = tau)
    
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
    'Arm0' = rep(x = obs_areas[1], times = 2),
    'Arm1' = rep(x = obs_areas[2], times = 2)
  )
  out$Contrast <- c('A1-A0', 'A1/A0')
  out$Estimate <- obs_stats
  out$L <- lower
  out$U <- upper
  out$P <- pval
  return(out)
}