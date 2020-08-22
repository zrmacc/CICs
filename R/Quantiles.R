# Purpose: Functions to perform bootstrap inference on the difference in event
# rates between two cumulative incidence curves.

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


#' Calculate Difference and Ratio in Event Rates. 
#' 
#' @param data Data.frame containing `time`, `status`, `arm`.
#' @param q Quantile probability.
#' @param return_quants Return the quantiles?
#' @importFrom cmprsk cuminc
#' @return If `return_quants = TRUE`, list containing:
#' \itemize{
#'   \item `rates` for each arm. 
#'   \item `stats`, including the difference and ratio of areas.
#' }
#'  If `return_quants = FALSE`, only `stats` is returned. 

Quantile.Stats <- function(data, q, return_quants = FALSE) {
  
  # Fit competing risks model.
  fit <- cuminc(ftime = data$time, fstatus = data$status, group = data$arm)
  
  # Cumulative incidence curves.
  q0 <- FindQuantile(fit$`0 1`$time, fit$`0 1`$est, q = q)
  q1 <- FindQuantile(fit$`1 1`$time, fit$`1 1`$est, q = q)
  quants <- c('q0' = q0, 'q1' = q1)
  
  # Difference and ratio
  diff <- q1 - q0
  ratio <- q1 / q0
  stats <- c('diff' = diff, 'ratio' = ratio)
  
  # Output
  if(return_quants){
    out <- list('quants' = quants, 'stats' = stats)
  } else {
    out <- stats
  }
  return(out)
}


#'  Bootstrap Inference for the Difference in Quantiles.
#'
#' Bootstrap inference for the difference and ratio of quantiles.
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param q Quantile probability.
#' @param reps Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom cmprsk cuminc
#' @import stats
#' @export
#' @return Data.frame containing these columns:
#' \describe{
#'   \item{Prob}{Quantile probability.}
#'   \item{Arm0}{Quantile for arm 0.}
#'   \item{Arm1}{Quantile for arm 1.}
#'   \item{Contrast}{'A1-A0' for the difference and 'A1/A0' for the ratio.}
#'   \item{Estimate}{The estimated difference and ratio of AUCs.}
#'   \item{L}{Confidence lower bound.}
#'   \item{U}{Confidence upper bound.}
#'   \item{P}{P-value.}
#' }

CompareQuants <- function(time, status, arm, q = 0.5, reps = 2000, alpha = 0.05) {

  # Data.
  data <- data.frame(time, status, arm)
  n <- nrow(data)
  
  # Split data.
  data1 <- data[data$arm == 1, ]
  data0 <- data[data$arm == 0, ]

  # Observed difference.
  obs <- Quantile.Stats(data = data, q = q, return_quants = TRUE)
  obs_quants <- obs$quants
  obs_stats <- obs$stats

  # Bootstrap function.
  aux <- function(b) {

    # Bootstrap data sets.
    boot1 <- BootData(data1)
    boot0 <- BootData(data0)
    boot <- rbind(boot1, boot0)
    
    # Bootstrap statistics.
    boot_stats <- Quantile.Stats(data = boot, q = q)
    
    # Permuted data.
    boot$arm <- boot$arm[sample(n, n, replace = FALSE)]
    perm_stats <- Quantile.Stats(data = boot, q = q)
    
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
    'Prob' = c(q, q), 
    'Arm0' = rep(x = obs_quants[1], times = 2),
    'Arm1' = rep(x = obs_quants[2], times = 2)
  )
  out$Contrast <- c('A1-A0', 'A1/A0')
  out$Estimate <- obs_stats
  out$L <- lower
  out$U <- upper
  out$P <- pval
  return(out)
}

