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


#' Generate Bootstrap Confidence Interval for the Difference in Quantiles.
#'
#' Bootstrap confidence interval for the difference in quantiles.
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param q Quantile probability.
#' @param B Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom cmprsk cuminc
#' @import stats
#' @export
#' @return Numeric vector containing the quantile probability `Quantile`, the
#'   estimated quantiles in arms 1 `Arm1` and 0 `Arm0`, the difference `Delta`,
#'   and the lower `L` and upper `U` confidence bounds.

QuantDiffCI <- function(time, status, arm, q = 0.5, B = 2000, alpha = 0.05) {

  # Data.
  data <- data.frame(time, status, arm)

  # Split data.
  data1 <- data[data$arm == 1, ]
  data0 <- data[data$arm == 0, ]

  # Fit competing risks model.
  fit <- cuminc(ftime = data$time, fstatus = data$status, group = data$arm)

  # Cumulative incidence curves.
  q0 <- FindQuantile(fit$`0 1`$time, fit$`0 1`$est, q = q)
  q1 <- FindQuantile(fit$`1 1`$time, fit$`1 1`$est, q = q)

  # Observed difference.
  dobs <- q1 - q0

  # Bootstrap function.
  aux <- function(b) {

    # Bootstrap data sets.
    boot1 <- BootData(data1)
    boot0 <- BootData(data0)
    boot <- rbind(boot1, boot0)

    # Estimate cumulative incidence curves.
    fit_boot <- cuminc(ftime = boot$time, fstatus = boot$status, group = boot$arm)

    # Bootstrap quantiles.
    qb0 <- FindQuantile(
      times = fit_boot$`0 1`$time,
      probs = fit_boot$`0 1`$est,
      q = q
    )
    qb1 <- FindQuantile(
      times = fit_boot$`1 1`$time,
      probs = fit_boot$`1 1`$est,
      q = q
    )

    # Bootstrap difference.
    dboot <- qb1 - qb0
    return(dboot)
  }

  boot <- lapply(seq(1:B), aux)
  boot <- do.call(c, boot)

  # Confidence interval.
  alpha2 <- alpha / 2
  L <- as.numeric(quantile(boot, alpha2, na.rm = TRUE))
  U <- as.numeric(quantile(boot, 1 - alpha2, na.rm = TRUE))

  # Output.
  out <- c("Quantile" = q, "Arm1" = q1, "Arm0" = q0, "Delta" = dobs, "L" = L, "U" = U)
  return(out)
}


#' Generate Permutation P-value for the Difference in Quantiles.
#'
#' Permutation p-value for the difference in quantiles.
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param q Quantile probability.
#' @param B Bootstrap replicates.
#' @importFrom cmprsk cuminc
#' @import stats
#' @export
#' @return Numeric vector containing the quantile probability `Quantile`, the
#'   estimated quantiles in arms 1 `Arm1` and 0 `Arm0`, the difference `Delta`,
#'   and the permutation p-value `P`.

QuantDiffP <- function(time, status, arm, B = 2000, q = 0.5) {

  # Data.
  data <- data.frame(time, status, arm)

  # Split data.
  data1 <- data[data$arm == 1, ]
  data0 <- data[data$arm == 0, ]

  # Fit competing risks model.
  fit <- cuminc(ftime = data$time, fstatus = data$status, group = data$arm)

  # Cumulative incidence curves.
  q0 <- FindQuantile(fit$`0 1`$time, fit$`0 1`$est, q = q)
  q1 <- FindQuantile(fit$`1 1`$time, fit$`1 1`$est, q = q)

  # Observed difference.
  dobs <- q1 - q0

  # Data to permute.
  data_perm <- data
  n <- nrow(data)

  # Function to bootstrap
  aux <- function(b) {

    # Permute treatment
    data_perm$arm <- data$arm[sample(n, n, replace = FALSE)]

    # Estimate cumulative incidence curves.
    fit_perm <- cuminc(
      ftime = data_perm$time,
      fstatus = data_perm$status,
      group = data_perm$arm
    )

    # Permutation quantiles.
    qb0 <- FindQuantile(
      times = fit_perm$`0 1`$time,
      probs = fit_perm$`0 1`$est,
      q = q
    )
    qb1 <- FindQuantile(
      times = fit_perm$`1 1`$time,
      probs = fit_perm$`1 1`$est,
      q = q
    )

    # Bootstrap difference.
    dperm <- qb1 - qb0

    # Bootstrap difference as or more extreme
    out <- 1 * (abs(dperm) >= abs(dobs))
    return(out)
  }

  perm <- lapply(seq(1:B), aux)
  perm <- do.call(c, perm)

  # Permutation p-value.
  perm <- c(1, perm)
  p <- mean(perm)

  # Output.
  out <- c("Quantile" = q, "Arm1" = q1, "Arm0" = q0, "Delta" = dobs, "P" = p)
  return(out)
}
