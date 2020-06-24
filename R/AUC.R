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


#' Generate Bootstrap Confidence Interval for the Difference in AUCs.
#'
#' Bootstrap confidence interval for the difference in areas under the
#' cumulative incidence curves.
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param tau Truncation time.
#' @param B Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom cmprsk cuminc
#' @import stats
#' @export
#' @return Numeric vector containing the truncation time `Time`, the
#'   estimated AUCs in arms 1 `Arm1` and 0 `Arm0`, the difference `Delta`,
#'   and the lower `L` and upper `U` confidence bounds.

AUCDiffCI <- function(time, status, arm, tau, B = 2000, alpha = 0.05) {

  # Data.
  data <- data.frame(time, status, arm)

  # Split data.
  data1 <- data[data$arm == 1, ]
  data0 <- data[data$arm == 0, ]

  # Fit competing risks model.
  fit <- cuminc(ftime = data$time, fstatus = data$status, group = data$arm)

  # Cumulative incidence curves.
  a0 <- FindAUC(times = fit$`0 1`$time, probs = fit$`0 1`$est, tau = tau)
  a1 <- FindAUC(times = fit$`1 1`$time, probs = fit$`1 1`$est, tau = tau)

  # Observed difference.
  dobs <- a1 - a0

  # Bootstrap function.
  aux <- function(b) {

    # Bootstrap data sets.
    boot1 <- BootData(data1)
    boot0 <- BootData(data0)
    boot <- rbind(boot1, boot0)

    # Estimate cumulative incidence curves.
    fit_boot <- cuminc(ftime = boot$time, fstatus = boot$status, group = boot$arm)

    # Bootstrap AUC.
    ab0 <- FindAUC(
      times = fit_boot$`0 1`$time,
      probs = fit_boot$`0 1`$est,
      tau = tau
    )
    ab1 <- FindAUC(
      times = fit_boot$`1 1`$time,
      probs = fit_boot$`1 1`$est,
      tau = tau
    )

    # Bootstrap difference.
    dboot <- ab1 - ab0
    return(dboot)
  }

  boot <- lapply(seq(1:B), aux)
  boot <- do.call(c, boot)

  # Confidence interval
  alpha2 <- alpha / 2
  L <- as.numeric(quantile(boot, alpha2, na.rm = TRUE))
  U <- as.numeric(quantile(boot, 1 - alpha2, na.rm = TRUE))

  # Output.
  out <- c("Time" = tau, "Arm1" = a1, "Arm0" = a0, "Delta" = dobs, "L" = L, "U" = U)
  return(out)
}


#' Generate Permutation P-value for the Difference in AUCs
#'
#' Permutation p-value for the difference in areas under the
#' cumulative incidence curves.
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param tau Truncation time.
#' @param B Permutations.
#' @importFrom cmprsk cuminc
#' @import stats
#' @export
#' @return Numeric vector containing the truncation time `Time`, the
#'   estimated AUCs in arms 1 `Arm1` and 0 `Arm0`, the difference `Delta`,
#'   and the permutation p-value `P`.

AUCDiffP <- function(time, status, arm, tau, B = 2000) {

  # Data.
  data <- data.frame(time, status, arm)

  # Split data.
  data1 <- data[data$arm == 1, ]
  data0 <- data[data$arm == 0, ]

  # Fit competing risks model.
  fit <- cuminc(ftime = data$time, fstatus = data$status, group = data$arm)

  # Cumulative incidence curves.
  a0 <- FindAUC(times = fit$`0 1`$time, probs = fit$`0 1`$est, tau = tau)
  a1 <- FindAUC(times = fit$`1 1`$time, probs = fit$`1 1`$est, tau = tau)

  # Observed difference.
  dobs <- a1 - a0

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

    # Permutation areas.
    ab0 <- FindAUC(
      times = fit_perm$`0 1`$time,
      probs = fit_perm$`0 1`$est,
      tau = tau
    )
    ab1 <- FindAUC(
      times = fit_perm$`1 1`$time,
      probs = fit_perm$`1 1`$est,
      tau = tau
    )

    # Bootstrap difference.
    dperm <- ab1 - ab0

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
  out <- c("Time" = tau, "Arm1" = a1, "Arm0" = a0, "Delta" = dobs, "P" = p)
  return(out)
}
