# Purpose: Functions to perform bootstrap inference on the difference in event rates
# between two cumulative incidence curves.

#' Find Event Rate at Time tau.
#'
#' @param times CIC times.
#' @param probs CIC probabilities.
#' @param tau Event time.
#' @importFrom stats stepfun

EventRate <- function(times, probs, tau) {
  g <- stepfun(x = times, y = c(0, probs), right = TRUE)
  return(g(tau))
}


#' Generate Bootstrap Confidence Interval for the Event Rate Difference.
#'
#' Bootstrap confidence interval for the event rate difference at time tau.
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param tau Event time.
#' @param B Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom cmprsk cuminc
#' @import stats
#' @export
#' @return Numeric vector containing the event time `Time`, the estimated rates
#'   in arms 1 `Arm1` and 0 `Arm0`, the difference `Delta`, and the lower `L`
#'   and upper `U` confidence bounds.

EventRateDiffCI <- function(time, status, arm, tau, B = 2000, alpha = 0.05) {

  # Data.
  data <- data.frame(time, status, arm)

  # Split data.
  data0 <- data[data$arm == 0, ]
  data1 <- data[data$arm == 1, ]

  # Fit competing risks model.
  fit <- cuminc(ftime = data$time, fstatus = data$status, group = data$arm)

  # Cumulative incidence curves.
  q0 <- EventRate(fit$`0 1`$time, fit$`0 1`$est, tau = tau)
  q1 <- EventRate(fit$`1 1`$time, fit$`1 1`$est, tau = tau)

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
    qb0 <- EventRate(
      times = fit_boot$`0 1`$time,
      probs = fit_boot$`0 1`$est,
      tau = tau
    )
    qb1 <- EventRate(
      times = fit_boot$`1 1`$time,
      probs = fit_boot$`1 1`$est,
      tau = tau
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
  out <- c("Time" = tau, "Arm1" = q1, "Arm0" = q0, "Delta" = dobs, "L" = L, "U" = U)
  return(out)
}


#' Generate Permutation P-value for the Event Rate Difference.
#'
#' Permutation p-value the difference in event rates at time time tau.
#'
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @param arm Arm, assumed to have two levels, coded 0/1.
#' @param tau Event time.
#' @param B Permutations.
#' @importFrom cmprsk cuminc
#' @import stats
#' @export
#' @return Numeric vector containing the event time `Time`, the estimated rates
#'   in arms 1 `Arm1` and 0 `Arm0`, the difference `Delta`, and the permutation
#'   p-value `P`.

EventRateDiffP <- function(time, status, arm, tau, B = 2000) {

  # Data.
  data <- data.frame(time, status, arm)

  # Split data.
  data1 <- data[data$arm == 1, ]
  data0 <- data[data$arm == 0, ]

  # Fit competing risks model.
  fit <- cuminc(ftime = data$time, fstatus = data$status, group = data$arm)

  # Cumulative incidence curves.
  q0 <- EventRate(fit$`0 1`$time, fit$`0 1`$est, tau = tau)
  q1 <- EventRate(fit$`1 1`$time, fit$`1 1`$est, tau = tau)

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
    qb0 <- EventRate(
      times = fit_perm$`0 1`$time,
      probs = fit_perm$`0 1`$est,
      tau = tau
    )
    qb1 <- EventRate(
      times = fit_perm$`1 1`$time,
      probs = fit_perm$`1 1`$est,
      tau = tau
    )

    # Bootstrap difference.
    dperm <- qb1 - qb0

    # Bootstrap difference as or more extreme
    out <- 1 * (abs(dperm) >= abs(dobs))
    return(out)
  }

  perm <- lapply(seq(1:B), aux)
  perm <- do.call(c, perm)

  # Permutation p-vaue.
  perm <- c(1, perm)
  p <- mean(perm)

  # Output.
  out <- c("Time" = tau, "Arm1" = q1, "Arm0" = q0, "Delta" = dobs, "P" = p)
  return(out)
}
