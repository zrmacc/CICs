# Purpose: Check CIC standard error calculation.
# Updated: 2024-08-03

library(CICs)
library(dplyr)
library(ggplot2)

#' Data Generating Process
#' @param n Sample size.
DGP <- function(n = n) {
  data <- GenData(
    n = n,
    event_rate = 0.50,
    death_rate = 0.25,
    censor_rate = 0.25,
    tau = 6
  )
}


#' Get Point Estimate
#' @param cic Output of `CalcCIC`.
#' @param tau Evaluation time.
GetPointEst <- function(cic, tau) {
  time <- max(cic$time[cic$time < tau])
  return(cic$cic_event[cic$time == time])
}


#' Get Asymptotic Variance
#' @param cic Output of `CalcCIC`.
#' @param tau Evaluation time.
GetAsympVar <- function(cic, tau) {
  time <- max(cic$time[cic$time < tau])
  return(cic$var_cic_event[cic$time == time])
}


#' Get Influence Function Variance
#' @param data Data.
#' @param tau Evaluation Time.
GetIFVar <- function(data, tau) {
  n <- nrow(data)
  time <- max(data$time[data$time < tau])
  inf <- InfluenceCIC(status = data$status, time = data$time, trunc_time = time)
  return(mean(inf^2) / n)
}



#' Simulation Instance
#' @param eval_times Evaluation times.
#' @param n Sample size.
SimInst <- function(eval_times, n) {
  
  data <- DGP(n = n)
  cic <- CalcCIC(status = data$status, time = data$time)
  
  results <- lapply(eval_times, function(tau) {
    out <- data.frame(
      time = tau,
      point = GetPointEst(cic, tau = tau),
      asymp = GetAsympVar(cic, tau = tau),
      inf = GetIFVar(data = data, tau = tau)
    )
    return(out)
  })
  results <- do.call(rbind, results)
  return(results)
}


#' Simulation
#' @param eval_times Evaluation times.
#' @param n Sample size.
#' @param reps Simulation replicates.
Sim <- function(eval_times, n, reps) {
  sim <- lapply(seq_len(reps), function(i) {
    result <- tryCatch(
      SimInst(eval_times = eval_times, n = n),
      error = function(cond) {return(NA)}
    )
    return(result)
  })
  sim <- do.call(rbind, sim)
  
  out <- sim %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      empirical_se = sd(point, na.rm = TRUE),
      asymptotic_se = sqrt(mean(asymp, na.rm = TRUE)),
      influence_se = sqrt(mean(inf, na.rm = TRUE))
    )
  return(out)
}


# -----------------------------------------------------------------------------

# Evaluation times.
eval_times <- seq(from = 0.5, to = 2.0, by = 0.5)

# Sample sizes.
ns <- c(50, 100, 250, 1000)

# Run simulation.
set.seed(1010)
out <- lapply(ns, function(n) {
  result <- Sim(eval_times = eval_times, n = n, reps = 2e3)
  result$n <- n
  return(result)
})
out <- do.call(rbind, out)


# -----------------------------------------------------------------------------

df <- out %>%
  tidyr::pivot_longer(
    empirical_se:influence_se,
    names_to = "estimator",
    values_to = "value"
  )
df$n <- glue::glue("N={df$n}")
df$n <- factor(
  x = df$n,
  levels = paste0("N=", ns),
  ordered = TRUE
)

# Plotting.
q <- ggplot(data = df) + 
  theme_bw() + 
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank()
  ) + 
  theme(
    axis.line = element_line(color = 'black')
  ) +
  geom_col(
    aes(x = time, y = value, fill = estimator),
    position = position_dodge()
  ) + 
  facet_wrap(~n) +
  ggsci::scale_fill_nejm(
    name = "Estimator",
    labels = c("Asymptotic", "Empirical", "Influence Function")
  ) + 
  scale_x_continuous(
    name = "Evaluation Time"
  ) + 
  scale_y_continuous(
    name = "Standard Error"
  )
  
show(q)  
if (FALSE) {
  ggsave(
    plot = q,
    file = "se_validation.pdf",
    width = 9.0,
    height = 4.5
  )
}



