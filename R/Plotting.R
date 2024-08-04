# Purpose: Plot the cumulative incidence curves comparing two treatment arms.
# Updated: 2021-05-02

# -----------------------------------------------------------------------------

#' Cumulative Incidence Curve.
#' 
#' Return a function that calculates the event probability for a single
#' treatment arm. Assumes status == 1 is the event of interest.
#' 
#' @param data Data.frame.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return stepfun.
#' @export
CICurve <- function(
  data,
  status_name = "status",
  time_name = "time"
) {
  
  # Prepare data.
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Fit cumulative incidence curve.
  fit <- CalcCIC(status = data$status, time = data$time)
  g <- stats::stepfun(x = fit$time, y = c(0, fit$cic_event))
  return(g)
}


#' Standard Error of Cumulative Incidence Curve.
#' 
#' Return a function that calculates the standard error of the event probability
#' for a single treatment arm. Assumes status == 1 is the event of interest.
#' 
#' @param data Data.frame.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return stepfun.
#' @export
SECurve <- function(
  data,
  status_name = "status",
  time_name = "time"
) {
  
  # Prepare data.
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Fit cumulative incidence curve.
  fit <- CalcCIC(status = data$status, time = data$time)
  g <- stats::stepfun(x = fit$time, y = c(0, fit$se_cic_event))
  return(g)
}


#' Number at Risk Curve
#' 
#' Return a function that calculates the number at risk for a single
#' treatment arm.
#' 
#' @param data Data.frame.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return stepfun.
#' @export
NARCurve <- function(
  data,
  status_name = "status",
  time_name = "time"
) {
  
  # Prepare data.
  df <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Fit cumulative incidence curve.
  fit <- CalcCIC(status = data$status, time = data$time)
  g <- stats::stepfun(x = fit$time, y = c(nrow(df), fit$nar))
  return(g)
}


# -----------------------------------------------------------------------------

#' Cumulative Incidence Curve Plotting Frame
#' 
#' Construct a data.frame containing the cumulative incidence of the 
#' status == 1 event across time for a single treatment arm.
#' 
#' @param data Data.frame.
#' @param eval_points Number of points at which to evaluate the curve.
#' @param status_name Name of status column.
#' @param tau Trunction time.
#' @param time_name Name of time column.
#' @return Data.frame.
CICPlotFrame1 <- function(
  data,
  eval_points = 1000,
  status_name = "status", 
  tau = NULL,
  time_name = "time"
) {
  
  # Prepare data.
  data <- data %>%
    dplyr::rename(
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Truncation time.
  if (is.null(tau)) {
    tau <- max(data$time)
  }
  
  # Prepare data.
  cic <- data %>% CICurve()
  times <- seq(from = 0, to = tau, length.out = eval_points)
  out <- data.frame(
    time = times,
    prob = cic(times)
  )
  return(out)
}


#' Cumulative Incidence Curve Plotting Frame
#' 
#' Construct a data.frame containing the cumulative incidence of the 
#' status == 1 event across time for a single treatment arm.
#' 
#' @param data Data.frame.
#' @param arm_name Name of arm column.
#' @param eval_points Number of points at which to evaluate the curve.
#' @param status_name Name of status column.
#' @param tau Trunction time.
#' @param time_name Name of time column.
#' @return Data.frame.
CICPlotFrame2 <- function(
  data,
  arm_name = "arm",
  eval_points = 1000,
  status_name = "status", 
  tau = NULL,
  time_name = "time"
) {
  
  # Prepare data.
  data <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Truncation time.
  if (is.null(tau)) {
    tau <- max(data$time)
  }
  
  # Prepare data.
  arm <- NULL
  df0 <- data %>% 
    dplyr::filter(arm == 0) %>%
    CICPlotFrame1(tau = tau) %>%
    dplyr::mutate(arm = 0)
  
  df1 <- data %>% 
    dplyr::filter(arm == 1) %>%
    CICPlotFrame1(tau = tau) %>%
    dplyr::mutate(arm = 1)
  
  out <- rbind(df0, df1)
  out$arm <- factor(out$arm, levels = c(0, 1), ordered = TRUE)
  return(out)
}


#' Number at Risk Plotting Frame
#' 
#' Numbers at risk for competing risks data.
#' 
#' @param data Data.frame.
#' @param x_breaks Time points at which to determine the NARs.
#' @param arm_name Name of arm column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @return Data.frame containing `time`, `nar_ctrl`, `nar_trt`.
NARPlotFrame <- function(
  data, 
  x_breaks, 
  arm_name = "arm",
  status_name = "status",
  time_name = "time"
) {
  
  # Prepare data.
  df <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # NAR functions.
  arm <- NULL
  g0 <- df %>% dplyr::filter(arm == 0) %>% NARCurve()
  g1 <- df %>% dplyr::filter(arm == 1) %>% NARCurve()
  
  # Output.
  out <- data.frame(
    time = x_breaks,
    nar_ctrl = g0(x_breaks),
    nar_trt = g1(x_breaks)
  )
  return(out)
}

# -----------------------------------------------------------------------------
# Plot cumulative incidence curve.
# -----------------------------------------------------------------------------

#' Plot Cumulative Incidence Curves
#' 
#' Plot the cumulative incidence curves comparing two treatment arms.
#'
#' @param data Data.frame.
#' @param arm_name Name of arm column.
#' @param color_labs Color labels.
#' @param legend_pos Legend position.
#' @param tau Truncation time.
#' @param ctrl_color Color for control arm.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param title Plot title.
#' @param trt_color Color for treatment arm.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis labels.
#' @param x_name X-axis name.
#' @param x_max X-axis upper limit, may differ from tau.
#' @param y_name Y-axis name.
#' @param y_lim Y-axis limits.
#' @return ggplot.
#' @export
PlotCICs <- function(
  data,
  arm_name = "arm",
  color_labs = c("Ctrl", "Trt"),
  ctrl_color = "#C65842",
  legend_pos = "top",
  status_name = "status",
  tau = NULL,
  time_name = "time",
  title = NULL,
  trt_color = "#6385B8",
  x_breaks = NULL,
  x_labs = NULL,
  x_name = "Time",
  x_max = NULL,
  y_name = "Cumulative Incidence",
  y_lim = c(0, 1)
) {
  
  # Prepare data.
  data <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    )
  
  # Set defaults.
  if (is.null(x_max)) {
    x_max <- max(data$time)
  }
  if (is.null(tau)) {
    tau <- x_max
  }
  
  # Plotting frame.
  df <- data %>% CICPlotFrame2(tau = tau)
  
  # Plotting.
  arm <- NULL
  prob <- NULL
  time <- NULL
  q <- ggplot2::ggplot() +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend_pos
    ) + 
    ggplot2::geom_step(
      data = df, 
      ggplot2::aes(x = time, y = prob, color = arm), 
      size = 1) + 
    ggplot2::scale_color_manual(
      name = NULL,
      values = c(ctrl_color, trt_color),
      labels = color_labs
    ) + 
    ggplot2::scale_x_continuous(
      name = x_name,
      breaks = x_breaks,
      labels = x_labs,
      limits = c(0, x_max)
    ) +
    ggplot2::scale_y_continuous(
      name = y_name,
      limits = y_lim
    ) + 
    ggplot2::ggtitle(
      label = title
    )
  
  # Output.
  return(q)
}


# -----------------------------------------------------------------------------

#' Plot AUCIC.
#' 
#' Plot area under the cumulative incidence curve for a single treatment arm.
#'
#' @param data Data including time, status, arm.
#' @param arm_label Label for the arm.
#' @param arm_name Name of arm column.
#' @param color Color.
#' @param legend_pos Legend position.
#' @param status_name Name of status column.
#' @param tau Truncation time for shading.
#' @param time_name Name of time column.
#' @param title Plot title.
#' @param which_arm Arm to plot.
#' @param x_breaks X-axis breaks.
#' @param x_labs X-axis labels.
#' @param x_max X-axis upper limit; may differ from tau.
#' @param x_name X-axis name.
#' @param y_name Y-axis name.
#' @param y_lim Y-axis limits.
#' @return ggplot.
#' @export
PlotAUCIC <- function(
  data,
  arm_label = "Ctrl",
  arm_name = "arm",
  color = "#C65842",
  legend_pos = c(0.2, 0.8),
  status_name = "status",
  time_name = "time",
  title = NULL,
  tau = NULL,
  which_arm = 0,
  x_breaks = NULL,
  x_labs = NULL,
  x_name = "Time",
  x_max = NULL,
  y_name = "Cumulative Incidence",
  y_lim = c(0, 1)
) {
  
  # Prepare data.
  arm <- NULL
  df <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    ) %>%
    dplyr::filter(arm == which_arm) %>%
    CICPlotFrame1() %>%
    dplyr::mutate(arm = factor(which_arm))
  
  # Set defaults.
  if (is.null(x_max)) {
    x_max <- max(data$time)
  }
  if (is.null(tau)) {
    tau <- x_max
  }

  # Plotting.
  prob <- NULL
  time <- NULL
  q <- ggplot2::ggplot(data = df) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = legend_pos
    ) + 
    ggplot2::geom_ribbon(
      ggplot2::aes(x = time, ymin = 0, ymax = prob, fill = arm),
      alpha = 0.5
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = color,
      labels = arm_label
    ) +
    ggplot2::geom_step(
      ggplot2::aes(x = time, y = prob), 
      color = color,
      size = 1
    ) + 
    ggplot2::scale_x_continuous(
      name = x_name,
      breaks = x_breaks,
      labels = x_labs,
      limits = c(0, x_max)
    ) +
    ggplot2::scale_y_continuous(
      name = y_name,
      limits = y_lim
    ) + 
    ggplot2::ggtitle(
      label = title
    )
  
  # Output.
  return(q)
}


# -----------------------------------------------------------------------------

#' Plot Numbers at Risk
#' 
#' @param data Data.frame.
#' @param x_breaks X-axis breaks.
#' @param arm_name Name of arm column.
#' @param status_name Name of status column.
#' @param time_name Name of time column.
#' @param x_labs X-axis labels.
#' @param x_max X-axis upper limit.
#' @param x_name X-axis name.
#' @param y_labs Y-axis labels.
#' @return ggplot.
#' @export
PlotNARs <- function(
  data,
  x_breaks,
  arm_name = "arm",
  status_name = "status",
  time_name = "time",
  x_labs = NULL,
  x_max = NULL,
  x_name = NULL,
  y_labs = c("Ctrl", "Trt")
) {
  
  # Defaults.
  if (is.null(x_labs)) {
    x_labs = x_breaks
  }
  if (is.null(x_max)) {
    x_max = max(x_breaks)
  }
  
  # Data prep.
  nar_ctrl <- NULL
  nar_trt <- NULL
  df <- data %>%
    dplyr::rename(
      arm = {{arm_name}},
      status = {{status_name}},
      time = {{time_name}}
    ) %>%
    NARPlotFrame(x_breaks) %>%
    tidyr::pivot_longer(
      cols = c(nar_ctrl, nar_trt),
      names_to = "arm",
      values_to = "nar"
    ) %>%
    dplyr::mutate(
      arm = factor(arm, c("nar_ctrl", "nar_trt"), y_labs)
    )
  
  # Plotting.
  arm <- NULL
  nar <- NULL
  time <- NULL
  q <- ggplot2::ggplot(data = df) +
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = time, y = arm, label = nar)
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      name = x_name,
      labels = x_labs,
      limits = c(0, x_max)
    ) + 
    ggplot2::scale_y_discrete(
      name = NULL,
      labels = y_labs
    )
  return(q)
}
