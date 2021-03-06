% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting.R
\name{PlotAUCIC}
\alias{PlotAUCIC}
\title{Plot AUCIC.}
\usage{
PlotAUCIC(
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
)
}
\arguments{
\item{data}{Data including time, status, arm.}

\item{arm_label}{Label for the arm.}

\item{arm_name}{Name of arm column.}

\item{color}{Color.}

\item{legend_pos}{Legend position.}

\item{status_name}{Name of status column.}

\item{time_name}{Name of time column.}

\item{title}{Plot title.}

\item{tau}{Truncation time for shading.}

\item{which_arm}{Arm to plot.}

\item{x_breaks}{X-axis breaks.}

\item{x_labs}{X-axis labels.}

\item{x_name}{X-axis name.}

\item{x_max}{X-axis upper limit; may differ from tau.}

\item{y_name}{Y-axis name.}

\item{y_lim}{Y-axis limits.}
}
\value{
ggplot.
}
\description{
Plot area under the cumulative incidence curve for a single treatment arm.
}
