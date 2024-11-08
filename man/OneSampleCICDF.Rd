% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting.R
\name{OneSampleCICDF}
\alias{OneSampleCICDF}
\title{One-Sample Cumulative Incidence Curve Plotting Frame}
\usage{
OneSampleCICDF(
  data,
  eval_points = 1000,
  status_name = "status",
  tau = NULL,
  time_name = "time"
)
}
\arguments{
\item{data}{Data.frame.}

\item{eval_points}{Number of points at which to evaluate the curve.}

\item{status_name}{Name of status column.}

\item{tau}{Trunction time.}

\item{time_name}{Name of time column.}
}
\value{
Data.frame.
}
\description{
Construct a data frame containing the cumulative incidence of the 
status == 1 event across time for a single treatment arm.
}
