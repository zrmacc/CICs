% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIC.R
\name{TabulateEvents}
\alias{TabulateEvents}
\title{Tabulate Events}
\usage{
TabulateEvents(time, status)
}
\arguments{
\item{time}{Event time.}

\item{status}{Status, coded as 1 for the recurrent event, 0 otherwise.}
}
\value{
Data.frame with the censorings, deaths, and events occurring
  at each distinct time point.
}
\description{
Tabulate the number at risk and the number of events at each unique
observation time.
}
