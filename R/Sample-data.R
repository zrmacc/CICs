#' Sample Competing Risks Data
#'
#' Synthetic time-to event data for two arms of 600 patients, formatted
#' as expected by this package. 
#'
#' @docType data
#' @usage data(cic_data)
#' @format A data.frame containing three fields:
#' \describe{
#'   \item{arm}{Treatment arm, 0 for reference, 1 for treatment.}
#'   \item{time}{Time-to event, integer between 1 and 30.}
#'   \item{status}{Event type, 0 for censoring, 1 for recovery, 2 for death.}
#' }
"cic_data"