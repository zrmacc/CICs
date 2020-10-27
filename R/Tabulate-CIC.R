# Purpose: Tabulate fitted CIC.

#' Tabulate CIC
#' 
#' Converts the list returned by \code{\link[cmprsk]{cuminc}} into a table.
#' 
#' @param time Observation time.
#' @param status Event status. The event coded as 1 is assumed to be the event
#'   of interest.
#' @importFrom cmprsk cuminc
#' @return Data.frame containing `Time`, the estimated `CIC`, the variance `Var_CIC` of the
#'   CIC, and event `Status`.

TabulateCIC <- function(time, status) {
  
  # Cumulative incidence curve.
  fit <- cuminc(ftime = time, fstatus = status)
  
  # Component names. 
  labs <- names(fit)
  
  # Convert to table. 
  aux <- function(i) {
    status <- gsub(pattern = '.+\\ ([0-9]+)', replacement = '\\1', x = labs[i])
    out <- data.frame(do.call(cbind, fit[[i]]))
    colnames(out) <- c('Time', 'CIC', 'Var_CIC')
    out$Status <- as.numeric(status)
    return(out)
  }

  # Output.
  out <- lapply(seq_along(fit), aux)
  out <- do.call(rbind, out)
  return(out)
}
