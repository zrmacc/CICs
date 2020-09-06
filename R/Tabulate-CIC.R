# Purpose: Tabulate fitted CIC.

#' Tabulate CIC
#' 
#' Converts from the list returned by \code{\link[cmprsk]{cuminc}}.
#' 
#' @param fit Fitted object of class "cuminc".
#' @return Data.frame containing `Time`, the estimated `CIC`, the variance `Var_CIC` of the
#'   CIC, and event `Status`.

TabulateCIC <- function (fit) {
  
  # Component names. 
  labs <- names(fit)
  
  # Convert to table. 
  aux <- function(i) {
    status <- gsub(pattern = '.+\\ ([1-9]+)', replacement = '\\1', x = labs[i])
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