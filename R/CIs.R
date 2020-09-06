#' Construct Confidence Interval
#' 
#' @param x Bootstrap replicates.
#' @param alpha Alpha level.
#' @return Numeric vector containing the `Alpha` level, and the 
#'   lower `L` and upper `U` confidence bounds. 

ConstructCI <- function (x, alpha = 0.05) {
  
  lower <- quantile(x = x, probs = alpha / 2)
  upper <- quantile(x = x, probs = 1 - alpha / 2)
  
  # Output.
  out <- c(
    "Alpha" = alpha,
    "L" = as.numeric(lower),
    "U" = as.numeric(upper)
  )
  return(out)
}