#' Generate Bootstrap Data.frame.
#'
#' @param data Data.frame.
#' @return Bootstrapped data.frame.

BootData <- function(data) {
  n <- nrow(data)
  key <- sample(x = n, size = n, replace = TRUE)
  out <- data[key, ]
  return(out)
}
