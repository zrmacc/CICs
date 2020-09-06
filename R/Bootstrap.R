#' Generate Bootstrap Data.frame
#'
#' @param data Data.frame.
#' @return Bootstrapped data.frame.

BootData <- function (data) {
  n <- nrow(data)
  key <- sample(x = n, size = n, replace = TRUE)
  out <- data[key, ]
  return(out)
}

#' Generate Bootstrap Data.frame
#'
#' @param data Data.frame containing `arm`.
#' @return Permuted data.frame.

PermData <- function (data) {
  n <- nrow(data)
  data$arm <- data$arm[sample(x = n, size = n, replace = TRUE)]
  return(data)
}