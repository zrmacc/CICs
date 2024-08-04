#' Generate Bootstrap Data.frame
#'
#' @param data Data.frame.
#' @return Bootstrapped data.frame.
BootData <- function(data) {
  n <- nrow(data)
  key <- sample(x = n, size = n, replace = TRUE)
  out <- data[key, ]
  return(out)
}


#' Stratified Bootstrap
#' 
#' @param data Data.frame containing `strata`
#' @return Data.frame bootstrapped within strata.
StratBoot <- function(data) {
  data_strata <- split(data, data$strata, drop = TRUE)
  boot_strata <- lapply(data_strata, BootData)
  out <- do.call(rbind, boot_strata)
  return(out)
}


# -----------------------------------------------------------------------------

#' Generate Permuted Data.frame
#'
#' @param data Data.frame containing `arm`.
#' @return Permuted data.frame.
PermData <- function(data) {
  n <- nrow(data)
  data$arm <- data$arm[sample(x = n, size = n, replace = TRUE)]
  return(data)
}


#' Stratified Permutation
#' 
#' @param data Data.frame containing `strata` and `arm`.
#' @return Data.frame bootstrapped within strata.
StratPerm <- function(data) {
  data_strata <- split(data, data$strata, drop = TRUE)
  perm_strata <- lapply(data_strata, PermData)
  out <- do.call(rbind, perm_strata)
  return(out)
}