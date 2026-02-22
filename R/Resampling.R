# Internal: apply fun within strata and rbind.
.stratified_apply <- function(data, fun) {
  strata_list <- split(data, data$strata, drop = TRUE)
  do.call(rbind, lapply(strata_list, fun))
}

#' Generate Bootstrap Data.frame
#'
#' @param data Data.frame.
#' @return Bootstrapped data.frame.
BootData <- function(data) {
  n <- nrow(data)
  data[sample.int(n, n, replace = TRUE), ]
}

#' Stratified Bootstrap
#'
#' @param data Data.frame containing `strata`
#' @return Data.frame bootstrapped within strata.
StratBoot <- function(data) {
  .stratified_apply(data, BootData)
}

#' Generate Permuted Data.frame
#'
#' @param data Data.frame containing `arm`.
#' @return Permuted data.frame.
PermData <- function(data) {
  n <- nrow(data)
  data$arm <- data$arm[sample.int(n, n, replace = FALSE)]
  data
}

#' Stratified Permutation
#'
#' @param data Data.frame containing `strata` and `arm`.
#' @return Data.frame with arm permuted within strata.
StratPerm <- function(data) {
  .stratified_apply(data, PermData)
}