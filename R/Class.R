#' Compare CICs Object
#'
#' Defines the object class returned by \code{\link{CompareCICs}}. 
#'
#' @slot CIs Confidence intervals. 
#' @slot Curves Cumulative incidence curves. 
#' @slot Pvals Bootstrap and permutation p-values.
#' @slot Reps Bootstrap and permutation statistics. 
#' @slot Stats Observed test statistics. 
#' @name compCICs-class
#' @rdname compCICs-class
#' @exportClass compCICs

setClass(
  Class = "compCICs",
  representation = representation(
   CIs = "data.frame",
   Curves = "data.frame",
   Pvals = "data.frame",
   Reps = "matrix",
   Stats = "data.frame"
  )
)

# -----------------------------------------------------------------------------
# Print Method
# -----------------------------------------------------------------------------

#' Print Method for Compre CICs Object.
#'
#' Print method for objects of class \code{compCICs}.
#'
#' @param x An object of class \code{compCICs}.
#' @param ... Unused.
#' @export

print.compCICs <- function (x, ...) {
  
  # Stats.
  cat('Stats:\n')
  show(x@Stats)
  cat('\n\n')
  
  # CIs.
  cat('CIs:\n')
  show(x@CIs)
  cat('\n\n')
  
  # P-values.
  cat('P-values:\n')
  show(x@Pvals)
  cat('\n\n')

}

# -----------------------------------------------------------------------------
# Show Method
# -----------------------------------------------------------------------------

#' Show Method for Compare CICs Object
#'
#' @param object An object of class \code{compCICs}.
#' @rdname fit-method
#' @importFrom methods show

setMethod(
  f = "show",
  signature = c(object = "compCICs"),
  definition = function (object) {print.compCICs(x = object)}
)

