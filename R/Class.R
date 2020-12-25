#' Compare CICs Object
#'
#' Defines the object class returned by \code{\link{CompareCICs}}. 
#'
#' @slot cis Confidence intervals. 
#' @slot cics Cumulative incidence curves. 
#' @slot pvals Bootstrap and permutation p-values.
#' @slot reps Bootstrap and permutation statistics. 
#' @slot marg_stats Observed test statistics. 
#' @slot stratum_stats Per-stratum weights and statistics. 
#' 
#' @name compCICs-class
#' @rdname compCICs-class
#' @exportClass compCICs

setClass(
  Class = "compCICs",
  representation = representation(
   cis = "data.frame",
   cics = "data.frame",
   pvals = "data.frame",
   reps = "list",
   marg_stats = "data.frame",
   stratum_stats = "data.frame"
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
  cat('Marginal Stats:\n')
  show(x@marg_stats)
  cat('\n\n')
  
  # CIs.
  cat('CIs:\n')
  show(x@cis)
  cat('\n\n')
  
  # P-values.
  cat('P-values:\n')
  show(x@pvals)
  cat('\n\n')
  
  # Weights.
  cat('Stratum Stats:\n')
  show(x@stratum_stats)
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

