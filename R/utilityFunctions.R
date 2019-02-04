#' Test a custom shape function
#'
#' Integrate a shape function from \code{tau_0} to \code{tau_n}.
#'
#' To generate transient storage zone statistics, shape functions (see, e.g.,
#' \code{\link{powerLaw}}) can be integrated numerically (passed to
#' \code{\link[stats]{integrate}}).  Occasionally, depending on the units of
#' time chosen for \code{tau} or for other reasons,
#' \code{\link[stats]{integrate}} will fail, reporting that the shape function
#' is "probabaly divergent" (see examples, below). \code{checkShapeFunction}
#' provides a means to test different integration options when numerical
#' integration fails.
#'
#' The function works by passing the shape function to
#' \code{\link[stats]{integrate}}.  Using this function, you can try different
#' units of \code{tau} (e.g., hours or days rather than seconds) or potentially
#' manipulate the optional parameters of \code{\link[stats]{integrate}} (e.g.,
#' \code{subdivisions}, \code{rel.tol}, etc.) using the \code{integrateArgs}
#' list. See the examples, below.
#'
#' @param shape A shape function or character string containing the name of the
#'   shape function.
#' @param tau_0 Lower time bound of the definite intergal of interest.
#' @param tau_n Upper time bound of the definite integral of interest.
#' @param ... Parameters required by <shape> over which the function will be
#'   vectorized.
#' @param MoreArgs A names list of parameters that will be passed directly to
#'   <shape> without vectorization.  Use this, e.g., if the <shape> file
#'   requires a list to prevent vectorization (and associated recyling) of the
#'   list.
#' @param integrationArgs A named list of optional parameters for
#'   \code{\link[stats]{integrate}}.
#' @param valuesOnly When T, return value will be a numeric vector of values.
#'   When F, return value will be a list of integration objects (see
#'   \code{\link[stats]{integrate}})
#' @examples checkShapeFunction("powerLaw", 60, 3600, alpha = 1.5)
#' checkShapeFunction("powerLaw", 60, 3600, alpha = c(1.3, 1.5, 1.7))
#'
#' # try integrating from 1 minute to 180 days:
#' # when units are seconds, numerical integration fails:
#' checkShapeFunction("powerLaw", 60, 15552000, alpha = 1.5) #FAIL
#' # numerical integration succeeds when units are hours
#' checkShapeFunction("powerLaw", 1/60, 4320, alpha = 1.5)
#'
#' # passing through optional parameters to integrate() function.
#' # Default values for subdivisions is 100 (see help for integrate() function).
#' # Set subdivisions unreasonably low, and integration fails:
#' checkShapeFunction("powerLaw", 60, 3600, alpha = 1.5, integrationArgs = list(subdivisions = 2))
#' # increase subdivisions somewhat and integration succeeds.
#' checkShapeFunction("powerLaw", 60, 3600, alpha = 1.5, integrationArgs = list(subdivisions = 10))
#' @export
checkShapeFunction = function(shape, tau, tau_n, ..., MoreArgs = NULL, integrateArgs = list(), valuesOnly = F) {
  fun = getShapeFunction(shape)
  if(is.null(fun)) {
    stop("Function '", shape, "' either was not found, or does not contain accept 'tau' and 'tau_n' as parameters.")
  }
  integrationValues(fun, lower = tau, upper = tau_n, tau_n = tau_n, ..., MoreArgs = MoreArgs, integrateArgs = integrateArgs, valuesOnly = valuesOnly)
}

# Looks for a function named <paste0(shape, suffix)>.  If found, checks to be
# sure the formals are as expected.  If not found, returns NULL.
getShapeFunction = function(shape, suffix = "", stopOnMissing = T) {
  fun = try(match.fun(paste0(shape, suffix)), silent = T)
  if(mode(fun) == "function") {
    if(suffix == "") {
      # expected formals for a shape function
      expectedFormals = c("tau", "tau_n")
    } else if(grepl("^Int", suffix)) {
      # expected formals for solved integration of PDF or CCDF of a shape function
      expectedFormals = c("tau_a", "tau_b", "tau_0", "tau_n")
    } else {
      # expected formals for solved PDF or CCDF
      expectedFormals = c("tau", "tau_0", "tau_n")
    }
    if(!all(expectedFormals %in% names(formals(fun)))) {
      stop("Function '", shape, suffix, "' doesn't take values '", paste0(expectedFormals, collapse = "', '"), "' as parameters")
    }
  } else {
    if(stopOnMissing) {
      stop("Function '", shape, suffix, "' was not found.")
    } else {
      fun = NULL
    }
  }
  return(fun)
}

# function that checks to be sure every value in each series formed by recycling
# parameters is sorted ascending.
checkOrder = function(...) {
  if(any(mapply(function(...) is.unsorted(c(...)), ...))) stop("Requested values of 'tau' were out of order (not ascending) or outside the specified minimum ('tau_0') or maximum ('tau_n') integration range.")
  return(T)
}

# print method for TSZ class
#' @export
print.TSZs = function(x) {
  class(x) = class(x)[-1]
  print(x)
  estimateIdx = which(sapply(attributes(x), function(attrib) "estimate" %in% names(attrib)))
  cat(
    paste0("Estimated '", names(attributes(x))[estimateIdx], "': ", format(attributes(x)[[estimateIdx]]$estimate, digits = 8),
           " | storage discrepancy: ", format(attr(x, "storage discrepancy"), digits = 8),
           " | q discrepancy: ", format(attr(x, "q discrepancy"), digits = 8))
  )
}

TSZsOrigin = function(TSZs) {
  if(!inherits(TSZs, "TSZs")) stop("Only a 'TSZs' object can be passed to originTSZs")
  remove = which(names(attributes(TSZs)) %in% c("names", "row.names", "class", "q discrepancy", "storage discrepancy"))
  origin = attributes(TSZs)[-remove]
  estimated = which(sapply(origin, function(attrib) "estimate" %in% names(attrib)))
  origin[estimated] = list(NULL)
  return(origin)
}
