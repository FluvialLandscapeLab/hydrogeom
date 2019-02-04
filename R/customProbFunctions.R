#' Example shape functions
#'
#' Transient storage zone statistics are calucated based on an assumed shape of
#' the hydrologic residence time distribution (RTD) of the hyporheic zone, which
#' is the probablity density function (PDF) that a water molecule will exit the
#' hyporheic zone at time \code{tau}.  \code{powerLaw} and \code{exponent} are
#' examples of functions that describe a basic shape for the RTD.
#'
#' Given an appropriate "shape function," \code{HydroGeom} can calculate
#' transient storage zone statistics for any assumed shape of the RTD. Shape
#' functions must have the signature \code{function(tau, tau_n, ...)} and
#' typically return zero when tau == tau_n.
#'
#' You can generate hyporheic TSZ statistics based on a RTD with a shape
#' corresponding to a custom function that you provide.  If you write a shape
#' function like \code{foo = function(tau, tau_n, chi, omega)}, where \code{chi}
#' and \code{omega} are parameters that influence the shape of the curve
#' described by your function (e.g., analagous to the way \code{alpha} serves as
#' an exponent that controls the specific shape of a power law curve in the
#' \code{powerLaw} function), you can then use \code{\link{TSZStats}} with
#' \code{shape = "foo"} to generate TSZ statistics using a custom shape for the
#' RTD of the hyporheic zone.  Implemented this way, results of \code{\link{TSZStats}}
#' are derived using numerical integration (see \code{\link{numericalSolution}}).
#' Solutions for integration of a custom function can also be provide (see \code{\link{fullSolution}} and
#' \code{\link{powerLawPDF}}).
#'
#' Technically speaking, \code{tau_n} serves as the upper limit of \code{tau}
#' used to calculate the normalizing constant of the PDF, which is defined as
#' the inverse of the definite integral of the shape function from some lower
#' limit (\code{tau_0}) to tau_n.  Pragmatically, \code{tau_n} represents the
#' maximum residence time of interest, which should approximate the maximum
#' subsurface residence time of hyporheic water.  Water that stays in the
#' subsurface longer than \code{tau_n} is assumed to enter the true groundwater
#' system.
#'
#' All paramerters for these fuctions correspond to the defintions found in
#' Poole et al. (In Press) "The Hydrolgic Geometry of Hyporheic Zones in
#' Expansive Coarse-Grained Alluvial Aquifers."
#'
#' @return Shape functions (including custom shape functions you write) must
#'   recycle values of \code{tau} and \code{tau_n}. Effectively, since the
#'   \code{tau_n} is typically length == 1, shape functions usually return a
#'   vector of numeric values equal in length to the \code{tau} vector passed to
#'   the function. Both \code{powerLaw} and \code{exponent} conform to these
#'   requirements. The \code{powerLaw} function returns \code{tau^(-alpha) -
#'   tau_n^(-alpha)}.  The \code{exponent} function returns
#'   \code{exp(-sigma*tau) - exp(-sigma * tau_n)}.  Implementation of these
#'   functions can be viewed by typing \code{print(powerLaw)} or
#'   \code{print(exponent)} into the console.
#'
#' @param tau A numeric vector of residence times.
#' @param tau_n A numeric vector (usually but not necessarily lenth == 1)
#'   representing the maximum \code{tau} of the hyporehic zone.
#' @param alpha The exponent of the power law; this value is negated within the
#'   \code{powerLaw} function.
#' @param sigma The hydrologic exchange rate for the \code{exponent} function;
#'   this value is negated within the \code{exponent} function.
#' @seealso \code{\link{checkShapeFunction}}
#' @export
powerLaw = function(tau, tau_n, alpha) {
  if(any(tau > tau_n)) stop("Values of tau must be less than values of tau_n")
  tau^(-alpha) - tau_n^(-alpha)
}

#' @rdname powerLaw
#' @export
exponent = function(tau, tau_n, sigma) {
  if(any(tau > tau_n)) stop("Values of tau must be less than values of tau_n")
  exp(-sigma*tau) - exp(-sigma * tau_n)
}

#' Calculate and itegrate probability functions for shape functions
#'
#' Given the name of a shape function and the specific calculation requested,
#' calculate either the PDF, CCDF, integration of the PDF or CCDF, or the
#' integration of tau*PDF or tau*CCDF.
#'
#' The function \code{numericalSolution} uses numerical integration of basic
#' shape function (e.g., \code{\link{powerLaw}, \link{exponent}}) to return
#' values of the probability density function (PDF) of the shape function,
#' complementary cumulative distribution funtion (CCDF) of the shape function,
#' finite integrals of the PDF or CCDF, or finite integrals of \code{tau*PDF} or
#' \code{tau*CCDF}.  Numerical integration to calculate the normalizing constant
#' for the PDF and CCDF, and for integrating these functions.
#'
#' The function \code{fullSolution} looks for a function with a name specified
#' by the concatination of \code{shape} and \code{suffix} and checks to be sure
#' any such function has the signature \code{function(tau, tau_0, tau_n, ...)}
#' if suffix is "PDF" or "CCDF", or \code{function(tau_a, tau_b, tau_0, tau_n,
#' ...} for other values of suffix.  Solutions for the "powerLaw" and "exponent"
#' shapes are provided as part of this package as follows:  PDF (e.g.
#' \code{\link{powerLawPDF}, \link{exponentPDF}}), CCDF (e.g.
#' \code{\link{powerLawCCDF}, \link{exponentCCDF}}), integral of the PFD (e.g.
#' \code{\link{powerLawIntPDF}, \link{exponentIntPDF}}), integral of the CCDF
#' (e.g. \code{\link{powerLawIntCCDF}, \link{exponentIntCCDF}}), integral of
#' tau*PDF (e.g. \code{\link{powerLawIntTau.PDF}, \link{exponentIntTau.PDF}}) or
#' integral of tau*CCDF (e.g. \code{\link{powerLawInt.TauCCDF},
#' \link{exponentIntTau.CCDF}}).  These solution functions can also be called
#' directly, but \code{fullSolution} is provided as a convenience wrapper that
#' calls the functions using the convention \code{shape} and \code{suffix}.
#'
#' The user can provide custom shape functions or solutions for other shapes,
#' however the names of the custom solution functions must follow the convention
#' of concatinating the shape name with the suffix as the name of the solution
#' function (e.g., the solution function for the integral of tua*PDF of a shape
#' called "foo" must be \code{fooIntTau.PDF}).  For more information, see the
#' documentation of \code{\link{powerLaw}} and \code{\link{powerLawPDF}}
#'
#' The function \code{autoSolution} looks first for a solution function that
#' matches the concatination of \code{shape} and \code{suffix} and, if found,
#' calls \code{fullSolution}. If the solution functions are not found,
#' \code{autoSolution} calls \code{numericSolution}.
#'
#' @return A vector of values representing either the PDF, CCDF, integration of
#'   the PDF, integration of the CCDF, integration of \code{tau} * PDF or the
#'   integration of \code{tau} * CCDF.
#' @param tau A vector of residence times for which PDF or CCDF values are
#'   requested, or the lower values of definate intergrals of the finite
#'   integral of the PDF, CCDF, tau*PDF, or tau*CCDF.
#' @param tau_0,tau_n Minimum and maximum residence times of the definite
#'   integration.
#' @param tau_b The upper values of definite integrals of the PDF, CCDF,
#'   tau*PDF, or tau*CCDF.  Ignored for calculating PDF or CCDF values.
#' @param ... Additional values required by the shape function.
#' @param shape A character vector of \code{length(shape) == 1} naming the shape
#'   function (e.g., \code{\link{powerLaw}}) to be integrated numerically if
#'   solutions are not available used.
#' @param suffix A character vector indicating the value to be returned; either
#'   "PDF", "CCDF", "IntPDF", "IntCCDF", "IntTau.PDF", or "IntTau.CCDF."
#' @param integrateArgs1,integrateArgs2 When numerical integration is used, a
#'   named list of optional values for \code{\link[stats]{integrate}} function
#'   (e.g., \code{subdivisions, rel.tol}, etc.) can be passed via this argument.
#'   Any values in integrateArgs1 list will be applied to the integration used
#'   within the calculations for the PDF or CCDF.  Values for integrateArgs2
#'   will be used when the PDF or CCDF are integrated.  Thus, integrateArgs1 are
#'   used regardless of the value of \code{suffix}.  integrateArgs2 are used
#'   only when \code{suffix} is "IntPDF", "IntCCDF", "IntTau.PDF", or
#'   "inteTau.CCDF".  When solutions functions are used (rather than numerical
#'   integration), these parameters are ignored.
#' @param forceNumeric When set to T, numerical integration will be used by
#'   \code{autoSolution} even if solutions to the PDF, CCDF, IntPDF, IntCCDF,
#'   IntTau.PDF, and IntTau.CCDF are available.
#' @export
autoSolution = function(suffix, shape, tau_0, tau_n, tau, tau_b = NULL, ..., MoreArgs = list(), integrateArgs1 = list(), integrateArgs2 = list(), forceNumeric = F) {
  fun = NULL
  # if user doesn't force numeric, look for solution functions
  if(!forceNumeric) fun = getShapeFunction(shape, suffix, stopOnMissing = F)
  # if user forces numeric or solution function isn't found, look for shape
  # function (getShapeFunction traps case where function is not found)
  if(is.null(fun)) {
    forceNumeric = T
    fun = getShapeFunction(shape)
  }

  # if forceNumeric, use the numerical integration.
  if(forceNumeric){
    result = numericalSolution(suffix = suffix, shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = tau, tau_b = tau_b, ...,
                               MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2,
                               forceNumeric = forceNumeric)
  # otherwise, call fullSolution.
  } else {
    result = fullSolution(suffix = suffix, shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = tau, tau_b = tau_b, ...,
                      MoreArgs = MoreArgs)

  }
  return(result)
}

#' @export
#' @rdname autoSolution
fullSolution = function(suffix, shape, tau_0, tau_n, tau, tau_b = NULL, ..., MoreArgs = list()) {
  # If tau_b is NULL, (typically because NULL is the default), removing tau_b
  # completely causes a very pleasant "tau_b not found" error to be thrown if
  # tau_b is required by the solution function.  This is easier than, after the
  # fact, trying to catch the error associated with passing NULL to a solution
  # function.
  if(is.null(tau_b)) {
    rm(tau_b)
    checkOrder(tau_0, tau, tau_n)
  } else {
    checkOrder(tau_0, tau, tau_b, tau_n)
  }
  fun = getShapeFunction(shape, suffix, stopOnMissing = T)
  # any function that requires "tau_a" will be an integration of the PDF or CCDF.
  if("tau_a" %in% names(formals(fun))) {
    argList = c(list(tau_a = tau, tau_b = tau_b, tau_0 = tau_0, tau_n = tau_n, ...), MoreArgs)
    # otherwise, it's a direct call to a PDF or CCDF function.
  } else {
    argList = c(list(tau = tau, tau_0 = tau_0, tau_n = tau_n, ...), MoreArgs)
  }
  result = do.call(fun, argList)
  return(result)
}


# A wrapper for autoSolution that allows use with the integrate() function
integratableSolution = function(tau, ...) {
  autoSolution(tau = tau, ...)
}

#' @rdname autoSolution
#' @export
numericalSolution = function(suffix, shape, tau_0, tau_n, tau, tau_b = NULL, ..., MoreArgs = list(), integrateArgs1 = list(), integrateArgs2 = list(), forceNumeric = T){
  # tau_b is not required for suffix = c("PDF", "CCDF").  It is required for all
  # other values of suffix.  By removing tau_b entirely if NULL, the code runs
  # when tau_b is not required, but throws a lovely "tau_b not found" error
  # otherwise.  Because tau_b has a default value, this is safe: tau_b will
  # always be local to this function, so this code won't remove a tau_b value
  # higher in the hierarchy of environments.
  if(is.null(tau_b)) {
    rm(tau_b)
    checkOrder(tau_0, tau, tau_n)
  } else {
    checkOrder(tau_0, tau, tau_b, tau_n)
  }
  fun = getShapeFunction(shape, stopOnMissing = F)
  formalNames = names(formals(fun))
  if (!all(c("tau", "tau_n") %in% formalNames)) {
    stop("The parameters of a shape function must include 'tau' and 'tau_n'.")
  }
  # if the user is asking for the PDF or CCDF value...
  if(grepl("(^P)|(^C{2})DF$", suffix)) {
    # calculate the normalizing constant
    normalizingConstants = 1 / integrationValues(fun, lower = tau_0, upper = tau_n, tau_n = tau_n, ..., MoreArgs = MoreArgs, integrateArgs = integrateArgs1)
    if(suffix == "PDF") {
      result = mapply(fun, tau = tau, tau_n = tau_n, ..., MoreArgs = MoreArgs) * normalizingConstants
    } else if (suffix == "CCDF") {
      result = integrationValues(fun, lower = tau, upper = tau_n, tau_n = tau_n, ..., MoreArgs = MoreArgs, integrateArgs = integrateArgs1) * normalizingConstants
    }
  # if user is asking for integration of PDF or CCDF
  } else if (grepl("^Int(Tau[.])?(P|C{2})DF$", suffix)){
#    if(is.null(tau_b)) stop("tau_b parameter must be specified for '", suffix, ".'")
    # if integration of PDF or CCDF:
    if (grepl("^Int(P|C{2})DF$", suffix)) {
      PDForCCDF = substr(suffix, 4, nchar(suffix))
      f = integratableSolution
    # otherwise, integration of tau*PDF or tau*CCDF
    } else if (grepl("^IntTau[.](P|C{2})DF$", suffix)) {
      PDForCCDF = substr(suffix, 8, nchar(suffix))
      f = function(tau, ...) {
        tau * integratableSolution(tau = tau, ...)
      }
    } else {
      stop("Unexpected suffix: '", suffix, "' was not trapped correctly.  Please send this error message to the package author.")
    }
      autoSolutionArgs =
        list(
          suffix = PDForCCDF,
          shape = shape,
          MoreArgs = MoreArgs,
          integrateArgs1 = integrateArgs1,
          forceNumeric = forceNumeric
        )
      result = integrationValues(f, lower = tau, upper = tau_b, tau_0 = tau_0, tau_n = tau_n, ..., MoreArgs = autoSolutionArgs, integrateArgs = integrateArgs2)
  } else {
    stop("Invalid suffix value '", suffix, "'; valid values are 'PDF', 'CCDF', 'IntPDF', 'IntCCDF', 'IntTau.PDF', and 'IntTau.CCDF'")
  }

  if(length(ls(pattern = "^result$")) == 0) stop("Unexpected suffix: '", suffix, "' was not trapped correctly.  Please send this error message to the package author.")

  return(result)
}

# a function that will recycle values of lower, upper, and anything in ... and
# return results from integrating each combination of values.  If valuesOnly is
# F, a list of integration objects is returned.  Otherwise, the function returns
# a vector of the $value attributes extracted from the list of integration
# objects.
integrationValues = function(f, lower, upper, ..., MoreArgs = NULL, integrateArgs = NULL, valuesOnly = T) {
  MoreArgs = append(list(f = f), MoreArgs)
  MoreArgs = append(MoreArgs, integrateArgs)
  multiIntegration = mapply(integrate, lower = lower, upper = upper, ..., MoreArgs = MoreArgs, SIMPLIFY = F)
  if(valuesOnly) {
    errors = sapply(multiIntegration, "[[", "abs.error")
    multiIntegration = sapply(multiIntegration, "[[", "value")
    fracErrors = errors/abs(multiIntegration)
    fracErrors[errors == 0 & multiIntegration == 0] = 0
    if(any(fracErrors > 0.01)) warning("Absolute error of numerical integration may be as much as ", round(max(fracErrors),3)*100, "% of estimated value.")
  }
  return(multiIntegration)
}
