

# NOTE: fun(tau, tau_0, tau_n, ...) must be a funtion that returns values of the
# PDF.  It will be integrated numerically.

# solution(tau_a, tau_b, tau_0, tau_n, ...) must be a function that returns the
# integration of the CCDF of the PDF.

# In either case, '...' are any other values that are required by fun or
# solution (e.g., alpha for power law).  The first of these extra values can be
# estimated if not passed to TSZStats, but storage and q are both provided.
# Thus, in the case of fun = aPDF(tau, tau_0, tau_n, x, y, z) or solution =
# anIntegratedCCDF(tau_a, tau_b, tau_0, tau_n, x, y, z) x can be estimated, but
# y and z must be provided (and named, e.g. y=3, z=2, when passed to TSZStats).

#' Generate transient storage statistics
#'
#' The exchange of water among and size of multiple transient storage zones
#' within the hyporehic zone are calculated according to the method describe in
#' Poole et al. (In press).  Calculations are made either by utilizing solutions
#' to various integrations of a shape function (e.g., \code{\link{powerLawPDF}}
#' and associated fuctions) or via numerical integration of a shape fuction
#' (e.g., \code{\link{powerLaw}}).
#'
#' When \code{factor == 0}, transient storage zone breakpoints will be
#' determined such that each transient storage zone has and equal amount of
#' water storage.  When factor is 1 or greater, \code{factor} is used as a
#' mutiplier to determine the residence times associated with each TSZ such that
#' the range of residence times contained in any TSZ is \code{<factor>} times
#' larger than the range of residence times contained in the previous TSZ.  So
#' if \code{factor == 2}, the range of residence times in each TSZ will be 2x
#' the range of the prior TSZ.
#'
#' If passed as NULL, the value for first parameter after \code{tau} and
#' \code{tau_n} in the formals of the function named by \code{shape} can be
#' estimated by \code{TSZStats}, however non-NULL values for both \code{storage}
#' and \code{q} must be provided.  For instance, for the function
#' \code{\link{powerLaw}}, a value for \code{alpha} will be estimated by
#' TSZStats if \code{alpha = NULL} is included in the call to \code{TSZStats}
#' along with non-NULL values for \code{storege} and \code{q}.
#'
#' \code{HyporheicBins} has been depricated and is now a wrapper for
#' \code{TSZStats}.  Use \code{TSZStats} instead.
#'
#' @param TSZs Either the number of transient storage zones requested, or a
#'   vector of residence times that serve as the breakpoints among transient
#'   stroage zones.
#' @param tau_0 The smallest residence time to be considered part of the
#'   hyporheic zone.
#' @param tau_n The largest residence time to be considered part of the
#'   hyporehic zone.
#' @param storage Amount of water stored in the hyporheic zone.
#' @param q Gross hyporheic exchange per time unit.  Exchange can be unit can be
#'   distance, area, or volume per unit time, but must be the same as used units
#'   used to describe storage.  Time unit must be the same as units used to
#'   describe tau_0 and tau_n.
#' @param .... Additional arguments required by shape function (e.g.,
#'   \code{alpha} for \code{\link{powerLaw}}).  Must be named.  First of these
#'   values can be passed as NULL and will be estimated if both \code{storage}
#'   and \code{q} are provided.  See details.
#' @param factor Controls the distribution of residence time for transient
#'   storage zones when \code{TSZs} is set equal to the number of transient
#'   storage zones. Value must be either 0 or >1. See details.  \code{factor} is
#'   ignored when TSZs is a vector of the residence time breakpoints between
#'   transient storage zones.
#' @param shape A character string naming the shape file (e.g.,
#'   \code{\link{powerLaw}}) to be used in the calculation of transient storage
#'   statistics.
#' @param MoreArgs A named list of additional parameter values to be passed to
#'   the shape function -- only valid when numerical integration is used;
#'   ignored when solutions to the shape function are used.  When passed this
#'   way, the values will not be vectorized with tau, tau_0, tau_n, tau_b, or
#'   arguments in '...'.
#' @param integrateArgs1,integrateArgs2 When numerical integration is used,
#'   these parameters are a named list of optional values for the
#'   \code{\link[stats]{integrate}} function.  For PDF and CCDF calculations,
#'   integrateArgs1 is used and integrateArgs2 is ignored.  For the integration
#'   of the PDF, CCDF, tau*PDF or tau*CCDF, integrateArgs1 will apply to the
#'   internal call to the PDF or CCDF functions and integrateArgs2 will apply to
#'   the subsequent integration of those functions.  In practise, use of these
#'   parameters is seldom required.  These parameters are useful only if and
#'   unexpected error message is returned by the \code{\link{integrate}}
#'   function.
#' @return A TSZ object, which is a data.frame with attributes that record the
#'   arguments (\code{TSZs, tau_0, tau_n, storage, q, factor, and forceNumeric},
#'   along with any shape specific parameters (e.g., \code{alpha} for
#'   \code{\link{powerLaw}})) used to calculate the values in the data.frame.
#'   Attributes also include the descrepancy between the request and calculated
#'   q and storage. Descrepancy values should be
#'
#' @export
TSZStats = function(TSZs, tau_0, tau_n, storage = NULL, q = NULL, ...,
                    factor = NULL, shape = "powerLaw", MoreArgs = list(),
                    integrateArgs1 = list(), integrateArgs2 = list(), forceNumeric = F,
                    optimizeInterval = NULL, optimizeTol = .Machine$double.eps^0.25)
  {

  if(any(!(sapply(list(tau_0, tau_n, storage, q, factor, optimizeTol, optimizeInterval), mode) %in% c("numeric", "NULL")))) stop("The parameters tau_0, tau_n, storage, q, factor, optimizeTol, optimizeInterval must be numeric or NULL")
  if(any(sapply(list(tau_0, tau_n, storage, q, factor, optimizeTol, ...), length) > 1)) stop("The parameters tau_0, tau_n, storage, q, factor, optimizeTol, and any value passed via ... must have a length() = 1.  If 'fun' requires a multi-value vector to be passed via ..., wrap the vector in a list of length = 1.")
  if(any(!(sapply(list(integrateArgs1, integrateArgs2, MoreArgs), mode) %in% c("list")))) stop("The parameters integrateArgs, and MoreArgs must be lists")

  if (any(sort(TSZs) != TSZs)) stop("TSZs must be sorted ascending.")

  # get the list of arguments
  argList = lapply(match.call()[-1], eval, envir = parent.frame())

  nullParamName = getNullParamName(storage = storage, q = q, ...)

  # if a parameter other than storage or q is specified as NULL, estimate the
  # parameter value and call TSZStats recursively with estimated value and q =
  # NULL
  if(!(nullParamName %in% c("storage", "q"))) {
    if(length(optimizeInterval) != 2) stop("The parameter optimizeInterval must be a vector of length = 2")


#    estimateFunctionValue <- function(tau_0, tau_n, storage = NULL, q = NULL, ..., shape = "powerLaw", MoreArgs = NULL, integrateArgs = NULL, optimizeInterval = c(1,2), optimizeTol = .Machine$double.eps^0.25) {

    # estimate the value
    value =
      estimateFunctionValue(
        tau_0 = tau_0, tau_n = tau_n, storage = storage, q = q, ...,
        shape = shape,
        MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2,
        optimizeInterval = optimizeInterval, optimizeTol = optimizeTol,
        forceNumeric = forceNumeric
      )

    # modify the arg list to use the estimated value while setting q = NULL
    argList[[nullParamName]] = value$minimum

    fraction = min(value$minimum - sort(optimizeInterval)[1], sort(optimizeInterval)[2] - value$minimum)/abs(diff(optimizeInterval))
    if(value$objective > 1 | fraction < 0.01) warning("Minimization routine to estimate '", nullParamName, "' may have failed or estimated value of '", nullParamName, "' was close to the edge of the specified optimization range.  True value of '", nullParamName, "' may be outside of specified optimization range.")
    argList["q"] = list(NULL)
    # call TSZStats again, using modified arguement list
    TSZs = do.call(TSZStats, argList)
    # store the estimate in the attributes list.
    names(value)[1] = "estimate"
    attr(TSZs, nullParamName) = value
    attr(TSZs, "q") = q
    attr(TSZs, "q discrepancy") = sum(TSZs$returning) - q
    attr(TSZs, "storage discrepancy") = sum(TSZs$waterStorage) - storage
  } else {
    # either q or storage is NULL..

    # calculate TSZ breaks
    if(!is.null(factor)) {
      if(factor >= 1) {
        TSZs = getBinBreaks(TSZs, factor, tau_n - tau_0) + tau_0
      } else if(factor == 0) {
        TSZs = getEvenBinBreaks(TSZs, tau_0, tau_n, ..., shape = shape, MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2, forceNumeric = forceNumeric)
      } else {
        stop("Value for 'factor' must be 1.0 or greater.")
      }
    } else {
      if(length(TSZs) < 2) stop("If 'factor' is specified as NULL, 'TSZs' must be a vector that includes the minimum and maximum residence times of interest, plus any residence times of intervening TSZ boundaries. Use a numeric vector of length > 1.")
    }

    # make a data frame with to and from columns -- one row for each bin
    TSZs = data.frame(from = TSZs[1:(length(TSZs)-1)], to = TSZs[2:length(TSZs)])

    # estimate total storage
    if(is.null(storage)) {
      storage = q * autoSolution(suffix = "IntCCDF", shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = tau_0, tau_b = tau_n, ...,
                                 MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2, forceNumeric = forceNumeric)
      argList$storage = list(estimate = storage, objective = NA)
    # estimate q
    } else if(is.null(q)) {
      q = storage/autoSolution(suffix = "IntCCDF", shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = tau_0, tau_b = tau_n, ...,
                               MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2, forceNumeric = forceNumeric)
      argList$q = list(estimate = q, objective = NA)
    } else {
      stop("Unexpected error.  Neither storage nor q is NULL")
    }

    # q times the CDF(tau = from) is flow entering a TSZ
    TSZs$entering = q * autoSolution(suffix = "CCDF", shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = TSZs$from, ..., MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, forceNumeric = forceNumeric)
    # q times the CDF(tau = to) is the flow that moves to the next TSZ
    TSZs$continuing = q * autoSolution(suffix = "CCDF", shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = TSZs$to, ..., MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, forceNumeric = forceNumeric)
    # difference is what upwells
    TSZs$returning = TSZs$entering - TSZs$continuing

    # integration of the CDF from 'from' to 'to' is useful for a few calculations
    integratedCDF = autoSolution(suffix = "IntCCDF", shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = TSZs$from, tau_b = TSZs$to, ..., MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2, forceNumeric = forceNumeric)
    integratedPDF = autoSolution(suffix = "IntPDF", shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = TSZs$from, tau_b = TSZs$to, ..., MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2, forceNumeric = forceNumeric)
    # Storage is q times integration of CDF.
    TSZs$waterStorage = q * integratedCDF
    # Mean water age is water age for a TSZ which equals the various water ages
    # in the bin weighted by the CDF.  So integrate the product of water age and
    # CDF, then divide by integrated CDF.
    TSZs$meanWaterAge = autoSolution(suffix = "IntTau.CCDF", shape = shape,tau_0 = tau_0, tau_n = tau_n,  tau = TSZs$from, tau_b = TSZs$to, ..., MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2, forceNumeric = forceNumeric)/integratedCDF
    # Mean residence time is the mean age of returning water from the bin.
    TSZs$meanResTime = autoSolution(suffix = "IntTau.PDF", shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = TSZs$from, tau_b = TSZs$to, ..., MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2, forceNumeric = forceNumeric)/integratedPDF
  #  TSZs$meanBinResTime = TSZs$meanWaterAge - TSZs$from
  #  TSZs$meanBinFlow = q * numericIntCDF(TSZs$from, TSZs$to, tau_0, tau_n, ..., solution = solvedIntCDF, fun = fun) / (TSZs$to - TSZs$from)
    attributes(TSZs) = c(attributes(TSZs), argList)
    class(TSZs) = c("TSZs", class(TSZs))
    attr(TSZs, "q discrepancy") = sum(TSZs$returning) - q
    attr(TSZs, "storage discrepancy") = sum(TSZs$waterStorage) - storage
  }
  return(TSZs)
}

getEvenBinBreaks = function(nBins, tau_0, tau_n, ..., shape, MoreArgs, integrateArgs1, integrateArgs2, forceNumeric) {

  storageIntegration = autoSolution(suffix = "IntCCDF", shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = tau_0, tau_b = tau_n, ...,
                                    MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2, forceNumeric = forceNumeric)
  target = storageIntegration/nBins
  binError = function(x, ...) {
    storageIntegration = autoSolution(suffix = "IntCCDF", shape = shape, tau_0 = tau_0, tau_n = tau_n, tau = tau_0, tau_b = x, ...,
                                      MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2, forceNumeric = forceNumeric)
    return((target * i - storageIntegration)^2)
  }
  TSZs = numeric(0)
  for (i in 1:(nBins-1)) {
    optimal = optimize(binError, ..., lower = tau_0, upper = tau_n)
    TSZs[i] = optimal$minimum
  }
  return(c(tau_0, TSZs, tau_n))
}

getBinBreaks = function(nbins, factor, maxRT) {
  # Divides a range into bins, where the cumulative bin size doubles with the addition of each bin
  # Returns a vector of bin breakpoints, including the range values as the first and last
  #   breakpoints
  #
  # nbins = number of bins to calculate
  # factor = size of bin relative to next smaller bin (e.g., 2.0 means each bin is double the prior)
  # range = vector with lowest value of first bin and highest value of last bin

  #nbins = 5

  range = c(0,maxRT)

  if(length(maxRT)!=1) stop("'maxRT' must be a single value")
  if(length(nbins)!=1) stop("'nbins' must be a single value")
  nbins = as.integer(nbins)
  range = as.numeric(range)
  factor = as.numeric(factor)
  if (nbins < 1)  stop("Requested number of bins must be a postive number")
  if (factor< 1) stop("'Factor' must be >= 1")
  if(nbins == 1) return(range)

  rangeSize = range[2] - range[1]
  relativeBinSizes = factor^(0:(nbins-1))
  #  relativeBinSizes = relativeBinSizes - c(0, relativeBinSizes[1:(length(relativeBinSizes)-1)])
  firstBinSize = rangeSize/(sum(relativeBinSizes))
  binSizes = relativeBinSizes*rangeSize/(sum(relativeBinSizes))
  binBreaks = c(0, sapply(1:nbins, function(x) sum(binSizes[1:x])))+range[1]
  return(binBreaks)
}

#' @export
estimateFunctionValue <- function(tau_0, tau_n, storage = NULL, q = NULL, ..., shape = "powerLaw", MoreArgs = list(), integrateArgs1 = list(), integrateArgs2 = list(), forceNumeric = F, optimizeInterval = c(1,2), optimizeTol = .Machine$double.eps^0.25) {
  if(any(!(sapply(list(tau_0, tau_n, storage, q, optimizeTol, optimizeInterval), mode) %in% c("numeric", "NULL")))) stop("The parameters tau_0, tau_n, storage, q, optimizeTol, optimizeInterval must be numeric or NULL")
  if(any(sapply(list(tau_0, tau_n, storage, q, optimizeTol), length) > 1)) stop("The parameters tau_0, tau_n, storage, and q should be either NULL or vectors of length = 1")
  if(length(optimizeInterval) != 2) stop("The parameter optimizeInterval must be a vector of length = 2")

  nullParamName = getNullParamName(storage = storage, q = q, ...)
  # error function for optimize; since storage is q * the integral of the CDF,
  # this function squares the difference between storage and q * numericIntCDF
  numericIntCDFArgs =
    list(
      shape = shape, suffix = "IntCCDF",
      tau_0 = tau_0, tau_n = tau_n, tau = tau_0, tau_b = tau_n,
      ...,
      MoreArgs = MoreArgs, integrateArgs1 = integrateArgs1, integrateArgs2 = integrateArgs2, forceNumeric = forceNumeric
    )

    bError = function(x) {
      numericIntCDFArgs[[nullParamName]] <<- x
      return((storage - q * do.call(autoSolution, numericIntCDFArgs))^2)
      #      return((storage - q * numericIntCDF(tau_0, tau_n, tau_0, tau_n, x, ..., solution = solvedIntCDF, fun = fun, integrateArgs = integrateArgs, CDFIntegrateArgs = CDFIntegrateArgs))^2)
    }
  # estimate the value using optimize
  do.call(optimize, c(list(f = bError, interval = optimizeInterval, tol = optimizeTol)))
}

# Returns the name of a parameter in ... that has a NULL value.
getNullParamName <- function(...) {
  params = list(...)
  # check to be sure that one of the key parameters is specified as NULL
  nullParams = sapply(params, is.null)
  if(sum(nullParams) != 1) stop("One (an only one) of the following parameters must be specified as NULL: '", paste(names(params), collapse = "', '"), "'")

  # get the name of the NULL parameter
  names(nullParams)[which(nullParams)]
}

#' @rdname TSZStats
#' @export
hyporheicBins = function(nbins, factor, minRT, maxRT, porosity, hyporheicSize, hyporheicExchange = NULL, b = NULL, integrateArgs = list()) {
  .Deprecated("TSZStats")
  if(is.null(hyporheicExchange)) {
    TSZs = TSZStats(TSZs = nbins, tau_0 = minRT, tau_n = maxRT, storage = hyporheicSize * porosity, alpha = -b, factor = factor)
  } else {
    TSZs = TSZStats(TSZs = nbins, tau_0 = minRT, tau_n = maxRT, storage = hyporheicSize * porosity, q = hyporheicExchange, alpha = NULL, factor = factor, optimizeInterval = c(1, 2))
  }
  saveAttr = attributes(TSZs)

  TSZs$aquiferStorage = TSZs$waterStorage/porosity
  TSZs = TSZs[c("from", "to", "entering", "returning", "continuing", "waterStorage", "aquiferStorage", "meanWaterAge")]
  saveAttr$names = attr(TSZs,"names")
  attributes(TSZs) = saveAttr
  b = ifelse(is.list(saveAttr$alpha), -saveAttr$alpha$estimate, -saveAttr$alpha)

  attr(TSZs, "a") = powerLawPDF(minRT, minRT, maxRT, -b) / powerLaw(minRT, maxRT, -b)
  attr(TSZs, "b") = b
  return(TSZs)
}
