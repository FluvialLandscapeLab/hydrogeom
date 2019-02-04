#' Solved probability functions
#'
#' \code{hydrogeom} can provide transient storage zone statistics based on
#' solutions to six integrations of a shape file. These functions are examples
#' of solutions for the \code{\link{powerLaw}} and \code{\link{exponent}} shape
#' functions.
#'
#' Calculation of transient storage zone (TSZ) statistics (see
#' \code{\link{TSZStats}}) requires six functions:
#'
#' 1) the probability density function (PDF) of the curve shape (i.e., the
#' probability that a water molecule will exit the hyporheic zone with a
#' residence time of \code{tau})
#'
#' 2) the complementary cumulative distribution function (CCDF) (i.e., the
#' probablity that a water molecule will still be in the hyporheic zone at
#' residence time \code{tau})
#'
#' 3) the definite integral of the PDF.
#'
#' 4) the definite integral of the CCDF.
#'
#' 5) the definite integral of \code{tau}*PDF.
#'
#' 6) the definite integral of \code{tau}*CCDF.
#'
#' These values can be derived with numeric integration of a shape function
#' (e.g., see \code{\link{powerLaw}}), but numerical stablity can be an issue
#' (e.g., see \code{\link{checkShapeFunction}}).  Providing solutions for the
#' six integrations, above, may yield more reliable stability in the
#' calculations.
#'
#' The solutions are published in Poole et al. (in Press) and their
#' implementation in this package can be viewed by printing the source code for
#' each of the following functions (e.g., type print(<functionName>) in the R
#' console).
#'
#' Solution functions for other shapes can be provided by the user.  However, they must
#' follow the naming conventions of the concatination of the shape name and
#' the specific suffixes: "PDF", "CCDF", "IntPDF", "IntCCDF", "IntTau.PDF", and "IntTau.CCDF".
#' For instance, for a shape called "foo", the function names would be "fooPDF", "fooCCDF", "fooIntPDF", etc.
#' Also, the PDF and CDF solutions must have the signature \code{function(tau, tau_0, tau_n, ...)} while the
#' remaining solutions must have the signature \code{function(tau, tau_0, tau_n, ...)} where "..." are
#' additional parameters required by the shape function.
#'
#' @param tau A residence time for which the PDF or CCDF value is desired.
#' @param tau_a,tau_b The lower and upper range for the definite integral of the PDF or CCDF.
#' @param tau_0,tau_n The total range of \code{tau} values considered.
#' @param alpha The exponent of the power law (this value is negated within the
#'   functions; entering 1.5 yields a power law with an exponent of -1.5).
#' @param sigma The exchange rate of groundwater with an exponential
#'   distribution.  (This value is negated in the calculations to simulate the
#'   decay of water volume with residence time.)
#' @return Returns a vector of numeric values.
#'
#'   PowerLawPDF and exponentPDF return the PDF value for \code{tau}.  This
#'   provides the probability density that a water molecule entering the
#'   hyporheic zone at time == \code{tau_0} will return to the stream at time ==
#'   \code{tau}.
#'
#'   PowerLawCCDF and exponentCCDF return the CCDF value for \code{tau}.  This
#'   provides the probablity that a water molecule entering the hyporehic zone
#'   at time == \code{tau_0} will still be in the hyporheic zone at time
#'   \code{tau} (e.g., will upwell after time == \code{tau})
#'
#'   PowerLawIntPDF and exponentIntPDF return the definite integral from \code{tau_a} to
#'   \code{tau_b} of the PDF.
#'
#'   PowerLawIntCCDF and exponentIntCCDF return the difinite integral from \code{tau_a} to
#'   \code{tau_b} of the CCDF.
#'
#'   PowerLawIntTau.PDF and exponentIntTau.PDF returns the definite integral from \code{tau_a} to
#'   \code{tau_b} of \code{tau}*PDF.
#'
#'   PowerLawIntTau.CCDF and exponentIntTau.CCDF return the difinite integral from \code{tau_a} to
#'   \code{tau_b} of \code{tau}*CCDF.
#' @export
powerLawPDF = function(tau, tau_0, tau_n, alpha) {
  numerator = (tau^(-alpha) - tau_n^(-alpha)) * (1-alpha)
  denominator = tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha)
  return(numerator/denominator)
}

#' @rdname powerLawPDF
#' @export
powerLawCCDF = function(tau, tau_0, tau_n, alpha) {
  numerator = tau_n^(1-alpha) - tau^(1-alpha) - tau_n^(-alpha)*(tau_n - tau)*(1-alpha)
  denominator = tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha)
  return(numerator/denominator)
}

#' @rdname powerLawPDF
#' @export
powerLawIntPDF = function(tau_a, tau_b, tau_0, tau_n, alpha) {
  numerator = tau_b^(1-alpha) - tau_a^(1-alpha) - tau_n^(-alpha)*(tau_b - tau_a)*(1-alpha)
  denominator = tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha)
  return(numerator/denominator)
}

#' @rdname powerLawPDF
#' @export
powerLawIntCCDF = function(tau_a, tau_b, tau_0, tau_n, alpha) {
  numerator =  (1-alpha)*tau_n^(-alpha)*(tau_b^2-tau_a^2)/2 - (tau_b^(2-alpha) - tau_a^(2-alpha))/(2-alpha) - (tau_b-tau_a)*(-alpha)*tau_n^(1-alpha)
  denominator = ((tau_n^(1-alpha) - tau_0^(1-alpha)) - tau_n^(-alpha) * (tau_n - tau_0) * (1-alpha))
  return(numerator/denominator)
}

#' @rdname powerLawPDF
#' @export
powerLawIntTau.CCDF = function(tau_a, tau_b, tau_0, tau_n, alpha) {
  intTau = (tau_b^2 - tau_a^2)/2
  intTau2 = (tau_b^3 - tau_a^3)/3
  numerator = intTau * tau_n^(1-alpha) - (tau_b^(3-alpha)-tau_a^(3-alpha))/(3-alpha) - (1-alpha)*tau_n^(-alpha)*(tau_n * intTau - intTau2)
  denominator = (tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha))
  return(numerator/denominator)
}

#' @rdname powerLawPDF
#' @export
powerLawIntTau.PDF = function(tau_a, tau_b, tau_0, tau_n, alpha) {
  normalizingConstant = (1-alpha) / (tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha))
  #integral of tau^(1-alpha)
  integral1 = (tau_b^(2-alpha) - tau_a^(2-alpha)) / (2-alpha)
  #integral of tau*tau_m^(1-alpha)
  integral2 = 0.5 * (tau_b^2 - tau_a^2) * tau_n^(-alpha)
  return(normalizingConstant * (integral1 - integral2))
}

#' @rdname powerLawPDF
#' @export
exponentPDF = function(tau, tau_0, tau_n, sigma) {
  if(any(tau > tau_n)) stop("Values of tau must be less than values of tau_n")
  normalizingConstant = sigma / (exp(-sigma * tau_n) * (exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1))
  return((exp(-sigma*tau) - exp(-sigma * tau_n)) * normalizingConstant)
}

#' @rdname powerLawPDF
#' @export
exponentCCDF = function(tau, tau_0, tau_n, sigma) {
  numerator =   exp(-sigma * (tau  -tau_n)) - sigma*(tau_n-tau  ) - 1
  denominator = exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1
  return(numerator/denominator)
}

#' @rdname powerLawPDF
#' @export
exponentIntPDF = function(tau_a, tau_b, tau_0, tau_n, sigma) {
  numerator = sigma*(tau_a - tau_b) + exp(-sigma*(tau_a - tau_n)) - exp(-sigma*(tau_b - tau_n))
  denominator = exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1
  return(numerator/denominator)
}

#' @rdname powerLawPDF
#' @export
exponentIntCCDF = function(tau_a, tau_b, tau_0, tau_n, sigma) {
  numerator1 = 0.5*sigma*(tau_a - tau_b)*(tau_a + tau_b - 2*tau_n)
  numerator2 = (1/sigma) * (exp(-sigma*(tau_a - tau_n)) - exp(-sigma*(tau_b - tau_n)))

  numerator = (-numerator1 + numerator2 + tau_a - tau_b)
  denominator = (exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1)

  return(numerator/denominator)
}

#' @rdname powerLawPDF
#' @export
exponentIntTau.PDF = function(tau_a, tau_b, tau_0, tau_n, sigma) {
  numerator1 = (tau_a*sigma + 1)*exp((tau_a - tau_n)*-sigma)
  numerator2 = (tau_b*sigma + 1)*exp((tau_b - tau_n)*-sigma)
  numerator3 = sigma^2 * (tau_b^2-tau_a^2)/2
  denominator = sigma * (exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1)

  #  denominator = sigma *(1 - (tau_n-tau_0)*sigma + exp(-sigma*(tau_0-tau_n)))
  return((numerator1 - numerator2 - numerator3)/denominator)
}

#' @rdname powerLawPDF
#' @export
exponentIntTau.CCDF = function(tau_a, tau_b, tau_0, tau_n, sigma) {
  numerator1a = (sigma*tau_a + 1)*exp(-sigma*(tau_a - tau_n))
  numerator1b = (sigma*tau_b + 1)*exp(-sigma*(tau_b - tau_n))
  numerator1 = (1/sigma^2)*(numerator1a - numerator1b)
  numerator2a = tau_a^2*(2*tau_a - 3*tau_n)
  numerator2b = tau_b^2*(3*tau_n - 2*tau_b)
  numerator2 = (sigma/6)*(numerator2a + numerator2b)
  numerator = numerator1 - numerator2 - (tau_b^2 - tau_a^2)/2
  denominator = (exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1)
  return(numerator/denominator)
}
