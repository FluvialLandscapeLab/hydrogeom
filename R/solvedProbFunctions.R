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
#' probability that a water molecule will still be in the hyporheic zone at
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
#' (e.g., see \code{\link{powerLaw}}), but numerical stability can be an issue
#' (e.g., see \code{\link{checkShapeFunction}}).  Providing solutions for the
#' six integrations, above, may yield more reliable stability in the
#' calculations.
#'
#' The solutions are published in Poole et al. (in Press) and their
#' implementation in this package can be viewed by printing the source code for
#' each of the following functions (e.g., type print(<functionName>) in the R
#' console).
#'
#' Solution functions for other shapes can be provided by the user.  However,
#' they must follow the naming conventions of the concatination of the shape
#' name and the specific suffixes: "PDF", "CCDF", "IntPDF", "IntCCDF",
#' "IntTau.PDF", and "IntTau.CCDF". For instance, for a shape called "foo", the
#' function names would be "fooPDF", "fooCCDF", "fooIntPDF", etc. Also, the PDF
#' and CDF solutions must have the signature \code{function(tau, tau_0, tau_n,
#' ...)} while the remaining solutions must have the signature
#' \code{function(tau_a, tau_b, tau_0, tau_n, ...)} where "..." are additional
#' parameters required by the shape function.
#'
#' @param tau A residence time for which the PDF or CCDF value is desired.
#' @param tau_a,tau_b The lower and upper range for the definite integral of the
#'   PDF or CCDF.
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
#'   provides the probability that a water molecule entering the hyporheic zone
#'   at time == \code{tau_0} will still be in the hyporheic zone at time
#'   \code{tau} (e.g., will upwell after time == \code{tau})
#'
#'   PowerLawIntPDF and exponentIntPDF return the definite integral from
#'   \code{tau_a} to \code{tau_b} of the PDF.
#'
#'   PowerLawIntCCDF and exponentIntCCDF return the definite integral from
#'   \code{tau_a} to \code{tau_b} of the CCDF.
#'
#'   PowerLawIntTau.PDF and exponentIntTau.PDF returns the definite integral
#'   from \code{tau_a} to \code{tau_b} of \code{tau}*PDF.
#'
#'   PowerLawIntTau.CCDF and exponentIntTau.CCDF return the definite integral
#'   from \code{tau_a} to \code{tau_b} of \code{tau}*CCDF.
#'
#' @details
#' The power law functions include special handling for limiting cases where
#' the general formula would involve division by zero:
#' \itemize{
#'   \item \code{powerLawPDF}, \code{powerLawCCDF}, \code{powerLawIntPDF},
#'     \code{powerLawIntCCDF}: Handle the alpha = 1 limiting case
#'   \item \code{powerLawIntCCDF}, \code{powerLawIntTau.PDF}: Handle the
#'     alpha = 2 limiting case
#'   \item \code{powerLawIntTau.CCDF}: Handles the alpha = 3 limiting case
#' }
#'
#' All functions validate input parameters and will stop with an error if:
#' \itemize{
#'   \item For power law functions: \code{tau_0 <= 0} (power law requires
#'     \code{tau_0 > 0} to avoid division by zero when computing
#'     \code{tau_0^(-alpha)})
#'   \item For PDF/CCDF functions: \code{tau} is outside the range
#'     \code{[tau_0, tau_n]}
#'   \item For integration functions: \code{tau_a < tau_0}, \code{tau_b > tau_n},
#'     or \code{tau_a > tau_b}
#' }
#' @export
powerLawPDF = function(tau, tau_0, tau_n, alpha) {
  if (any(tau_0 <= 0)) {
    stop("tau_0 must be greater than 0 for power law functions")
  }
  if (any(tau < tau_0 | tau > tau_n)) {
    stop("Values of tau must satisfy tau_0 <= tau <= tau_n")
  }

  # General formula (produces NaN/Inf when alpha = 1)
  result = (tau^(-alpha) - tau_n^(-alpha)) /
    ((tau_n^(1-alpha) - tau_0^(1-alpha))/(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0))

  # Handle alpha = 1 edge case; standard calc fails with alpha = 1
  alt_calc = abs(alpha - 1) < 1e-10

  if(any(alt_calc)) {
    # Recycle alt_calc to match result length (in case alpha, which determines alt_calc, is shorter than result)
    alt_calc = alt_calc == rep(TRUE, length(result))
    # Compute limiting formula and replace affected elements
    alt_result = (tau^(-1) - tau_n^(-1)) /
      (log(tau_n/tau_0) - tau_n^(-1)*(tau_n - tau_0))
    result[alt_calc] = alt_result[alt_calc]
  }

  return(result)
}
# powerLawPDF = function(tau, tau_0, tau_n, alpha) {
#   numerator = (tau^(-alpha) - tau_n^(-alpha)) * (1-alpha)
#   denominator = tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha)
#   return(numerator/denominator)
# }

#' @rdname powerLawPDF
#' @export
powerLawCCDF = function(tau, tau_0, tau_n, alpha) {
  if (any(tau_0 <= 0)) {
    stop("tau_0 must be greater than 0 for power law functions")
  }
  if (any(tau < tau_0 | tau > tau_n)) {
    stop("Values of tau must satisfy tau_0 <= tau <= tau_n")
  }

  result = (tau_n^(1-alpha) - tau^(1-alpha) - tau_n^(-alpha)*(tau_n - tau)*(1-alpha)) /
    (tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha))

  # See powerLawPDF comments for explanation of alt_calc
  alt_calc = abs(alpha - 1) < 1e-10
  if(any(alt_calc)) {
    alt_calc = alt_calc == rep(TRUE, length(result))
    alt_result = (log(tau_n/tau) - tau_n^(-1)*(tau_n - tau)) /
      (log(tau_n/tau_0) - tau_n^(-1)*(tau_n - tau_0))
    result[alt_calc] = alt_result[alt_calc]
  }

  return(result)
}

# powerLawCCDF = function(tau, tau_0, tau_n, alpha) {
#   numerator = tau_n^(1-alpha) - tau^(1-alpha) - tau_n^(-alpha)*(tau_n - tau)*(1-alpha)
#   denominator = tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha)
#   return(numerator/denominator)
# }

#' @rdname powerLawPDF
#' @export
powerLawIntPDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, alpha) {
  if (any(tau_0 <= 0)) {
    stop("tau_0 must be greater than 0 for power law functions")
  }
  if (any(tau_a < tau_0 | tau_b > tau_n | tau_a > tau_b)) {
    stop("Values must satisfy tau_0 <= tau_a <= tau_b <= tau_n")
  }

  result = (tau_b^(1-alpha) - tau_a^(1-alpha) - tau_n^(-alpha)*(tau_b - tau_a)*(1-alpha)) /
    (tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha))

  # See powerLawPDF comments for explanation of alt_calc
  alt_calc = abs(alpha - 1) < 1e-10
  if(any(alt_calc)) {
    alt_calc = alt_calc == rep(TRUE, length(result))
    alt_result = (log(tau_b/tau_a) - tau_n^(-1)*(tau_b - tau_a)) /
      (log(tau_n/tau_0) - tau_n^(-1)*(tau_n - tau_0))
    result[alt_calc] = alt_result[alt_calc]
  }

  return(result)
}

# powerLawIntPDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, alpha) {
#   numerator = tau_b^(1-alpha) - tau_a^(1-alpha) - tau_n^(-alpha)*(tau_b - tau_a)*(1-alpha)
#   denominator = tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha)
#   return(numerator/denominator)
# }

#' @rdname powerLawPDF
#' @export
powerLawIntCCDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, alpha) {
  if (any(tau_0 <= 0)) {
    stop("tau_0 must be greater than 0 for power law functions")
  }
  if (any(tau_a < tau_0 | tau_b > tau_n | tau_a > tau_b)) {
    stop("Values must satisfy tau_0 <= tau_a <= tau_b <= tau_n")
  }

  result = ((1-alpha)*tau_n^(-alpha)*(tau_b^2 - tau_a^2)/2 -
      (tau_b^(2-alpha) - tau_a^(2-alpha))/(2-alpha) +
      alpha*tau_n^(1-alpha)*(tau_b - tau_a)) /
    (tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha))

  # See powerLawPDF comments for explanation of alt_calc
  alt_calc = abs(alpha - 1) < 1e-10
  if(any(alt_calc)) {
    alt_calc = alt_calc == rep(TRUE, length(result))
    # th
    alt_result = (tau_b*log(tau_n/tau_b) - tau_a*log(tau_n/tau_a) + (tau_b^2 - tau_a^2)/(2*tau_n)) /
      (log(tau_n/tau_0) - tau_n^(-1)*(tau_n - tau_0)) * rep(1, length(alpha)) # force recycling to to length of alpha
    result[alt_calc] = alt_result[alt_calc]
  }

  alt_calc = abs(alpha - 2) < 1e-10
  if(any(alt_calc)) {
    alt_calc = alt_calc == rep(TRUE, length(result))
    alt_result = (-tau_n^(-2)*(tau_b^2 - tau_a^2)/2 - log(tau_b/tau_a) + 2*tau_n^(-1)*(tau_b - tau_a)) /
      (tau_n^(-1) - tau_0^(-1) + tau_n^(-2)*(tau_n - tau_0))  * rep(1, length(alpha)) # force recycling to to length of alpha
    result[alt_calc] = alt_result[alt_calc]
  }

  return(result)
}

# powerLawIntCCDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, alpha) {
#   numerator =  (1-alpha)*tau_n^(-alpha)*(tau_b^2-tau_a^2)/2 - (tau_b^(2-alpha) - tau_a^(2-alpha))/(2-alpha) - (tau_b-tau_a)*(-alpha)*tau_n^(1-alpha)
#   denominator = ((tau_n^(1-alpha) - tau_0^(1-alpha)) - tau_n^(-alpha) * (tau_n - tau_0) * (1-alpha))
#   return(numerator/denominator)
# }

#' @rdname powerLawPDF
#' @export
powerLawIntTau.CCDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, alpha) {
  if (any(tau_0 <= 0)) {
    stop("tau_0 must be greater than 0 for power law functions")
  }
  if (any(tau_a < tau_0 | tau_b > tau_n | tau_a > tau_b)) {
    stop("Values must satisfy tau_0 <= tau_a <= tau_b <= tau_n")
  }

  intTau = (tau_b^2 - tau_a^2)/2
  intTau2 = (tau_b^3 - tau_a^3)/3

  result = (intTau * tau_n^(1-alpha) - (tau_b^(3-alpha) - tau_a^(3-alpha))/(3-alpha) -
      (1-alpha)*tau_n^(-alpha)*(tau_n * intTau - intTau2)) /
    (tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha))

  # See powerLawPDF comments for explanation of alt_calc
  alt_calc = abs(alpha - 3) < 1e-10
  if(any(alt_calc)) {
    alt_calc = alt_calc == rep(TRUE, length(result))
    alt_result = (intTau * tau_n^(-2) - log(tau_b/tau_a) + 2*tau_n^(-3)*(tau_n * intTau - intTau2)) /
      (tau_n^(-2) - tau_0^(-2) + 2*tau_n^(-3)*(tau_n - tau_0)) * rep(1, length(alpha)) # force recycling to length of alpha
    result[alt_calc] = alt_result[alt_calc]
  }

  return(result)
}

# powerLawIntTau.CCDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, alpha) {
#   intTau = (tau_b^2 - tau_a^2)/2
#   intTau2 = (tau_b^3 - tau_a^3)/3
#   numerator = intTau * tau_n^(1-alpha) - (tau_b^(3-alpha)-tau_a^(3-alpha))/(3-alpha) - (1-alpha)*tau_n^(-alpha)*(tau_n * intTau - intTau2)
#   denominator = (tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha))
#   return(numerator/denominator)
# }

#' @rdname powerLawPDF
#' @export
powerLawIntTau.PDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, alpha) {
  if (any(tau_0 <= 0)) {
    stop("tau_0 must be greater than 0 for power law functions")
  }
  if (any(tau_a < tau_0 | tau_b > tau_n | tau_a > tau_b)) {
    stop("Values must satisfy tau_0 <= tau_a <= tau_b <= tau_n")
  }

  result = ((1-alpha)*(tau_b^(2-alpha) - tau_a^(2-alpha))/(2-alpha) -
      tau_n^(-alpha)*(tau_b^2 - tau_a^2)/2*(1-alpha)) /
    (tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha))

  # See powerLawPDF comments for explanation of alt_calc
  alt_calc = abs(alpha - 2) < 1e-10
  if(any(alt_calc)) {
    alt_calc = alt_calc == rep(TRUE, length(result))
    alt_result = (log(tau_b/tau_a) - tau_n^(-2)*(tau_b^2 - tau_a^2)/2) /
      (tau_0^(-1) - tau_n^(-1) - tau_n^(-2)*(tau_n - tau_0)) * rep(1, length(alpha)) # force recycling to length of alpha
    result[alt_calc] = alt_result[alt_calc]
  }

  return(result)
}

# powerLawIntTau.PDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, alpha) {
#   normalizingConstant = (1-alpha) / (tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha))
#   #integral of tau^(1-alpha)
#   integral1 = (tau_b^(2-alpha) - tau_a^(2-alpha)) / (2-alpha)
#   #integral of tau*tau_m^(1-alpha)
#   integral2 = 0.5 * (tau_b^2 - tau_a^2) * tau_n^(-alpha)
#   return(normalizingConstant * (integral1 - integral2))
# }

# NOTE on the exponent* functions: every exponential below is written as
# exp(-sigma * (x - tau_0)) for some x in [tau_0, tau_n].  That is the general
# formula with exp(-sigma * tau_0) factored out of numerator and denominator.
# Because x >= tau_0, every exponent is <= 0 and every exp() is in [0, 1], so
# nothing can overflow regardless of sigma, tau_0, or tau_n.  When
# sigma * (tau_n - tau_0) is large, exp(-sigma * (tau_n - tau_0)) underflows to
# 0 and the formulas reduce to the infinite-tail (asymptotic) solutions on their
# own.  No special-case branches are needed.

#' @rdname powerLawPDF
#' @export
exponentPDF <- function(tau, tau_0, tau_n, sigma) {
  if (any(tau < tau_0 | tau > tau_n)) {
    stop("Values of tau must satisfy tau_0 <= tau <= tau_n")
  }

  (sigma * (exp(-sigma * (tau - tau_0)) - exp(-sigma * (tau_n - tau_0)))) /
    (1 - (1 + sigma * (tau_n - tau_0)) * exp(-sigma * (tau_n - tau_0)))
}
# exponentPDF = function(tau, tau_0, tau_n, sigma) {
#   if(any(tau > tau_n)) stop("Values of tau must be less than values of tau_n")
#   normalizingConstant = sigma / (exp(-sigma * tau_n) * (exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1))
#   return((exp(-sigma*tau) - exp(-sigma * tau_n)) * normalizingConstant)
# }

#' @rdname powerLawPDF
#' @export
exponentCCDF <- function(tau, tau_0, tau_n, sigma) {
  if (any(tau < tau_0 | tau > tau_n)) {
    stop("Values of tau must satisfy tau_0 <= tau <= tau_n")
  }

  (exp(-sigma * (tau - tau_0)) - (1 + sigma * (tau_n - tau)) * exp(-sigma * (tau_n - tau_0))) /
    (1 - (1 + sigma * (tau_n - tau_0)) * exp(-sigma * (tau_n - tau_0)))
}

# exponentCCDF = function(tau, tau_0, tau_n, sigma) {
#   numerator =   exp(-sigma * (tau  -tau_n)) - sigma*(tau_n-tau  ) - 1
#   denominator = exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1
#   return(numerator/denominator)
# }

#' @rdname powerLawPDF
#' @export
exponentIntPDF <- function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, sigma) {
  if (any(tau_a < tau_0 | tau_b > tau_n | tau_a > tau_b)) {
    stop("Values must satisfy tau_0 <= tau_a <= tau_b <= tau_n")
  }

  (exp(-sigma * (tau_a - tau_0)) - exp(-sigma * (tau_b - tau_0)) -
      sigma * (tau_b - tau_a) * exp(-sigma * (tau_n - tau_0))) /
    (1 - (1 + sigma * (tau_n - tau_0)) * exp(-sigma * (tau_n - tau_0)))
}
# exponentIntPDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, sigma) {
#   numerator = sigma*(tau_a - tau_b) + exp(-sigma*(tau_a - tau_n)) - exp(-sigma*(tau_b - tau_n))
#   denominator = exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1
#   return(numerator/denominator)
# }

#' @rdname powerLawPDF
#' @export
exponentIntCCDF <- function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, sigma) {
  if (any(tau_a < tau_0 | tau_b > tau_n | tau_a > tau_b)) {
    stop("Values must satisfy tau_0 <= tau_a <= tau_b <= tau_n")
  }

  ((1/sigma) * (exp(-sigma * (tau_a - tau_0)) - exp(-sigma * (tau_b - tau_0))) +
      (tau_b - tau_a) * (sigma/2 * (tau_a + tau_b - 2*tau_n) - 1) * exp(-sigma * (tau_n - tau_0))) /
    (1 - (1 + sigma * (tau_n - tau_0)) * exp(-sigma * (tau_n - tau_0)))
}
# exponentIntCCDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, sigma) {
#   numerator1 = 0.5*sigma*(tau_a - tau_b)*(tau_a + tau_b - 2*tau_n)
#   numerator2 = (1/sigma) * (exp(-sigma*(tau_a - tau_n)) - exp(-sigma*(tau_b - tau_n)))
#
#   numerator = (-numerator1 + numerator2 + tau_a - tau_b)
#   denominator = (exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1)
#
#   return(numerator/denominator)
# }

#' @rdname powerLawPDF
#' @export
exponentIntTau.PDF <- function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, sigma) {
  if (any(tau_a < tau_0 | tau_b > tau_n | tau_a > tau_b)) {
    stop("Values must satisfy tau_0 <= tau_a <= tau_b <= tau_n")
  }

  ((tau_a + 1/sigma) * exp(-sigma * (tau_a - tau_0)) -
      (tau_b + 1/sigma) * exp(-sigma * (tau_b - tau_0)) -
      sigma * (tau_b^2 - tau_a^2) / 2 * exp(-sigma * (tau_n - tau_0))) /
    (1 - (1 + sigma * (tau_n - tau_0)) * exp(-sigma * (tau_n - tau_0)))
}
# exponentIntTau.PDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, sigma) {
#   numerator1 = (tau_a*sigma + 1)*exp((tau_a - tau_n)*-sigma)
#   numerator2 = (tau_b*sigma + 1)*exp((tau_b - tau_n)*-sigma)
#   numerator3 = sigma^2 * (tau_b^2-tau_a^2)/2
#   denominator = sigma * (exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1)
#
#   #  denominator = sigma *(1 - (tau_n-tau_0)*sigma + exp(-sigma*(tau_0-tau_n)))
#   return((numerator1 - numerator2 - numerator3)/denominator)
# }

#' @rdname powerLawPDF
#' @export
exponentIntTau.CCDF <- function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, sigma) {
  if (any(tau_a < tau_0 | tau_b > tau_n | tau_a > tau_b)) {
    stop("Values must satisfy tau_0 <= tau_a <= tau_b <= tau_n")
  }

  ((1/sigma) * ((tau_a + 1/sigma) * exp(-sigma * (tau_a - tau_0)) -
                   (tau_b + 1/sigma) * exp(-sigma * (tau_b - tau_0))) +
        (sigma * (tau_b^3 - tau_a^3) / 3 -
         (1 + sigma * tau_n) * (tau_b^2 - tau_a^2) / 2) * exp(-sigma * (tau_n - tau_0))) /
    (1 - (1 + sigma * (tau_n - tau_0)) * exp(-sigma * (tau_n - tau_0)))
}
# exponentIntTau.CCDF = function(tau_a, tau_b, tau_0 = tau_a, tau_n = tau_b, sigma) {
#   numerator1a = (sigma*tau_a + 1)*exp(-sigma*(tau_a - tau_n))
#   numerator1b = (sigma*tau_b + 1)*exp(-sigma*(tau_b - tau_n))
#   numerator1 = (1/sigma^2)*(numerator1a - numerator1b)
#   numerator2a = tau_a^2*(2*tau_a - 3*tau_n)
#   numerator2b = tau_b^2*(3*tau_n - 2*tau_b)
#   numerator2 = (sigma/6)*(numerator2a + numerator2b)
#   numerator = numerator1 - numerator2 - (tau_b^2 - tau_a^2)/2
#   denominator = (exp(-sigma * (tau_0-tau_n)) - sigma*(tau_n-tau_0) - 1)
#   return(numerator/denominator)
# }
