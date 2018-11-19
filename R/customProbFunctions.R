powerLaw = function(tau, tau_0, tau_n, alpha) {
  tau^(-alpha) - tau_n^(-alpha)
}

customPDF = function(tau, tau_0, tau_n, ..., solution = powerLawPDF, fun = NULL, returnIntegrateClass = F) {
  checkOrder(c(tau_0 = tau_0, tau = tau, tau_n = tau_n))
  if(is.null(fun)) {
    result = solution(tau = tau, tau_0 = tau_0, tau_n = tau_n, ...)
  } else {
    fun = match.fun(fun)
    # get the function value associated with tau
    funResult = fun(tau, tau_0, tau_n, ...)
    # determine the normalizing constant by integrating the function
    normalizingConstant = integrate(fun, tau_0, tau_n, tau_0 = tau_0, tau_n = tau_n, ...)
    # if the normalizing constant is calculated alright...
    if(normalizingConstant$message == "OK") {
      #check to see if the use wants an S3 "integrate" object.
      if(returnIntegrateClass) {
        # make a copy of the normalizing constant (because it is an S3 integrate object)
        result = normalizingConstant
        # apply the normalizing constant to the function result
        result$value = funResult/normalizingConstant$value
        # carry thru the error in the normalizing constant by dividing the
        # function result by the upper and lower estimate of the normalizing
        # constant.
        interval = funResult/c(normalizingConstant$value + normalizingConstant$abs.error, normalizingConstant$value - normalizingConstant$abs.error)
        result$abs.error = max(abs(result$value - interval))
      } else {
        # just apply the normalizing constant to the function result
        result = funResult/normalizingConstant$value
      }
    }
  }
  return(result)
}

customCDF = function(tau, tau_0, tau_n, ..., solution = powerLawCDF, fun = NULL, returnIntegrateClass = F) {
  checkOrder(c(tau_0 = tau_0, tau = tau, tau_n = tau_n))
  if(is.null(fun)) {
    result = solution(tau = tau, tau_0 = tau_0, tau_n = tau_n, ...)
  }else{
    # integrate the function from tau to tau_n
    integratedFun = integrate(fun, tau, tau_n, tau_0 = tau_0, tau_n = tau_n, ...)
    if(integratedFun$message != "OK") {
      integratedFun$message = paste("While integrating from tau to tau_n:", integratedFun$message)
    } else {
      # store the prior integration for use as the numerator...
      integratedNumerator = integratedFun
      # and calculate the normalizing constant for the denominator...
      integratedFun = integrate(fun, tau_0, tau_n, tau_0 = tau_0, tau_n = tau_n, ...)
      if(integratedFun$message != "OK") {
        integratedFun$message = paste("While integrating from tau_0 to tau_n:", integratedFun$message)
      } else {
        # store the denominator
        integratedDenominator = integratedFun
        if(returnIntegrateClass) {
          # calculates the error intervals for numerator and denominator
          intervals =
            lapply(
              list(integratedNumerator, integratedDenominator),
              function(x) {
                c(x$value + x$abs.error, x$value - x$abs.error)
              }
            )
          # get best estimate of CDF
          integratedFun$value = integratedNumerator$value/integratedDenominator$value
          # calculate maximum error
          integratedFun$abs.error = max(abs(integratedFun$value - sapply(intervals[[1]], function(x) x/intervals[[2]])))
          # return both subdivisions
          integratedFun$subdivisions = c(numerator = integratedNumerator$subdivisions, denominator = integratedDenominator$subdivisions)
          # return both calls
          integratedFun$call = list(numerator = integratedNumerator$call, denominator = integratedDenominator$call)
          result = integratedFun
        } else {
          result = integratedNumerator$value/integratedDenominator$value
        }
      }
    }
  }
  return(result)
}

customIntCDF = function(tau_a, tau_b, tau_0, tau_n, ..., solution = powerLawIntCDF, fun = NULL, returnIntegrateClass = F) {
  checkOrder(c(tau_0 = tau_0, tau_a = tau_a, tau_b = tau_b, tau_n = tau_n))
  if(is.null(fun)) {
    result = solution(tau_a = tau_a, tau_b = tau_b, tau_0 = tau_0, tau_n = tau_n, ...)
  } else {
    # set up some variables
    messages = character(0)
    abs.errors = numeric(0)
    # a function that returns CDF values when a vector of tau values is passed
    # in.  This will allow integration of the CDF.
    justCDFValues = function(tau, bound) {
      # calculate CDF values, returned as a list of integrate classes
      resultList = lapply(tau, customCDF, fun = fun, tau_0 = tau_0, tau_n = tau_n, ..., returnIntegrateClass = T)
      # store all the messages in the parent environment
      messages <<- c(messages, sapply(resultList, `[[`, "message"))
      # extract all the error estimates from the list...
      current.error = sapply(resultList, `[[`, "abs.error")
      # and store them in the parent environments
      abs.errors <<- c(abs.errors, current.error)
      # return the CDF estimates.
      return(sapply(resultList, `[[`, "value") + (current.error * bound))
    }
    # integreate using best estimates of CDF
    if(returnIntegrateClass) {
      integratedFun = integrate(justCDFValues, tau_a, tau_b, bound = 0)
      # integrate using upper error values of CDF
      integratedUpper = integrate(justCDFValues, tau_a, tau_b, bound = 1)
      # integrate using lower error values of CDF
      integratedLower = integrate(justCDFValues, tau_a, tau_b, bound = -1)
      if(any(messages != "OK"))stop("Error when calculating CDF during numerical integration.  Consider using the 'customCDF()' function with various arguments to find out when calculating the CDF of your custom function fails.")
      # report the worst error
      integratedFun$abs.error = max(abs(integratedFun$value - c(integratedUpper$value, integratedLower$value))) + integratedFun$abs.error
      result = integratedFun
    } else {
      result = integrate(justCDFValues, tau_a, tau_b, bound = 0)$value
    }
  }
  return(result)
}

checkOrder = function(aNamedVector) {
  if(!identical(aNamedVector, sort(aNamedVector))) stop("The following condition was not met: ", paste0(names(aNamedVector), collapse = " <= "))
  return(T)
}

useSolution = function(solution, ...) {
  makeIntegrateClass(
    do.call(solution, list(...))
  )
}

makeIntegrateClass = function(value) {
  structure(
    list(
      value = value,
      abs.error = 0,
      subdivisions = NA,
      message = "OK",
      call = NA
    ),
    class = "integrate"
  )
}
