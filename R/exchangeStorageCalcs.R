exchangeRateFromGeom = function(tau_0, tau_n, alpha, A, V, phi, solution = powerLawIntCDF, fun = NULL) {
  exchangeRateFromStorage(tau_0, tau_n, alpha, V*phi/A, fun, solution)
}

exchangeRateFromStorage = function(tau_0, tau_n, s, ..., fun = NULL, solution = powerLawIntCDF) {
  if(is.null(fun)) {
    result = makeIntegrateClass(s / useSolution(solution, tau_a = tau_0, tau_b = tau_n, tau_0 = tau_0, tau_n = tau_n, ...)$value)
  } else {
    result = customIntCDF(tau_0, tau_n, tau_0, tau_n, ..., fun = fun)
    xRate = s / result$value
    result$abs.error = max(abs(xRate - (s / (result$value + (result$abs.error * c(1, -1))))))
    result$value = xRate
  }
  return(result)
}

storageFromExchangeRate = function(tau_0, tau_n, q, ..., solution = powerLawIntCDF, fun = NULL) {
  if(is.null(fun)) {
    result = makeIntegrateClass(q * useSolution(solution, tau_a = tau_0, tau_b = tau_n, tau_0 = tau_0, tau_n = tau_n, ...)$value)
  } else {
    result = customIntCDF(tau_0, tau_n, tau_0, tau_n, ..., fun = fun)
    xRate = q * result$value
    result$abs.error = max(abs(xRate - (q * (result$value + (result$abs.error * c(1, -1))))))
    result$value = xRate
  }
  return(result)
}

estimateAlpha = function(tau_0, tau_n, s, q, interval = c(1.1, 1.9), solution = powerLawIntCDF, fun = NULL) {
  result = estimateFunctionParam(tau_0, tau_n, s, q, interval = interval, solution = solution, fun = fun)
  names(result)[1] = "alpha"
  return(result)
}

estimateFunctionParam = function(tau_0, tau_n, s, q, ..., interval, solution=NULL, fun=NULL) {

  errorFunction = function(x, bound, ...) {
    sEst = storageFromExchangeRate(tau_0 = tau_0, tau_n = tau_n, q = q, x, ..., solution = solution, fun = fun)
    if(sEst$message != "OK") stop("Error message from 'storageFromExchange' function.  You can troubleshoot by passing the same values of tau_0, tau_n, s, q, solution, fun, etc. to 'storageFromExchange' manually.  Error message was: '", error$message,"'")
    error = (s - (sEst$value + (sEst$abs.error * bound)))^2
    return(error)
  }

  estimate = optimize(errorFunction, interval, bound = 0, ...)
  upperEstimate = optimize(errorFunction, interval, bound = 1, ...)
  lowerEstimate = optimize(errorFunction, interval, bound = -1, ...)
  estimate$abs.error = max(abs(estimate$minimum - c(upperEstimate$minimum, lowerEstimate$minimum)))
  names(estimate)[1] = "estimate"
  return(estimate)
}
