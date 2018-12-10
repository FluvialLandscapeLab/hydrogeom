powerLaw = function(tau, tau_n, alpha) {
  if(any(tau > tau_n)) stop("Values of tau must be greater than values of tau_n")
  tau^(-alpha) - tau_n^(-alpha)
}

checkCustomFunction = function(fun, tau, tau_0, tau_n, ...){
  integrate(fun, lower = tau, upper = tau_n, tau_n = tau_n, ...)
}

getNormalizingConstant = function(fun, tau_0, tau_n, ..., integrateArgs = list()) {
  paramList = apply(cbind(tau = tau_0, tau_0, tau_n, ...), 1, function(x) c(list(fun = fun), as.list(x), integrateArgs))
  sapply(paramList, function(x) do.call(checkCustomFunction, x)$value)
}

customPDF = function(tau, tau_0, tau_n, ..., solution = powerLawPDF, fun = NULL, integrateArgs = list()) {
  checkOrder(tau_0, tau, tau_n)
  if(is.null(fun)) {
    solution = match.fun(solution)
    result = solution(tau = tau, tau_0 = tau_0, tau_n = tau_n, ...)
  } else {
    # get the function values
    fun = match.fun(fun)
    funResult = fun(tau, tau_n, ...)
    result = funResult/getNormalizingConstant(fun, tau_0, tau_n, ..., integrateArgs = integrateArgs)
  }
  return(result)
}

customCDF = function(tau, tau_0, tau_n, ..., solution = powerLawCDF, fun = NULL, integrateArgs = list()) {
  checkOrder(tau_0, tau, tau_n)
  if(is.null(fun)) {
    solution = match.fun(solution)
    result = solution(tau = tau, tau_0 = tau_0, tau_n = tau_n, ...)
  }else{
    # integrate the custom function from tau to tau_n
    paramList = apply(cbind(tau, tau_0, tau_n, ...), 1, function(x) c(list(fun = fun), as.list(x), integrateArgs))
    integratedFun = sapply(paramList, function(x) do.call(checkCustomFunction, x)$value)
    result = integratedFun/getNormalizingConstant(fun, tau_0, tau_n, ..., integrateArgs = integrateArgs)
  }
  return(result)
}

customIntCDF = function(tau_a, tau_b, tau_0, tau_n, ..., solution = powerLawIntCDF, fun = NULL, integrateArgs = list(), CDFIntegrateArgs = list()) {
  checkOrder(tau_0, tau_a, tau_b, tau_n)
  if(is.null(fun)) {
    result = solution(tau_a = tau_a, tau_b = tau_b, tau_0 = tau_0, tau_n = tau_n, ...)
  } else {
    paramList =
      apply(
        cbind(lower = tau_a, upper = tau_b, tau_0, tau_n, ...),
        1,
        function(x) c(list(f = customCDF), as.list(x), list(fun = fun, integrateArgs = CDFIntegrateArgs), integrateArgs)
      )
    result = sapply(paramList, function(x) do.call(integrate, x)$value)
  }
  return(result)
}

customIntTau.CDF = function(tau_a, tau_b, tau_0, tau_n, ..., solution = powerLawIntTau.CDF, fun = NULL, integrateArgs = list(), CDFIntegrateArgs = list()) {
  checkOrder(tau_0, tau_a, tau_b, tau_n)
  if(is.null(fun)) {
    result = solution(tau_a = tau_a, tau_b = tau_b, tau_0 = tau_0, tau_n = tau_n, ...)
  } else {
    tau.CDF = function(tau, tau_0, tau_n, ...) {
      tau * customCDF(tau, tau_0, tau_n, ..., fun = fun, integrateArgs = CDFIntegrateArgs)
    }
    paramList =
      apply(
        cbind(lower = tau_a, upper = tau_b, tau_0, tau_n, ...),
        1,
        function(x) c(list(f = tau.CDF), as.list(x), integrateArgs)
      )
    result = sapply(paramList, function(x) do.call(integrate, x)$value)
  }
  return(result)
}

checkOrder = function(...) {
  aMatrix = cbind(...)
  sortedRows = t(apply(aMatrix, 1, sort))
  if(!all(aMatrix == sortedRows)) stop("Requested values of 'tau' were out of order (not ascending) or outside the specified minimum ('tau_0') or maximum ('tau_n') integration range.")
  return(T)
}
