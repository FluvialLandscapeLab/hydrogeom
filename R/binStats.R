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

TSZStats = function(TSZs, tau_0, tau_n, storage = NULL, q = NULL, ..., factor = NULL, solvedCDF = powerLawCDF, solvedIntCDF = powerLawIntCDF, fun = NULL, interval = NULL) {
  arguments = as.list(match.call())
  arguments = sapply(arguments[2:length(arguments)], eval)
  argumentTypes = sapply(arguments, mode)
  invalidArguments = !(argumentTypes %in% c("numeric", "NULL"))
  if(any(invalidArguments)) {
    stop(
      "All arguments passed to 'TSZStats' must be numeric, however '",
      paste0(
        names(arguments[invalidArguments]),
        "' was type '",
        sapply(argumentTypes[invalidArguments], mode),
        collapse = "'; "),
      "'."
    )
  }

  if (any(sort(TSZs) != TSZs)) stop("TSZs must be sorted ascending.")

  if(sum(is.null(storage), is.null(q)) == 2) stop("At least one of the parameters 'storage' or 'q' must be specified")

  # getStorage = function(...) {
  #   storageIntegration = customIntCDF(tau_0, tau_n, tau_0, tau_n, ..., solution = solvedIntCDF, fun = fun)
  #   if(storageIntegration$message != "OK") {
  #     stop("Numerical integration error.  Message was '", storageIntegration$message, "'")
  #   }
  #   return(storageIntegration)
  # }

  if(is.null(storage)) {
    storage = q * customIntCDF(tau_0, tau_n, tau_0, tau_n, ..., solution = solvedIntCDF, fun = fun)$value
  } else if(is.null(q)) {
    q = storage/customIntCDF(tau_0, tau_n, tau_0, tau_n, ..., solution = solvedIntCDF, fun = fun)$value
  } else {
    # try to use optimize() to estimate a value from the integration function
    bError = function(x, ...) {
#      storageIntegration = getStorage(x, ...)
#      return((storage - storageIntegration$value*q)^2)
      return((storage - q * customIntCDF(tau_0, tau_n, tau_0, tau_n, x, ..., solution = solvedIntCDF, fun = fun)$value)^2)
    }
    value = optimize(f = bError, interval = interval)$minimum
  }

  if(!is.null(factor)) {
    if(factor >= 1) {
      TSZs = getBinBreaks(TSZs, factor, tau_n - tau_0) + tau_0
    } else if(factor == 0) {
      TSZs = getEvenAreaBinBreaks(TSZs, tau_0, tau_n, ..., solution = solvedIntCDF, fun = fun)
    } else {
      stop("Value for 'factor' must be 1.0 or greater.")
    }
  } else {
    if(length(TSZs) < 2) stop("The 'TSZs' vector must include the minimum and maximum residence times of interest, plus any residence times of intervening TSZ boundaries. Use a numeric vector of length > 1.")
  }

  TSZs = data.frame(from = TSZs[1:(length(TSZs)-1)], to = TSZs[2:length(TSZs)])

  flowList =
    lapply(
      TSZs,
      sapply,
      function(f) {
        flow = customCDF(tau = f, tau_0 = tau_0, tau_n = tau_n, ..., solution = solvedCDF, fun = fun)
        if(flow$message != "OK") stop("Integration error in calculating CDF. Message: ", flow$message)
        return(q * flow$value)
      }
    )

  TSZs$entering = flowList$from
  TSZs$returning = flowList$from - flowList$to
  TSZs$continuing = flowList$to

  TSZs$waterStorage =
    mapply(
      function(ta, tb, ...) {
        storage = customIntCDF(tau_a = ta, tau_b = tb, ...)
        if(storage$message != "OK") stop("Error integrating CDF. Message:", storage$message)
        return(q * storage$value)
      },
      ta = TSZs$from,
      tb = TSZs$to,
      MoreArgs =
        list(
          tau_0 = tau_0,
          tau_n = tau_n,
          ...,
          solution = solvedIntCDF,
          fun = fun
        )
    )

  return(TSZs)
}


getEvenAreaBinBreaks = function(nBins, tau_0, tau_n, ..., solution = powerLawIntCDF, fun = NULL) {

  storageIntegration = customIntCDF(tau_0, tau_n, tau_0, tau_n, ..., solution = solution, fun = fun)
  if(storageIntegration$message != "OK") {
    stop("Numerical integration error.  Message was '", storageIntegration$message, "'")
  }
  target = storageIntegration$value/nBins
  binError = function(x, ...) {
    storageIntegration = customIntCDF(tau_0, x, tau_0, tau_n, ..., solution = solution, fun = fun)
    if(storageIntegration$message != "OK") {
      stop("Numerical integration error.  Message was '", storageIntegration$message, "'")
    }
    return((target * i - storageIntegration$value)^2)
  }
  TSZs = numeric(0)
  for (i in 1:(nBins-1)) {
    optimal = optimize(binError, ..., lower = tau_0, upper = tau_n)
    TSZs[i] = optimal$minimum
  }
  return(c(tau_0, TSZs, tau_n))
}







hyporheicBins = function(nbins, factor, minRT, maxRT, porosity, hyporheicSize, hyporheicExchange = NULL, b = NULL) {

  hyporheicSize = hyporheicSize*porosity

    bError = function(b) {
      sizeEstimate = hyporheicExchange * integrate(function(x) PDF(x, b, minRT, maxRT) * (x - minRT), minRT, maxRT)$value  #storageFactor(b, minRT, maxRT, minRT, maxRT)
      return((hyporheicSize - sizeEstimate)^2)
    }
  # find b associated with user-specified storage:exchange ratio
    if(is.null(b)) {
      b = optimize(f = bError, interval = c(-5, -1))$minimum
    } else {
      hyporheicExchange = hyporheicSize / integrate(function(x) PDF(x, b, minRT, maxRT) * (x - minRT), minRT, maxRT)$value
    }
    a = hyporheicExchange * PDFaValue(b, minRT, maxRT)

  binBreaks = getBinBreaks(nbins, factor, maxRT - minRT) + minRT
  from = 1:nbins
  to = from + 1

  inFlow = hyporheicExchange * CDF(binBreaks[from], b, minRT, maxRT)
  continuedFlow = hyporheicExchange * CDF(binBreaks[to], b, minRT, maxRT)
  returnFlow = inFlow - continuedFlow

  storage = hyporheicExchange * storageFactor(b, minRT, maxRT, binBreaks[from], binBreaks[to])

  meanBinStatList = meanBinStats(b, minRT, maxRT, binBreaks[from], binBreaks[to])

  rtdStats = list(a=a, b=b)

  results = data.frame(
    from = binBreaks[from],
    to = binBreaks[to],
    entering = inFlow,
    returning = returnFlow,
    continuing = continuedFlow,
    waterStorage = storage,
    aquiferStorage = storage/porosity,
    meanWaterAge = meanBinStatList$meanWaterAge,
    meanBinResTime = meanBinStatList$meanWaterAge - binBreaks[from],
    meanBinFlow = hyporheicExchange * meanBinStatList$meanBinFlow
  )
  attributes(results) = c(attributes(results), rtdStats)
  return(results)


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


meanRTError = function(meanRT, b, startRT, endRT, maxRT, integrateQ, infiniteAquifer = T) {
  if(infiniteAquifer) {
    getQ = function(b, x, m) x^b
  } else {
    getQ = function(b, x, m) (x^b) * (m-x) / m
  }
  return((getQ(b, meanRT, maxRT) - integrateQ(b, startRT, endRT, maxRT) / (endRT - startRT))^2)
}

# nbins = number of partitions to use in HZ
# factor = controls relative size of bins; 2 is a doubling of residence time for each bin.
# minRT = minimum residence time to consider in the HZ
# maxRT = maximum residence time to consider in the HZ
# aquiferVol = volume of HZ, including sediment
# hyporheicQ = gross exchange rate across stream bed for flow paths with total residence time between minRT and maxRT
# porosity = porosity (fraction void space) of aquifer.
# infiniteAqulfer = Assume aquifer is infinite if T, assume linear decrease to aquifer volume = 0 at maxRT if F.
# getHyporheicBinStats3 = function(nbins, factor, minRT, maxRT, aquiferVol, hyporheicQ, porosity, infiniteAquifer = T) {
#
#   if(infiniteAquifer) {
#     integrateQ = function(b,s,e,m) (e^(b+1) - s^(b+1))/(b+1)
#     integrateV = function(b,s,e,m) (e^(b+2) - s^(b+2))/(b+2)
#   } else {
#     integrateQ = function(b,s,e,m) (e^(b+1)/(b+1)-e^(b+2)/((b+2)*m)) - (s^(b+1)/(b+1)-s^(b+2)/((b+2)*m))
#     integrateV = function(b,s,e,m) (e^(b+2)/(b+2)-e^(b+3)/((b+3)*m)) - (s^(b+2)/(b+2)-s^(b+3)/((b+3)*m))
#   }
#
#   # find b where user-specified hyporheicQ/aquiferWaterStorage = integrationOfQ(minRT-->maxRT)/integrationOfS(minRT-->maxRT)
#   b = optimize(f = bError, hyporheicQ = hyporheicQ, aquiferVol = aquiferVol, porosity = porosity, minRT = minRT, maxRT = maxRT, integrateQ = integrateQ, integrateV = integrateV, lower = -1.999999999, upper = -1.00000001)$minimum
#
#   rtBinBreaks = getBinBreaks(nbins, factor, maxRT - minRT) + minRT
#   L = length(rtBinBreaks)
#   binRTStarts = rtBinBreaks[1:(L-1)]
#   binRTEnds = rtBinBreaks[2:L]
#   binRTs = binRTEnds - binRTStarts
#
#   binReturnFlows = mapply(integrateQ, b, binRTStarts, binRTEnds, maxRT)
#   binVols = mapply(integrateV, b, binRTStarts, binRTEnds, maxRT) / porosity
#
#   binMeanFlowPathResTime = mapply(optimize, startRT = binRTStarts, endRT = binRTEnds, lower = binRTStarts, upper = binRTEnds, MoreArgs = list(f = meanRTError, b = b, maxRT = maxRT, integrateQ = integrateQ, infiniteAquifer = infiniteAquifer))
#   binMeanBinResTime = unlist(binMeanFlowPathResTime["minimum",]) - binRTStarts
#   a = hyporheicQ / sum(binReturnFlows)
#   binReturnFlows = a * binReturnFlows
#   binVols = a * binVols
#
#   rtdStats = list(a=a, b=b)
#   results = data.frame(
#     from = binRTStarts,
#     to = binRTEnds,
#     volume = binVols,
#     meanFlowPathResTime = unlist(binMeanFlowPathResTime["minimum",]),
#     meanBinResTime = binMeanBinResTime,
#     returnFlow = binReturnFlows
#   )
#   attributes(results) = c(attributes(results), rtdStats)
#   return(results)
# }
