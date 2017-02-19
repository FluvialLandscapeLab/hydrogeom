hyporheicBins = function(nbins, factor, minRT, maxRT, hyporheicExchange, porosity, hyporheicSize) {

  hyporheicSize = hyporheicSize*porosity

  bError = function(b) {
    sizeEstimate = hyporheicExchange * integrate(function(x) PDF(x, b, minRT, maxRT) * (x - minRT), minRT, maxRT)$value  #storageFactor(b, minRT, maxRT, minRT, maxRT)
    return((hyporheicSize - sizeEstimate)^2)
  }

  # find b associated with user-specified storage:exchange ratio
  b = optimize(f = bError, interval = c(-5, -1))$minimum
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
getHyporheicBinStats3 = function(nbins, factor, minRT, maxRT, aquiferVol, hyporheicQ, porosity, infiniteAquifer = T) {

  if(infiniteAquifer) {
    integrateQ = function(b,s,e,m) (e^(b+1) - s^(b+1))/(b+1)
    integrateV = function(b,s,e,m) (e^(b+2) - s^(b+2))/(b+2)
  } else {
    integrateQ = function(b,s,e,m) (e^(b+1)/(b+1)-e^(b+2)/((b+2)*m)) - (s^(b+1)/(b+1)-s^(b+2)/((b+2)*m))
    integrateV = function(b,s,e,m) (e^(b+2)/(b+2)-e^(b+3)/((b+3)*m)) - (s^(b+2)/(b+2)-s^(b+3)/((b+3)*m))
  }

  # find b where user-specified hyporheicQ/aquiferWaterStorage = integrationOfQ(minRT-->maxRT)/integrationOfS(minRT-->maxRT)
  b = optimize(f = bError, hyporheicQ = hyporheicQ, aquiferVol = aquiferVol, porosity = porosity, minRT = minRT, maxRT = maxRT, integrateQ = integrateQ, integrateV = integrateV, lower = -1.999999999, upper = -1.00000001)$minimum

  rtBinBreaks = getBinBreaks(nbins, factor, maxRT - minRT) + minRT
  L = length(rtBinBreaks)
  binRTStarts = rtBinBreaks[1:(L-1)]
  binRTEnds = rtBinBreaks[2:L]
  binRTs = binRTEnds - binRTStarts

  binReturnFlows = mapply(integrateQ, b, binRTStarts, binRTEnds, maxRT)
  binVols = mapply(integrateV, b, binRTStarts, binRTEnds, maxRT) / porosity

  binMeanFlowPathResTime = mapply(optimize, startRT = binRTStarts, endRT = binRTEnds, lower = binRTStarts, upper = binRTEnds, MoreArgs = list(f = meanRTError, b = b, maxRT = maxRT, integrateQ = integrateQ, infiniteAquifer = infiniteAquifer))
  binMeanBinResTime = unlist(binMeanFlowPathResTime["minimum",]) - binRTStarts
  a = hyporheicQ / sum(binReturnFlows)
  binReturnFlows = a * binReturnFlows
  binVols = a * binVols

  rtdStats = list(a=a, b=b)
  results = data.frame(
    from = binRTStarts,
    to = binRTEnds,
    volume = binVols,
    meanFlowPathResTime = unlist(binMeanFlowPathResTime["minimum",]),
    meanBinResTime = binMeanBinResTime,
    returnFlow = binReturnFlows
  )
  attributes(results) = c(attributes(results), rtdStats)
  return(results)
}
