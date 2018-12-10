
hyporheicBinsOld = function(nbins, factor, minRT, maxRT, porosity, hyporheicSize, hyporheicExchange = NULL, b = NULL) {

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
