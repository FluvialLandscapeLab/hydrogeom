PDF = function(k, x, xMin){
  -((k+1)/xMin^(k+1)) * x^k
}


# returns the CDF of a PDF of t^k for any time t
# relative to tmin. (e.g. assumes CDF(tmin) = 1)
CDF = function(k, x, xMin) {
  (1/xMin^(k+1)) * x^(k+1)
}

#definite integral of CDF, above
iCDF = function(k, xMin, xMax, CDFxMin = xMin) {
  return((1/CDFxMin^(k+1)) * (xMax^(k+2) - xMin^(k+2))/(k+2))

  #(1/(xMin^(k+1) - xMax^(k+1))) * (((tos-xMax^(k+1))^(k+2)) - ((froms - xMax^(k+1))^(k+2)))/(k+2)
}


storage_Ratios = function(k, xMin, xMax, CDFmin, CDFmax) {

  sratios = (iCDF(k, xMin, xMax, CDFmin) - CDF(k, CDFmax, CDFmin)*(xMax - xMin)) /
            (iCDF(k, CDFmin, CDFmax) - CDF(k, CDFmax, CDFmin)*(CDFmax - CDFmin))

  return(sratios)
}

meanWaterAge = function(k, xMin, xMax, CDFmin) {
  pMin = CDF(k, xMax, CDFmin)
  pMax = CDF(k, xMin, CDFmin)

  kstar = (k+2)/(k+1)
  return(CDFmin * (k+1) * (pMax^kstar - pMin^kstar) / ((k+2)*(pMax - pMin)))
}
