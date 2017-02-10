PDF = function(k, x, xmin){
  -((k+1)/xmin^(k+1)) * x^k
}


# returns the CDF of a PDF of t^k for any time t
# relative to tmin. (e.g. assumes CDF(tmin) = 1)
CDF = function(k, x, xmin) {
  (1/xmin^(k+1)) * x^(k+1)
}

#definite integral of CDF, above
iCDF = function(k, froms, tos, xmin = froms) {
  return (1/xmin^(k+1)) * (tos^(k+2) - froms^(k+2))/(k+2)
}


storage_Ratios = function(k, froms, tos, tmin, tmax) {

  sratios = (iCDF(k, froms, tos) - CDF(k, tmax, tmin)*(tos - froms)) /
            (iCDF(k, tmin, tmax) - CDF(k, tmax, tmin)*(tmax - tmin))

  return(sratios)
}
