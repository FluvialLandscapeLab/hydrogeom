PDF = function(k, x, xmin){
  -((k+1)/xmin^(k+1)) * x^k
}


# returns the CDF of a PDF of t^k for any time t
# relative to tmin. (e.g. assumes CDF(tmin) = 1)
CDF = function(k, x, xmin) {
  (1/xmin^(k+1)) * x^(k+1)
}

#definite integral of CDF, above
iCDF = function(k, xmin, xmax, CDFxmin = xmin) {
  return((1/CDFxmin^(k+1)) * (xmax^(k+2) - xmin^(k+2))/(k+2))

  #(1/(xmin^(k+1) - xmax^(k+1))) * (((tos-xmax^(k+1))^(k+2)) - ((froms - xmax^(k+1))^(k+2)))/(k+2)
}


storage_Ratios = function(k, xmin, xmax, CDFmin, CDFmax) {

  sratios = (iCDF(k, xmin, xmax, CDFmin) - CDF(k, CDFmax, CDFmin)*(xmax - xmin)) /
            (iCDF(k, CDFmin, CDFmax) - CDF(k, CDFmax, CDFmin)*(CDFmax - CDFmin))

  return(sratios)
}
