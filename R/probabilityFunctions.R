PDF = function(x, k, n, m) {
  bad = x<n | x>m
  result = (x^k - m^k) * PDFaValue(k, n, m)
  result[bad] = NA
  return(result)
}

PDFaValue = function(k, n, m) {
  bad = n>=m
  result = (k+1) / (m^(k+1) - n^(k+1) - (m - n) * m^k * (k+1))
  result[bad] = NA
  return(result)
}

definiteIntegralPDF = function(x, k, m) {
  bad = x<=0 | x>m
  result = (m^(k+1) - x^(k+1)) / (k+1) - (m-x) * m^k
  result[bad] = NA
  return(result)
}

CDF = function(x, k, n, m) {
  bad = x<n | x>m
  result = definiteIntegralPDF(x, k, m) / definiteIntegralPDF(n, k, m)
  result[bad] = NA
  return(result)
}

storageFactor = function(k, n, m, b1, b2, numerical = F){

  #!!!!! NOTE: This is just the integral of the CDF from b1 to b2!!!!
#  integratedCDFNumerator = ((b1^(k+2) - b2^k * ((b2-b1) * k * ((k+1)*(b2-b1)/2 - b1) + b1^2))/(k+2))


  integratedCDFNumerator =  (k+1)*m^k*(b2^2-b1^2)/2 - (b2^(k+2) - b1^(k+2))/(k+2) - (b2-b1)*k*m^(k+1)
  CDFDenom = ((m^(k+1) - n^(k+1)) - m^k * (m - n) * (k+1))
  #integrate (x^k - m^k) * (k+1) * (x - n) / (m^(k+1) - n^(k+1) - (m - n) * m^k * (k+1)) dx from n to m for k=-1.5, n = 10, m=100
  bad = n>=m

  if(numerical == T) {
    sF = function(x, b) PDF(x, k, n, m) * (x - b)
    result = rep(numeric(0), length(b1))
    for(i in 1:length(b1)) {
      result[i] = integrate(sF, b1[i], b2[i], b = b1[i])$value + CDF(b2[i], k, n, m) * (b2[i] - b1[i])
    }
  } else {
    result = integratedCDFNumerator/CDFDenom
  }
#  part1 = b1^2 * (k*m^k*(3*m-b1) + b1*(-m^k + 6*b1^k / (6+5*k+k^2)))
#  part2 = b2 * (k*m^k * (-3*m*b2 + 2*b2^2 + 6*m*b1 - 3*b2*b1) + (b2*((6+5*k+k^2)*m^k*(2*b2-3*b1)-6*b2^k*((2+k)*b2-(3+k)*b1)))/((2+k)*(3+k)))
#  result = (-part1+part2)/(6*(k+1))
#  result = (k+1) * (n*m^k*b2 - m^k*b2^2/2 - n*b2^(k+1)/(k+1) + b2^(k+2)/(k+2) - n*m^k*b1 + m^k*b1^2/2 + n*b1^(k+1)/(k+1) - b1^(k+2)/(k+2)) / (m^(k+1) - n^(k+1) - (m-n)*m^k*(k+1))
#  result = (m^k * ( (m-n)*k * (k*(m-n)/2 + ((m-n)/2 - n)) + n^2 ) - n^(k+2) ) / ( (k+2) * ( m^k * ((m-n)*k - n) + n^(k+1) ) )
  result[bad] = NA
  return(result)
}

meanBinStats = function(k, n, m, b1, b2) {
  bad = b1>=b2
  result = list(meanWaterAge = rep(numeric(0), length(b1)), meanBinFlow = rep(numeric(0), length(b1)))
  for (i in 1:length(b1)) {
    integratedCDF = integrate(function(x) CDF(x, k, n, m), b1[i], b2[i])$value
#    lessArea = CDF(b2[i], k, n, m) * (b2[i]-b1[i])
#    diffCDF = CDF(b1[i], hyporhk, n, m) - CDF(b2[i], k, n, m)
#    result$meanWaterAge[i] = b1[i] + (integratedCDF - lessArea) / diffCDF
    result$meanWaterAge[i] = integrate(function(x) CDF(x, k, n, m) * x, b1[i], b2[i])$value / integratedCDF
    result$meanBinFlow[i] = integratedCDF / (b2[i] - b1[i])
  }
  result$meanWaterAge[bad] = NA
  result$meanBinFlow[bad] = NA
  return(result)
}

meanFlow = function(k, n, m, b1, b2) {
  bad = b1>=b2
  result = rep(numeric(0), length(b1))
  for (i in 1:length(b1)) {
    integratedCDF = integrate(function(x) CDF(x, k, n, m), b1[i], b2[i])$value
  }
}


linPDF = function(x, slope, int, n, m) {
  lpdf = function(x) slope * x + int
  return(lpdf(x)/integrate(lpdf, n, m)$value)
}

linCDF = function(x, slope, int, n, m) {
  return(sapply(x, function(x) integrate(ilpdf, x, m)$value))
}

ilpdf = function(x) {
  linPDF(x, slope, int, n, m)
}
