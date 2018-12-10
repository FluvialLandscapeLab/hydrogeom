powerLawPDF = function(tau, tau_0, tau_n, alpha) {
  numerator = (tau^(-alpha) - tau_n^(-alpha)) * (1-alpha)
  denominator = tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha)
  return(numerator/denominator)
}

powerLawCDF = function(tau, tau_0, tau_n, alpha) {
  numerator = tau_n^(1-alpha) - tau^(1-alpha) - tau_n^(-alpha)*(tau_n - tau)*(1-alpha)
  denominator = tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha)
  return(numerator/denominator)
}

powerLawIntCDF = function(tau_a, tau_b, tau_0, tau_n, alpha) {
  numerator =  (1-alpha)*tau_n^(-alpha)*(tau_b^2-tau_a^2)/2 - (tau_b^(2-alpha) - tau_a^(2-alpha))/(2-alpha) - (tau_b-tau_a)*(-alpha)*tau_n^(1-alpha)
  denominator = ((tau_n^(1-alpha) - tau_0^(1-alpha)) - tau_n^(-alpha) * (tau_n - tau_0) * (1-alpha))
  return(numerator/denominator)
}

powerLawIntTau.CDF = function(tau_a, tau_b, tau_0, tau_n, alpha) {
  intTau = (tau_b^2 - tau_a^2)/2
  intTau2 = (tau_b^3 - tau_a^3)/3
  numerator = intTau * tau_n^(1-alpha) - (tau_b^(3-alpha)-tau_a^(3-alpha))/(3-alpha) - (1-alpha)*tau_n^(-alpha)*(tau_n * intTau - intTau2)
  denominator = (tau_n^(1-alpha) - tau_0^(1-alpha) - tau_n^(-alpha)*(tau_n - tau_0)*(1-alpha))
  return(numerator/denominator)
}
