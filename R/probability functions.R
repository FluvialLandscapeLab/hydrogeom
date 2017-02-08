# returns the fraction of flux remaining in the aquifer at time t given.
# Assumes that fraction of flux = 1.0 at t = tmin
CDF_Flux = function(k, t, tmin) {
  t^(k+1)/tmin^(k+1)
}

