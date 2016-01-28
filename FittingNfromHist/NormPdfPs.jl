function NormPdfPs(y,mu=0,s2=1)
#NormPdfPs    Returns pdf value of normally distributed (univariate) variables
#
#
#
#  Usage:    pdfy = NormPdfPs(y[,mu[,s2]])
#
#  Input:    y          mxn matrix of values of y which is N(mu,s^2)
#            mu         scalar or 1xn, mean of y, defaults to 0
#            s2         scalar or 1xn, variance of y, defaults to 1
#
#
#  Output:   pdfy       mxn matrix of pdf values
#
#
#
#
#  Paul.Soderlind@unisg.ch, Nov 2002, to Julia Jan 2016
#------------------------------------------------------------------------------

  s = sqrt(s2)

  z    = (y .- mu)./s
  pdfy = exp(-0.5*z.^2)./(sqrt(2*pi)*s)   #pdf of y ~ N(mu,s^2)

  return pdfy

end
#------------------------------------------------------------------------------

