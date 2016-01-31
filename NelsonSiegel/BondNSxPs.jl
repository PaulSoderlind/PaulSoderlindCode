function BondNSxPs(m,b0,b1,b2,tau,b3=0,tau2=1)
#BondNS2xPs    Extended Nelson and Siegel (1987) spot rate, forward rate, and discount function
#
#
#
#
#  Usage:    (s,f,d) = BondNSxPs(m,b0,b1,b2,tau[,b3[,tau2]])
#
#  Input:    m         NxK matrix  times to maturity in years (don't use less than 1e-8 or so)
#            b0        scalar, beta0 parameter in NS
#            b1        scalar, beta1 parameter in NS
#            b2        scalar, beta2 parameter in NS
#            tau       scalar, tau   parameter in NS
#            b3        (optional) scalar, beta2 parameter in extended NS [0]
#            tau2      (optional) scalar, tau2  parameter in extended NS [1]
#
#  Output:   s         NxK matrix,  spot rate (continously compounded)
#            f         NxK matrix,  forward rate (continously compounded)
#            d         NxK matrix,  discount function
#
#
#
#  Note:     For more details, see
#            (a) Svensson (1995), "Estimating Forward Interest Rates with
#                the Extended Nelson & Siegel Method," Quarterly Review,
#                Sveriges Riksbank, 1995:3, 13-26.
#            (b) Soderlind and Svensson (1997), "New Techniques to Extract
#                Market Expectations from Financial Instruments,"
#                Journal of Monetary Economics 40, 383-429.
#
#
#
#  Paul.Soderlind@unisg.ch, April 2002, to Julia 2015
#----------------------------------------------------------------------------*/

                                                 #forward rate
  f =  b0 + b1*exp(-m/tau) + b2*(m/tau).*exp(-m/tau) + b3*(m/tau2).*exp(-m/tau2)

  s =  b0 + b1*(1 - exp(-m/tau))./(m/tau) +
            b2*((1 - exp(-m/tau)) ./(m/tau)  - exp(-m/tau)) +
            b3*((1 - exp(-m/tau2))./(m/tau2) - exp(-m/tau2))        #spot rate

  d = exp(-s.*m)             #discount function

  return s,f,d

end
#------------------------------------------------------------------------------
