function BondPrice3Ps(y,c,t,FaceValue=1)
#BondPrice3Ps   Calculates price from the yield to maturity and the payment
#               stream (arbitrary periods and coupons)
#
#
#
#
#
#  Usage:     [Q,dQ_dy] = BondPrice3Ps(y,c,t[,FaceValue]);
#
#  Input:
#             y          n vector, effective yield to maturity, e.g. 0.07
#             c          scalar or n vector or nxm matrix , coupon rate and face value,
#                        e.g. 0.06 or 0.09/2
#             t          m vector, dates of coupon payments. Principal
#                        is paid at the same time as the last coupon. Dates
#                        should be expressed as fractions of the period.
#             FaceValue  (optional) face value (default: 1)
#
#  Output:    Q      n vector, bond price (eg. 1.01)
#             dQ_dy  n vector, derivative of bond price wrt. yield
#
#
#
#  Note: (1) the following is calculated:
#
#             Q =  c/(1+y)^t(1) + c/(1+y)^t(2) + ... + (1+c)/(1+y)^t(m)
#
#        (2) With semi-annual coupons after 2,8, and 14 months, use
#            t = [0.333;0.333+1;0.333+2] to calculate semiannual yield.
#
#
#
#  Example:        t = (1:2)'
#                  c = [0.09 0.09;
#                       0.10 0.10];
#                  y = [0.0626;0.07];  gives Q = [1.05;1.05]
#
#
#
#   Paul Soderlind (Paul.Soderlind@unisg.ch), June 2009, to Julia 2015
# ----------------------------------------------------------------------------

  y = collect(y)
  t = collect(t)'
  n = length(y)
  m = length(t)

  y = repmat(y,1,m)                         #n -> nxm
  t = repmat(t,n,1)                         #m -> nxm matrix

  if length(c) == 1                         #scalar c
    c = fill(c,(n,m))
  elseif length(c) == n                     #one c for each bond
    c = repmat(c,1,m)
  end
  if length(FaceValue) == 1
    FaceValue = fill(FaceValue,(n,1))
  end

  c[:,end] = c[:,end] + vec(FaceValue)     #add face value to last coupon payment

  cfac = c./((1+y).^t)                     #c/(1+y)^t1 + c/(1+y)^t2 + ...+ c/(1+y)^m
  Q    = sum(cfac,2)                       #theoretical price

  dcfac = cfac .* (-t./(1+y))      #derivative wrt ytm
  dQ_dy = sum(dcfac,2)

  return Q,dQ_dy

end
#-----------------------------------------------------------------------------

