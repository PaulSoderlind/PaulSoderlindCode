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
#             y          nx1 vector, effective yield to maturity, e.g. 0.07
#             c          scalar or n vector or mxn matrix , coupon rate and face value,
#                        e.g. 0.06 or 0.09/2
#             t          mx1 vector, dates of coupon payments. Principal
#                        is paid at the same time as the last coupon. Dates
#                        should be expressed as fractions of the period.
#             FaceValue  (optional) face value (default: 1)
#
#  Output:    Q      nx1 vector, bond price (eg. 1.01)
#             dQ_dy  nx1 vector, derivative of bond price wrt. yield
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
#  Example:        t = (1:2)';
#                  c = [0.09,0.10;
#                       0.09,0.10];
#                  y = [0.0626;0.07];  gives Q = [1.05;1.05]
#
#
#
#   Paul Soderlind (Paul.Soderlind@unisg.ch), June 2009, to Julia 2015
# ----------------------------------------------------------------------------

  t = collect(t)
  y = collect(y)
  m = length(t)
  n = length(y)

  y = repmat(y',m,1)                      #nx1 -> mxn
  t = repmat(t ,1,n)                      #mx1 -> mxn matrix

  if length(c) == 1                         #scalar c
    c = fill(c,(m,n))
  elseif length(c) == n                     #one c for each bond
    c = repmat(vec(c)',m,1)
  end
  if length(FaceValue) == 1
    FaceValue = fill(FaceValue,(1,n))
  end

  c[end,:] = c[end,:] + vec(FaceValue)'     #add face value to last coupon payment

  cfac = c./((1+y).^t)                     #c/(1+y)^t1 + c/(1+y)^t2 + ...+ c/(1+y)^m
  Q    = sum(cfac,1)                       #theoretical price
  Q    = vec(Q)

  dcfac = cfac .* (-t./(1+y))      #derivative wrt ytm
  dQ_dy = sum(dcfac,1)
  dQ_dy = vec(dQ_dy)

  return Q,dQ_dy

end
#-----------------------------------------------------------------------------

