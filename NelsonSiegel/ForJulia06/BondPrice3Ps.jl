function BondPrice3Ps(y,c0,t,FaceValue=1)
#BondPrice3Ps   Calculates price from the yield to maturity and the payment
#               stream (arbitrary periods and coupons)
#
#
#  Usage:     (Q,dQ_dy) = BondPrice3Ps(y,c0,t[,FaceValue]);
#
#  Input:
#             y          n vector, effective yield to maturity, e.g. 0.07
#             c0         1xn vector or mxn matrix, coupon rate and face value,
#                        e.g. 0.06 or 0.09/2
#             t          m vector, dates of coupon payments. Principal
#                        is paid at the same time as the last coupon. Dates
#                        should be expressed as fractions of the period.
#             FaceValue  (optional) scalar or n vector, face value (default: 1)
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
#  Example:        t = 1:2
#                  c = [0.09 0.09;
#                       0.10 0.10];
#                  y = [0.0626 0.07];  gives Q = [1.05;1.05]
#
#
#
#   Paul Soderlind (Paul.Soderlind@unisg.ch), June 2009, to Julia 2015
# ----------------------------------------------------------------------------

  n = length(y)                          #n is determined by data on y
  if length(c0) == 1 && n > 1
    error("length(c0)=1 requires n=1")
  end

  m = length(t)
  y = repeat(vecPs(y)',outer=(m,1))      #n -> mxn
  t = repeat(vecPs(t),outer=(1,n))       #m -> mxn matrix

  c = deepcopy(c0)                          #breaks link with c0 argument
  if length(c) == 1                         #scalar c, single bond
    c = fill(c[1],(m,1))
  elseif size(c,1,2) == (1,n)
    c = repeat(c,outer=(m,1))
  end
  if length(FaceValue) == 1
    c[end,:] = c[end:end,:] + FaceValue[1]    #add face value to last coupon payment
  else
    c[end,:] = c[end:end,:] + vecPs(FaceValue)'
  end
  #println("c and t")
  #printmat(c)
  #printmat(t)

  cfac  = c./((1+y).^t)                     #c/(1+y)^t1 + c/(1+y)^t2 + ...+ c/(1+y)^m
  Q     = vec(sum(cfac,1))                  #theoretical price

  dcfac = cfac .* (-t./(1+y))               #derivative wrt ytm
  dQ_dy = vec(sum(dcfac,1))
  if n == 1                                 #1x1 array to scalar
    (Q,dQ_dy) = (Q[1],dQ_dy[1])
  end

  return Q,dQ_dy

end
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------

function BondYieldToMat3Ps(Q,c,t,Method=1,yLH=[-0.1;0.5],tol=1e-7,FaceValue=1)
#BondYieldToMat3Ps    Calculates yield to maturity from bond price (several methods available).
#                     Can also be used for general IRR calculations. Works with arbitrary
#                     coupon periods.
#
#
#  Usage:     y = BondYieldToMat3Ps(Q,c,t[,Method[,yLH[,tol[,FaceValue]]]])
#
#
#  Input:     Q          n vector, bond prices (for instance, 1.01)
#             c          1xn vector or mxn matrix, c rate, e.g. 0.06, or 0.09/2
#             t          m vector, dates of coupon payments. Principal
#                        is paid at the same time as the last c. Dates
#                        should be expressed as fractions of the period.
#                        Example: coupons after 2,8, and 14 months with semi-annual coupons->
#                        t = [0.333;0.333+1;0.333+2]
#             Method     optional, 1: Newton-Raphson; 2: bisection [1]
#             yLH        optional, if Method==1: scalar, yLH[1] is initial guess of roots
#                                  if Method==2: 2x1 vector, yLH[1] is lower boundary and yLH[2] upper boundary
#             tol        optional, scalar, convergence criterion for y, [1e-7]
#             FaceValue  optional, scalar or n vector, face value, default [1]
#
#  Output:    y          n vector, effective yield to maturity (per period)
#
#
#
#  Note:      (1)  Method 1 (Newton-Raphson) is pretty fast. The Newton-Raphson iterations
#                  are based on
#                  Q = F(y0) + J0*Dy, where Dx = y1-y0 and J0 is the Jacobian at y0.
#                  We choose Dy to make this exact, that is as  (Q - F(y0))./J0 = Dy
#
#             (2)  Method 2 (bisection) is a bit slow, but very robust.
#
#             (3)  The formula is
#                   Q =  c/y1^t(1) + c/y1^t(2) + ... + (1+c)/y1^t(m)#
#                   with y1 = 1+y, where y is the yield to maturity,
#                   c is per period coupon and y1 the gross yield (per period)
#                   to maturity
#
#             (4) To use for IRR calculations, set Q=0 and FaceValue=0
#
#             (5) For Newton-Raphson, we could use (c+(FaceValue-B)/n)/((FaceValue+B)/2)
#                 as a good starting value
#
#             (6) This code is slower than needed since the t and c arrays
#                 are recreated in each iteration. Good enough for occasional
#                 computations, but should be changed in case execution time
#                 becomes important.
#
#
#  Calls on:  BondPrice3Ps
#
#
#  Paul Soderlind (Paul.Soderlind@unisg.ch), June 2009, to Julia 2015
#------------------------------------------------------------------------------

  Q = vecPs(Q)                    #to make sure is an n vector (as F0 is)
  n = length(Q)

  if Method == 1                   #Newton-Raphson

    y = fill(yLH[1],n)
    Dy = 1e+198
    while maximum(abs.(Dy)) > tol
      (F0,J0) = BondPrice3Ps(y,c,t,FaceValue)  #F(y0) and derivative, nx1
      Dy = (Q - F0)./J0                         #nx1
      y  = y + Dy
    end

  elseif Method == 2                   #bisection

    yL = fill(yLH[1],n)                #lower boundary for yield
    yH = fill(yLH[2],n)                #upper boundary
    if any(yL .> yH)
      warn("Lower bound greater than upper bound")
    end
    F0, = BondPrice3Ps(yL,c,t,FaceValue)     #create starting value, so [yL,yH] brackets the roots
    while any(F0 .< Q)                       #as long as theoretical < actual price
      vv = F0 .< Q
      yL[vv] = yL[vv] - 0.01
      #println(yL)
      F0, = BondPrice3Ps(yL,c,t,FaceValue)
    end
    F0, = BondPrice3Ps(yH,c,t,FaceValue)
    while any(F0 .> Q)
      vv = F0 .> Q
      yH[vv] = yH[vv] + 0.01
      F0, = BondPrice3Ps(yH,c,t,FaceValue)
    end

    y = (yL + yH)/2
    while any((yH-yL) .> tol)              #iteration loop
      y = (yL + yH)/2                      #mid point for yield
      F0, = BondPrice3Ps(y,c,t,FaceValue)  #price at guessed yield
      vvH     = F0 .>= Q                   #logical where F0 >= Q  (decreasing function)
      yL[vvH] = y[vvH]                     # => root must be higher than y
      vvL     = F0 .< Q                    #logical where F0 < Q
      yH[vvL] = y[vvL]                     # -> root must be lower than y
    end

  end                                      #end different methods

  if n == 1                                 #1x1 array to scalar
    y = y[1]
  end

  return y

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function BondDurationMacaulay3Ps(Q,c0,t,Method=1,yLH=[-0.1;0.5],tol=1e-7,FaceValue=1)
#BondDurationMacaulay3Ps   Calculates bond durations. See also BondYieldToMat3Ps
#
#
#
#  Usage:     (D,Da,Dmac,ytm) = BondDurationMacaulay3Ps(Q,c0,t[,Method[,yLH[,tol[,FaceValue]]]])
#
#  Input:     Q         n vector, bond prices (for instance, 1.01)
#             c0        s1nn vector or mxn matrix, c rate, e.g. 0.06, or 0.09/2
#             t         mx1 vector, dates of coupon payments. Principal
#                       is paid at the same time as the last c. Dates
#                       should be expressed as fractions of the period
#                          Example: coupons after 2,8, and 14 months ->
#                       t = [0.333;0.333+1;0.333+2], to calculate semi-annual yield
#             Method    (optional) 1: Newton-Raphson; 2: bisection [1]
#             yLH       (optional) if Method==1: scalar, yLH(1) is initial guess of roots
#                                  if Method==2: 2x1 vector, yLH(1) is lower boundary and yLH(2) upper boundary
#             tol       (optional) scalar, convergence criterion for y, [1e-7]
#             FaceValue (optional), scalar or n vector, face value (default [1])
#
#  Output:    D        nx1 vector, (dollar) duration
#             Da        ""       , adjusted (or modified) duration
#             Dmac      ""       , Macaulays duration
#             ytm
#
#
#  Note:  (a)     The bond price satisfies (with yield expressed per period)
#
#                 Q =  c/(1+y)^t(1) + c/(1+y)^t(2) + ... + (1+c)/(1+y)^t(m)
#
#                 where y is the yield to maturity, c is per period coupon
#                 and Macaulay's duration is
#
#          D =  [t(1)*c/(1+y)^t(1) + t(2)*c/(1+y)^t(2) + ... + t(m)*(1+c)/(1+y)^t(m)]/Q
#
#
#
#
#  Reference: Campbell. Lo, MacKinlay, p 403.
#
#
#  Calls on:  BondYieldToMat3Ps
#
#  Paul Soderlind (Paul.Soderlind@unisg.ch), June 2009
#------------------------------------------------------------------------------

  Q = vecPs(Q)

  n = length(Q)
  m = length(t)

  c = deepcopy(c0)                          #breaks link with c0 argument
  if length(c) == 1                         #scalar c, single bond
    c = fill(c[1],(m,1))
  elseif size(c,1,2) == (1,n)
    c = repeat(c,outer=(m,1))
  end

  ytm = BondYieldToMat3Ps(Q,c,t,Method,yLH,tol,FaceValue)  #yield to maturity, effective interest rate
  y   = repeat(vecPs(ytm)',outer=(m,1))            #1xn -> mxn
  t   = repeat(vecPs(t),outer=(1,n))               #m -> mxn matrix

  if length(FaceValue) == 1
    c[end,:] = c[end:end,:] + FaceValue[1]    #add face value to last coupon payment
  else
    c[end,:] = c[end:end,:] + vecPs(FaceValue)'
  end

  cfac  = c.*t./((1+y).^(t+1))         #c/(1+y)^2 + 2c/(1+y)^3 + ...+ Tc/(1+y)^(T+1)

  D    = vec(sum(cfac,1))              #Duration
  Da   = D ./ Q                        #adjusted duration
  Dmac = D .* vec(1+y[1,:])./ Q        #Macaulays duration

  if n == 1                                 #1x1 array to scalar
    (D,Da,Dmac,ytm) = (D[1],Da[1],Dmac[1],ytm[1])
  end

  return D,Da,Dmac,ytm

end
#------------------------------------------------------------------------------
