"""
    BondPricePs(y,c0,t,FaceValue=1)

Calculates coupon bond price from the spot rate curve and the payment stream
(arbitrary periods and coupons), m periods


# Input
- ```y::Vector or Number```:    m vector or scalar, effective spot rates, e.g. [0.07;0.08], m periods
- ```c0::Vector or Number```:   m vector or scalar, coupon rate, e.g. 0.06 or 0.09/2, m periods
- ```t::Vector```:              m vector, time to the cash flows
- ```FaceValue::Number```:   (optional) face value, scalar [1]

# Output
- ```Q::Number```:  scalar, bond price (eg. 1.01)

# Notice
The face value is paid at the same time as the last coupon.
Dates (t) should be expressed as fractions of the period.

The following is calculated:

Q =  c[1]/(1+y(1))^t(1) + c[2]/(1+y(2))^t(2) + ... + (1+c[m])/(1+y(m))^t(m)

Notice that the coupons and interest rates depend on time to payment (optionally)


# Example
y = [0.0626;0.07],
c = [0.09; 0.09],
t = 1:2 gives Q = 1.0367


Paul Soderlind (Paul.Soderlind@unisg.ch)

"""
function BondPricePs(y,c0,t,FaceValue=1)

  isa(c0,Number) ? (c = fill(c0,length(t))) : (c = copy(c0))
  c[end] = c[end] + FaceValue

  cdisc  = c./((1.0.+y).^t)        #c/(1+y(1))^t1 + c/(1+y(2))^t2 + ...+ c/(1+y(m))^tm
  Q      = sum(cdisc)              #price (sum(.,1)) to handle also m=1

  return Q

end
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
"""
    BondPriceQDyPs(y,c0,t,FaceValue=1)

Calculates coupon bond price and derivative.

# Input
- see BondPricePs

# Output
- ```Q::Number```:  scalar, bond price (eg. 1.01)
- ```dQ_dy::Number```:  scalar, derivative of bond price wrt y (valid only when y is a scalar)

"""
function BondPriceQDyPs(y,c0,t,FaceValue=1)

  isa(c0,Number) ? (c = fill(c0,length(t))) : (c = copy(c0))
  c[end] = c[end] + FaceValue

  cdisc  = c./((1.0.+y).^t)        #c/(1+y(1))^t1 + c/(1+y(2))^t2 + ...+ c/(1+y(m))^tm
  Q      = sum(cdisc)              #price (sum(.,1)) to handle also m=1

  dcdisc = cdisc .* (-t./(1.0.+y))        #derivative wrt ytm
  dQ_dy  = sum(dcdisc)                    #only strictly valid when y is a scalar

  return Q, dQ_dy

end
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
"""
    BondYieldToMatPs(Q,c0,t,FaceValue=1,method=1,yLH=[-0.1;0.5],tol=1e-7)

Calculate yield to maturity from bond price (several methods available).
Can also be used for general IRR calculations. Works with arbitrary coupon periods.

# Input
- ```Q::Number```:          scalar
- ```c::Number```:          see BondPricePs
- ```t::Number```:          BondPricePs
- ```FaceValue::Number```:  optional, scalar or n vector, face value, default [1]
- ```method::Number```:     optional, 1: Newton-Raphson; 2: bisection [1]
- ```yLH::Number```:        optional, if method==1: scalar, yLH[1] is initial guess of roots
                            if method==2: 2 vector, yLH[1] is lower boundary and yLH[2] upper boundary
- ```tol::Number```:        optional, scalar, convergence criterion for y, [1e-7]

# Output
- ```y::Number```           effective yield to maturity (per period)

# Notice
- method 1 (Newton-Raphson) is pretty fast.
- method 2 (bisection) is a bit slow, but very robust.
- To use for IRR calculations, set Q=0 and FaceValue=0
- For Newton-Raphson, we could use (c+(FaceValue-B)/n)/((FaceValue+B)/2) as a good starting value


# Requires
- BondPricePs

"""
function BondYieldToMatPs(Q,c0,t,FaceValue=1,method=1,yLH=[-0.1;0.5],tol=1e-7)

  if method == 1                   #Newton-Raphson

    y = yLH[1]
    Dy = 1e+198
    while abs(Dy) > tol
      (F0,J0) = BondPriceQDyPs(y,c0,t,FaceValue)    #F(y0), derivative
      Dy = (Q - F0)/J0
      y  = y + Dy
    end

  elseif method == 2                   #bisection

    (yL,yH) = (yLH[1],yLH[2])                #lower and upper boundary for yield
    if yL > yH
      println("warning: lower bound greater than upper bound")
    end
    F0 = BondPricePs(yL,c0,t,FaceValue)     #create starting value, so [yL,yH] brackets the roots
    while F0 < Q                             #as long as theoretical < actual price
      yL = yL .- 0.01
      F0 = BondPricePs(yL,c0,t,FaceValue)
    end
    F0 = BondPricePs(yH,c0,t,FaceValue)
    while F0 > Q
      yH = yH .+ 0.01
      F0 = BondPricePs(yH,c0,t,FaceValue)
    end

    y = (yL + yH)/2
    while (yH-yL) > tol                    #iteration loop
      y = (yL + yH)/2                      #mid point for yield
      F0 = BondPricePs(y,c0,t,FaceValue)  #price at guessed yield
      if F0 >= Q                           #logical where F0 >= Q  (decreasing function)
        yL = y                             # => root must be higher than y
      else
        yH = y                             # -> root must be lower than y
      end
    end

  end                                      #end different methods

  return y

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
"""
    BondDurationPs(Q,c0,t,FaceValue=1,method=1,yLH=[-0.1;0.5],tol=1e-7)

Calculate bond durations. See also BondYieldToMatPs


# Input
- see BondYieldToMatPs

# Output
- ```D::Number```:        scalar, (dollar) duration
- ```Da::Number```:       scalar, adjusted (or modified) duration
- ```Dma::Numberc```:     scalar, Macaulays duration
- ```ytm::Number```:      scalar


# Notice
The bond price satisfies (with yield expressed per period)

Q =  c/(1+y)^t(1) + c/(1+y)^t(2) + ... + (1+c)/(1+y)^t(m)

where y is the yield to maturity, c is per period coupon
and Macaulay's duration is

D =  [t(1)*c/(1+y)^t(1) + t(2)*c/(1+y)^t(2) + ... + t(m)*(1+c)/(1+y)^t(m)]/Q


# Requires
- BondPricePs, BondYieldToMatPs

"""
function BondDurationPs(Q,c0,t,FaceValue=1,method=1,yLH=[-0.1;0.5],tol=1e-7)

  y = BondYieldToMatPs(Q,c0,t,FaceValue,method,yLH,tol)  #ytm

  (_,Dy) = BondPriceQDyPs(y,c0,t,FaceValue)
  D      = - Dy                                   #duration
  Da     = D / Q                                  #adjusted duration
  Dmac   = D*(1+y)/ Q                             #Macaulays duration

  return D, Da, Dmac, y

end
#------------------------------------------------------------------------------
