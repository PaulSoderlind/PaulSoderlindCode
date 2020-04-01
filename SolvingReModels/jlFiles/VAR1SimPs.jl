"""
VAR1SimPs(A,epsilon,T,x0=0.0)

Calculate impulse response function of a VAR(1) system
x(t) = A * x(t-1) +  epsilon(t), where x(t) is nx1

# Input
- `A::Matrix`:                  nxn VAR(1) matrix, see above
- `epsilon::Number or Vector`:  n-vector of shocks in inital period, or Txn matrix with shocks in all periods
- `T::Number`:                  scalar, last period to calculate for
- `x0::Number or Vector`:       n-vector with starting values, optional

# Output
- `xM::Matrix`:               Txn matrix, impulse response function

Paul.Soderlind@unisg.ch, to Julia Nov 2015

"""
function VAR1SimPs(A,epsilon,T,x0=0.0)

  n = size(A,1)

  isa(x0,Number)         && (x0      = fill(x0,n))           #if scalar
  isa(epsilon,Number)    && (epsilon = fill(epsilon,n))
  (length(epsilon) == n) && (epsilon = vcat(vec(epsilon)',zeros(T-1,n)))

  xM      = fill(NaN,(T,n))                        #to put results in
  xM[1,:] = A*vec(x0) + epsilon[1,:]
  for t = 2:T                                      #loop over time periods
    xM[t,:] = A*xM[t-1,:] + epsilon[t,:]
  end

  return xM

end
#-----------------------------------------------------------------------
