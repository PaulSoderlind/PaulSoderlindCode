function VAR1SimPs(A,epsilon,T,x0=0.0)
#VAR1SimPs   Calculates impulse response function of a VAR(1) system.
#
#             x(t) = A * x(t-1) +  epsilon(t), where x(t) is nx1
#
#  Usage:     xM = VAR1SimPs(A,epsilon,T,x0) or
#                = VAR1SimPs(A,epsilon,T)
#
#  Input:     A             nxn VAR(1) matrix, see above
#             epsilon       nx1 or 1xn vector of shocks in inital period, or Txn
#                           matrix with shocks in all periods
#             T             scalar, last period to calculate for
#             x0            nx1 or 1xn vector with starting values, optional
#
#  Output:    xM            T x n matrix, impulse response function
#
#
#
#  Paul.Soderlind@unisg.ch, to Julia Nov 2015
#-----------------------------------------------------------------------

  n = size(A,1)

  if isa(x0,Number)                               #if scalar
    x0 = fill(x0,n)
  end

  if isa(epsilon,Number)
    epsilon = fill(epsilon,n)
  end

  if length(epsilon) == n
    epsilon = vcat(vec(epsilon)',zeros(T-1,n))
  end

  x1_t_1 = vec(x0)                                #starting vector
  xM     = fill(NaN,(T,n))                        #to put results in
  for t = 1:T                                     #loop over time periods
    x1      = A*x1_t_1 + epsilon[t,:]             #[t,:] gives column vec in 0.5+
    xM[t,:] = x1
    x1_t_1  = copy(x1)
  end

  return xM

end
#-----------------------------------------------------------------------
