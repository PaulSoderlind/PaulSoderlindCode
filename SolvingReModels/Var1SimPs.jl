function Var1SimPs(A,epsilon,T,x0=0.0)
#Var1SimPs   Calculates impulse response function of a VAR(1) system.
#
#             x(t) = A * x(t-1) +  epsilon(t), where x(t) is nx1
#
#  Usage:     xM = Var1SimPs(A,epsilon,T,x0) or
#                = Var1SimPs(A,epsilon,T)
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
#  Paul.Soderlind@unisg.ch, May 2001, to Julia Nov 2015
#-----------------------------------------------------------------------

  n = size(A,1)

  if length(epsilon) == n
    epsilon2       = zeros(T,n)                   #non-zero only in first period
    epsilon2[1,:]  = vec(epsilon)'                #to accomodate T=1
    epsilon        = epsilon2
  end
  if isa(x0,Number)                               #if scalar
    x0 = repmat([x0],n,1)
  end

  x1_t_1 = vec(collect(x0))                       #starting vector
  xM     = fill(NaN,(T,n))                        #to put results in
  for t = 1:T                                     #loop over time periods
    x1      = A*x1_t_1 + epsilon[t:t,:]'          #[t] to keep as row vector
    xM[t,:] = x1'
    x1_t_1  = deepcopy(x1)
  end

  return xM

end
#-----------------------------------------------------------------------
