function NewEst3Ps(g0,m)
#NewEst3Ps    Calculates covariance matrix of sqrt(T)*sample average.
#
#
#  Usage:     Shat = NewEst3Ps(hhat,m)
#
#  Input:     g0          TxK matrix of moment functions
#             m           order of autoregression
#
#  Output:    Shat        covariance matrix of sum( hhat/sqrt(T) ),
#                         that is the covariance matrix of sqrt(T)*mean(hhat)
#
#
#  Note: The CLT typically say that sqrt(T)*sample average ->d N(mu,Shat),
#        where Shat is what this code estimates. Clearly, Var(sample average) =
#        Shat/T.
#
#
#  Reference: Newey and West, 1987, Econometrica
#             Newey, 1985, Journal of Econometrics
#
#
#  Paul Soderlind (Paul.Soderlind@unisg.ch), to Julia Oct 2015
#-----------------------------------------------------------------------

  T = size(g0,1)                     #g is Txq
  m = min(m,T-1)                     #number of lags

  g = g0 .- mean(g0,1)               #Normalizing to Eg=0

  S = g'g/T                          #(qxT)*(Txq)
  for s = 1:m
    Omega_s = g[s+1:T,:]'g[1:T-s,:]/T   #same as Sum[g(t)*g(t-s)',t=s+1,T]
    S       = S + (1 - s/(m+1))*(Omega_s + Omega_s')
  end

  return S

end
#-----------------------------------------------------------------------
