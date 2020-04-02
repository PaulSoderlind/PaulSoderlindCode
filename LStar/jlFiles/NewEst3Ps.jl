"""
    NewEst3Ps(g0,m=0,HHQ=false)

Calculate covariance matrix of sqrt(T)*sample average.


# Input
- `g0::Array`:  Txq array or T-element vector of q (1) moment conditions
- `m::Int`:     scalar, number of lags to use
- `HHQ::Bool`:  bool, true: do Hansen-Hodrick, else Newey-West

# Output
- `S::Array`: qxq covariance matrix of sum( g0/sqrt(T) )

 # Notice
- The CLT typically say that sqrt(T)*sample average ->d N(mu,Shat), where Shat is
what this code estimates. Clearly, Var(sample average) = Shat/T.
- Rule-of-thumb for m are floor(0.75*T^(1/3)) and floor(4*(T/100)^(2/9)).

# Reference
Newey and West, 1987 (Econometrica),  Newey, 1985 (Journal of Econometrics)

"""
function NewEst3Ps(g0,m=0,HHQ=false)

  (T,q) = (size(g0,1),size(g0,2))        #g is Txq
  m     = min(m,T-1)                     #number of lags

  g = g0 .- mean(g0,dims=1)              #normalizing to mean(g)=0

  S = g'g/T                              #(qxT)*(Txq)
  isa(S,Number) && (S=fill(S,1,1))       #to [1,1] matrix if S is a scalar

  for s = 1:m
    Omega_s = g[s+1:T,:]'g[1:T-s,:]/T      #same as Sum[g(t)*g(t-s)',t=s+1,T]
    HHQ ? w = 1 : w = (1 - s/(m+1))        #Hansen-Hodrick or Newey-West
    S       = S + w*(Omega_s + Omega_s')   #( ) is a matrix, S should be too
  end

  (q == 1) && (S = S[1,1])                 #to scalar

  return S

end
#-----------------------------------------------------------------------
