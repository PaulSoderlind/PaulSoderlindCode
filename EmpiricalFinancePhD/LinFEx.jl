#------------------------------------------------------------------------------
#  LinFEx.jl
#
#
#
#
#
#
#
#  Paul.Soderlind@unisg.ch   Oct 2015
#------------------------------------------------------------------------------

include("jlFiles/OlsFn.jl")
include("jlFiles/NWFn.jl")
#------------------------------------------------------------------------------

xx   = readdlm("Data/FFmFactorsPs.csv",',',header=true)
x    = xx[1]
Rme  = x[:,2]
RSMB = x[:,3]                #small minus big firms
RHML = x[:,4]                #high minus low book-to-market ratio
Rf   = x[:,5]                    #interest rate


x = readdlm("Data/FF25Ps.csv",',')  #no header line: x is matrix
R  = x[:,2:end]                  #returns for 25 FF portfolios
Re = R - repmat(Rf,1,size(R,2))  #excess returns for the 25 FF portfolios

(T,n) = size(Re)                 #no. obs and  no. test assets
#------------------------------------------------------------------------------

                                           #FF, testing alphas
cf        = [ones(T,1) Rme RSMB RHML]       #factors, constant and 3 FF factors
(b,epsM,) = OlsFn(Re,cf)
alfaM     = b[1:1,:]'                       #nx1

#g_    = HDirProdFn(cf,epsM)
g = fill(NaN,(T,size(cf,2)*n))            #moment conditions, regressors*residual
for t = 1:T
  g[t,:] = kron(cf[t:t,:],epsM[t:t,:])
end
S0   = NWFn(g,0)                          #ACov(sqrt(T)*gbar)
Sxx  = cf'cf/T
D0_1 = -kron(inv(Sxx),eye(n))
V    = D0_1*S0*D0_1'
WaldStat  = alfaM'inv(V[1:n,1:n]/T)*alfaM   #test if alpha=0

println("\nTesting alpha = 0, test statistic, df, and 10% critical value of Chi-square(25)")
println(round([WaldStat length(alfaM) 34.38],3))
#------------------------------------------------------------------------------

                                          #FF, testing ERe = lambda'beta
cf        = [ones(T,1) Rme RSMB RHML]     #factors, constant and 3 FF factors
K         = size(cf,2) - 1                #number of factors, excluding constant
(b,epsM,) = OlsFn(Re,cf)
bM        = b'                            #nx(1+K)

ERe    = mean(Re,1)'
betaM  = bM[:,2:end]
theta  = betaM'
lambda = (theta*betaM)\(theta*ERe)       #cross-sectional estimate of price of factor risk

#g_    = [HDirProdFn(cf,epsM),Re - repmat(lambda'*betaM',T,1)]
g = fill(NaN,(T,size(cf,2)*n))           #moment conditions, regressors*residual
for t = 1:T
  g[t,:] = kron(cf[t:t,:],epsM[t:t,:])
end
g = [g (Re - repmat(lambda'*betaM',T,1))]
p = length(bM) + length(lambda)             #no. parameters
q = size(g,2)                               #no. moment conditions

gbar = mean(g,1)'
S0   = NWFn(g,0)

Sxx  = cf'cf/T
D0LL = kron([0 lambda'],eye(n))             #lower left of D0
D0   = - [kron(Sxx,eye(n)) zeros(n*(1+K),K);
          D0LL              betaM           ]

A    = [eye(n*(1+K))     zeros(n*(1+K),n);     #A*gbar, traditional approach
        zeros(K,n*(1+K)) theta ]

Psia = eye(q) - D0*inv(A*D0)*A
Psi3 = Psia*S0*Psia'                       #Cov[sqrt(T)*gbar], rank q-p

WaldStat  = gbar'pinv(Psi3/T)*gbar            #test of moment conditions
println("\nTesting ERe = lambda'beta, test statistic, df and 10% critical value of Chi-square(22)")
println(round([WaldStat (q-p) 30.81],3))
#------------------------------------------------------------------------------
