#------------------------------------------------------------------------------
#  GmmOptEx.jl
#
#
#
#
#
#
#
#
#
#
#
#
#  Paul.Soderlind@unisg.ch   Oct 2015
#------------------------------------------------------------------------------


include("jlFiles/Gmm4MomFn.jl")
include("jlFiles/NWFn.jl")

using Optim
#------------------------------------------------------------------------------

xx   = readdlm("Data/FFmFactorsPs.csv",',',header=true)      
x    = xx[1]
ym   = x[:,1]              #[yearmonth]
x    = x[:,2]              #excess market returns
#------------------------------------------------------------------------------

T  = size(x,1)

mu = mean(x)                      #same as setting A*gbar=0
s2 = var(x)*(T-1)/T               #var() uses 1/(T-1) formula 

par_a = [mu;s2]
k     = length(par_a)
println("\n","parameters and traditional std(parameters)")
println(round([par_a [sqrt((s2/T));sqrt(2*s2^2/T)]],4))


(g,gbar) = Gmm4MomFn(par_a,x)        #Tx4, moment conditions
q = size(g,2)
A = [1 0 0 0;                       #A in A*gbar=0 (here: all weight on first two moments)
     0 1 0 0]
println("\n","Checking if mean of A*g_t = 0")
println(round(A*gbar,4))
D  = [-1                  0;                #Jacobian
      -2*mean(x-mu)      -1;
      -3*mean((x-mu).^2)   0;
      -4*mean((x-mu).^3)  -6*s2]
S  = NWFn(g,1)
V3 = inv(A*D)*A*S*A'inv(A*D)'
println("\n","parameter, std(parameters)")
println(round([par_a sqrt(diag(V3/T))],4))
#------------------------------------------------------------------------------
                                          #gbar'W*gbar
W     = diagm([1;1;0;0])                  #weighting matrix
x1    = optimize(par->Gmm4MomLossFn(par,x,W),par_a)
par_b = x1.minimum
g,    = Gmm4MomFn(par_b,x)              #Tx4, moment conditions, evaluated at point estimates
S     = NWFn(g,1)                         #variance of sqrt(T)"gbar, NW with 1 lag
V2    = inv(D'W*D)*D'W*S*W'D*inv(D'W*D)
println("\n","parameter, std(parameters)")
println(round([par_b sqrt(diag(V2/T))],4))
#------------------------------------------------------------------------------

#W,S, par_c should be initialized outside loop to make them visible after it         
par_c = par_b + 0.0
Dpar  = 1
i     = 1
println("\n","iterating over W")
while (Dpar > 1e-3) || (i < 2)    #require at least one iteration
  println("iteration  ", i, ": old and new parameters")
  par_b           = par_c + 0.0     #important, par_b=par_c would make them always identical
  W               = inv(S)                           
  x1              = optimize(par->Gmm4MomLossFn(par,x,W),par_b)   #use last estimates as starting point
  par_c           = x1.minimum
  g,              = Gmm4MomFn(par_c,x)
  S               = NWFn(g,1)
  Dpar            = maximum(abs(par_c-par_b))
  tst1 = (Dpar > 1e-5)
  println(round([par_b par_c],4)," ",round(Dpar,4)," ",tst1)
  i = i + 1
end

V2 = inv(D'W*D)*D'W*S*W'D*inv(D'W*D)
V1 = inv(D'inv(S)*D)
println("\n","parameter, std_version2(parameters), std_version1(parameters)")
println(round([par_c sqrt(diag(V2/T)) sqrt(diag(V1/T))],4))
#----------------------------------------------------------------------------
