#------------------------------------------------------------------------------
#  SRConf2Ex.jl
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

include("jlFiles/NWFn.jl")
include("jlFiles/SRFn.jl")
#------------------------------------------------------------------------------

xx   = readdlm("Data/FFmFactorsPs.csv",',',header=true)      
x    = xx[1]
Rme  = x[:,2]/100
T    = size(Rme,1)
#------------------------------------------------------------------------------

mu    = mean(Rme,1)                  #GMM to estimate mean and 2nd moment
mu_2  = mean(Rme.^2,1)
println("\n","[mu mu_2]",round([mu mu_2],4))

g    = [(Rme .- mu) (Rme.^2 .- mu_2)]    #moment conditions
T    = size(g,1)
gbar = mean(g,1)
println("\n","Sample moment conditions, gbar: ",round(gbar,4))
S = NWFn(g,1)                     #Var[sqrt(T)*gbar], Newey-West
D = -1
V = inv(D*inv(S)*D')              #Var[sqrt(T)*(mu,mu_2)]
println("Cov(params)")
println(round(V,6))
                                  #delta method
df = [ (mu_2/(mu_2 - mu.^2).^(3/2)) (-mu/(2*(mu_2 - mu.^2).^(3/2)))]
println("\n","Analytical derivatives",round(df,3))

par0 = [mu;mu_2]
SR0  = SRFn(par0)       #Sharpe ratio, function

Delta = 1e-7            #numerical derivatives
k     = 2               #number of parameters
Jac   = fill(NaN,(1,k))
for j = 1:k              #loop over columns (parameters)
  btilde    = par0 + 0.0 #break the link between par0 and btilde, now not identical
  btilde[j] = btilde[j] + Delta
  Jac[1,j]  = (SRFn(btilde)- SR0)/Delta
end
println("\n","Numerical derivatives, just for comparison",round(Jac,3))

Std_SR = sqrt(df*V*df'/T)
println("\n","SR and its Std",round([SR0 Std_SR],4))
println("SR and 90% conf band",round([SR0 (SR0-1.65*Std_SR) (SR0+1.65*Std_SR)],4))
#----------------------------------------------------------------------------
