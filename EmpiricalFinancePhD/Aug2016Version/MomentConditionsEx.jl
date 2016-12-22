#------------------------------------------------------------------------------
#  MomentConditionsEx.jl
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
#------------------------------------------------------------------------------

xx = readdlm("Data/FFmFactorsPs.csv",',',header=true)   #start on line 2, column 1
x  = xx[1]
ym = x[:,1]         #[yearmonth]
x  = x[:,2]/100     #excess market returns
T  = size(x,1)
#------------------------------------------------------------------------------

muHat = mean(x,1)                     #it is clear that this is the GMM estimate

g    = x .- muHat                     #moment condition, .- 'broadcasts muHat' to x
gbar = mean(g,1)
println("Sample moment condition at estimate: ",round(gbar,4))

T   = size(g,1)
S   = var(x)
D   = -1
V   = inv(D'inv(S)*D)

println("\n[muHat Std(muHat)]")
println(round([muHat sqrt(V/T)],4))
#----------------------------------------------------------------------------

println("\nAs an alternative, use NW for S")

Sb = NWFn(g,1)            #Newey-West covariance matrix

Vb = inv(D'inv(Sb)*D)

println("\n[muHat Std(muHat)] according to NW")
println(round([muHat sqrt(Vb/T)],4))
#----------------------------------------------------------------------------
