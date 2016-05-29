#------------------------------------------------------------------------------
#  LStarTest.jl
#
#  Test of LSTAR code used in "Carry Trade" paper.
#
#
#
#
#
#
#
#
#  Paul.Soderlind@unisg.ch   18 January 2013, to Julia Nov 2015
#------------------------------------------------------------------------------


include("readdlmFixPs.jl")            #CHANGE THESE PATHS to where you put the files
include("excise.jl")
include("OlsPs.jl")
include("NewEst3Ps.jl")
include("NumJac3Ps.jl")
include("OlsLStar3Ps.jl")

using Optim
#------------------------------------------------------------------------------


x     = readdlm("CTData.csv",',')      #date, FXV, Return_CT, SP, Ty
x     = readdlmFixPs(x,0)              #missing -> 0 (as in Matlab example)
FXV   = x[:,2]                         #FX volatility
Re_CT = x[:,3]                         #carry trade return
SP    = x[:,4]                         #S&P return
Ty    = x[:,5]                         #treasury return

T = size(Re_CT,1)

gM = collect(linspace(1,4,15))
cM = collect(linspace(-1,1,17))

z = (FXV - mean(FXV))/std(FXV)          #standardised regime variable
#------------------------------------------------------------------------------

gKeep = [NaN NaN]                       #set to NaN if estimated in NLS, otherwise imposed
fnOutput = OlsLStar3Ps(Re_CT,[ones(T,1) SP Ty],Array{Float64}(T,0),true,z,gM,cM,gKeep)
#Any[sse,theta,Stdtheta,Covtheta,b,Stdb_ols,R2a,Gquant,gc,sseM,yHat,slopeDiff,yHatLH]

theta     = fnOutput[2]
Stdtheta  = fnOutput[3]
Covtheta  = fnOutput[4]
b         = fnOutput[5]
Stdb_ols  = fnOutput[6]
R2a       = fnOutput[7]
slopeDiff = fnOutput[12]



println("\n","theta is [g;c;b_low;b_high;slopes without regimes]")
println("[theta Stdtheta]")
println(round([theta Stdtheta],4))

println("\n","difference of slope (high minus low state), t-stat")
println(round(slopeDiff,4))
println("-----------------------------------------------------")
#error()
