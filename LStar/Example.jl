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


using Statistics, DelimitedFiles, LinearAlgebra, Optim, ForwardDiff

include("FindNoNaNPs.jl")
include("NewEst3Ps.jl")
include("OlsPs.jl")
include("OlsLStar3Ps.jl")
#------------------------------------------------------------------------------

x     = readdlm("CTData.csv",',')      #date, FXV, Return_CT, SP, Ty

FXV   = x[:,2]                         #FX volatility
Re_CT = x[:,3]                         #carry trade return
SP    = x[:,4]                         #S&P return
Ty    = x[:,5]                         #treasury return

T = size(Re_CT,1)

gM = range(1,stop=4,length=15)
cM = range(-2,stop=3,length=21)

z = (FXV .- mean(FXV))/std(FXV)        #standardised regime variable
#------------------------------------------------------------------------------

gcKeep = [NaN;NaN]                       #(gamma,c) set to NaN if estimated in NLS, otherwise imposed
#gcKeep = [3.0;NaN]                      #try to uncomment this to restrict gamma to 3
#gcKeep = [NaN;1.3]                      #or this


w = zeros(T,0)
(theta,Stdtheta,fnOutput) = OlsLStar3Ps(Re_CT,[ones(T) SP Ty],w,true,z,gM,cM,gcKeep)

if all(isnan.(gcKeep))
  println("\n","theta is [g;c;b_low;b_high;(slopes without regimes, if any)]")
elseif isnan.(gcKeep) == [true;false]
  println("\n","theta is [g;b_low;b_high;(slopes without regimes, if any)]")
elseif isnan.(gcKeep) == [false;true]
  println("\n","theta is [c;b_low;b_high;(slopes without regimes, if any)]")
end
println("[theta Stdtheta]")
display(round.([theta Stdtheta],digits=3))

println("\n","difference of slope (b, high minus low state), t-stat")
display(round.(fnOutput.slopeDiff,digits=3))
println("-----------------------------------------------------")

#------------------------------------------------------------------------------

#Plot of the loss function. Comment out this if you do not have PyPlot installed.

using PyPlot
close("all")

PyPlot.matplotlib[:rc]("mathtext",fontset="stix")
PyPlot.matplotlib[:rc]("font",family="STIXGeneral",style="normal",size=18)
PyPlot.PyObject(PyPlot.axes3D)          #to activate mplot3
fig = figure(figsize=(16,16/1.2))
ax = gca(projection="3d")
  surf(gM,cM,copy(fnOutput.sseM'),rstride=2,cstride=2,cmap=ColorMap("summer"),alpha=0.8)
  scatter(fnOutput.gcHat[1],fnOutput.gcHat[2],zs=minimum(fnOutput.sseM),s=200,color="k")
  ax[:view_init](elev=20.0, azim=10)
  title("Sum of squared residuals (optimum at dot)")
  xlabel(L"$\gamma$")
  ylabel(L"$c$")
#------------------------------------------------------------------------------

