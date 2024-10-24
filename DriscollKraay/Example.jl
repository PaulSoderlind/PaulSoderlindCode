#------------------------------------------------------------------------------
#  Example.jl
#
#
#  This file imports the data from csv files with Hoechle's data and then
#  runs the regressions.
#
#  Compare the results to Table 2 in Hoechle, 2007, "Robust Standard Errors...",
#  Stata Journal, 7, 281-312.
#
#
#
#
#  Paul.Soderlind@unisg.ch
#------------------------------------------------------------------------------

using Printf, DelimitedFiles, LinearAlgebra, Statistics

include("jlFiles/printmat.jl")
include("jlFiles/FindNN.jl")
include("jlFiles/PanelOls.jl")
#------------------------------------------------------------------------------

BA    = readdlm("Data/BA.csv",',')
TRMS  = readdlm("Data/TRMS.csv",',')
TRMS2 = readdlm("Data/TRMS2.csv",',')
aVol  = readdlm("Data/aVol.csv",',')
Size  = readdlm("Data/Size.csv",',')

(T,N) = size(BA)
#------------------------------------------------------------------------------

println("\nPooled OLS")
x   = [vec(aVol') vec(Size') vec(TRMS2') vec(TRMS') vec(ones(T,N)')]
y   = vec(BA')           #stack rows in column vector




vv  = FindNN(y,x)
(y,x) = (y[vv],x[vv,:])    #prune NaNs

b    = x\y
res  = y - x*b
Covb = inv(x'x)*var(res)
Stdb = sqrt.(diag(Covb))

println("OLS, assuming iid errors")
colNames = ["coef","std","t-stat"]
rowNames = ["aVol","Size","TRMS2","TRMS","c"]
printmat(b,Stdb,b./Stdb;colNames,rowNames,prec=4)
#------------------------------------------------------------------------------

x0 = cat(aVol,Size,TRMS2,TRMS,ones(T,N),dims=3)
x  = permutedims(x0,[1,3,2])               #TxNxK -> TxKxN

fnO = PanelOls(BA,x,0;FixNaNQ=true)       #no lags in covariance estimation
θ     = fnO.theta                           #coefs
stdLS = sqrt.(diag(fnO.CovLS))          #standard errors from LS, White, DK
stdW  = sqrt.(diag(fnO.CovW))

fn8 = PanelOls(BA,x,8;FixNaNQ=true)       #8 lags in covariance estimation
stdNW = sqrt.(diag(fn8.CovW))
stdDK = sqrt.(diag(fn8.CovDK))


println("panel regression, LS standard errors")
printmat(θ,stdLS,θ./stdLS;colNames,rowNames,prec=4)

println("panel regression, White's standard errors")
printmat(θ,stdW,θ./stdW;colNames,rowNames,prec=4)

println("panel regression, Newey-West's standard errors")
printmat(θ,stdNW,θ./stdNW;colNames,rowNames,prec=4)

println("panel regression, DK's standard errors")
printmat(θ,stdDK,θ./stdDK;colNames,rowNames,prec=4)

println("Compare the results to Table 2 in Hoechle, 2007, 'Robust Standard Errors...'")
#------------------------------------------------------------------------------
