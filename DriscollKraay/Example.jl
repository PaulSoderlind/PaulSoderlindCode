#------------------------------------------------------------------------------
#  Example.jl
#
#
#  This file imports the data from a jld file with Hoechle's data and then
#  runs the regressions.
#
#  Compare the results to Table 2 in Hoechle, 2007, "Robust Standard Errors...",
#  Stata Journal, 7, 281-312.
#
#
#
#
#
#
#
#  Paul.Soderlind@unisg.ch   Apr 2016
#------------------------------------------------------------------------------

using Printf, DelimitedFiles, LinearAlgebra, Statistics

include("FindNoNaNPs.jl")
include("HDirProdPs.jl")
include("iterationPrintPs.jl")
include("HszDk5dwPs.jl")
#------------------------------------------------------------------------------

#using JLD
#xTN = load("bidaskspread.jld","xTN")          #load data, from Hoechle
#(BA,TRMS,TRMS2,aVol,Size)    = [xTN[:,:,i] for i=4:8]

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

vv  = FindNoNaNPs(1,y,x)
(y,x) = (y[vv],x[vv,:])    #prune NaNs

b    = x\y
res  = y - x*b
Covb = inv(x'x)*var(res)
Stdb = sqrt.(diag(Covb))
println("b, Std and tstat from OLS")
display(round.([b Stdb b./Stdb],digits=4))
#------------------------------------------------------------------------------

println("\nHszDk5dwPs, effectively weighting all obs equally")
x = ones(T)
z = cat(aVol,Size,TRMS2,TRMS,ones(T,N),dims=3)
fnO = HszDk5dwPs(BA,x,z,false,8,0)
theta    = fnO.theta
stdDKj   = fnO.stdDKj

println("theta, std and t-stat")
display(round.([theta stdDKj theta./stdDKj],digits=4))

println("Compare the results to Table 2 in Hoechle, 2007, 'Robust Standard Errors...'")
#------------------------------------------------------------------------------

println("\nHszDk5dwPs, effectively weighting all periods equally, irrespective of number of obs")
x = ones(T,1)
z = cat(aVol,Size,TRMS2,TRMS,ones(T,N),dims=3)
fnO = HszDk5dwPs(BA,x,z,false,8,true)
theta    = fnO.theta
stdDKj   = fnO.stdDKj

println("theta, std and t-stat")
display(round.([theta stdDKj theta./stdDKj],digits=4))
#------------------------------------------------------------------------------
