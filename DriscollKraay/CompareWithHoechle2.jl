#------------------------------------------------------------------------------
#  CompareWithHoechle2.m
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

using JLD                            #install the JLD package by Pkg.add("JLD")
include("excisePs.jl")
include("HDirProdPs.jl")
include("HszDk5dwPs.jl")
#------------------------------------------------------------------------------

xTN = load("bidaskspread.jld","xTN")          #load data, from Hoechle

BA    = xTN[:,:,4]
TRMS  = xTN[:,:,5]
TRMS2 = xTN[:,:,6]
aVol  = xTN[:,:,7]
Size  = xTN[:,:,8]

(T,N) = size(BA)
#------------------------------------------------------------------------------

println("\nPooled OLS")
x   = [vec(aVol') vec(Size') vec(TRMS2') vec(TRMS') vec(ones(T,N)')]
y   = vec(BA')           #stack rows in column vector
yx, = excisePs([y x])     #prune NaNs
y   = yx[:,1]
x   = yx[:,2:end]

b    = x\y
res  = y - x*b
Covb = inv(x'x)*var(res,1)[1]
Stdb = sqrt(diag(Covb))
println("b, Std and tstat from OLS")
println(round([b Stdb b./Stdb],4))
#------------------------------------------------------------------------------

println("\nHszDk5dwPs, effectively weighting all periods equally, irrespective of number of obs")
x = ones(T,1)
z = cat(3,aVol,Size,TRMS2,TRMS,ones(T,N))
fnOutput = HszDk5dwPs(BA,x,z,false,8,1)
theta    = fnOutput[1]
stdDKj   = fnOutput[9]

println("theta, std and t-stat")
println(round([theta stdDKj theta./stdDKj],4))
#------------------------------------------------------------------------------

println("\nHszDk5dwPs, effectively weighting all obs equally")
x = ones(T,1)
z = cat(3,aVol,Size,TRMS2,TRMS,ones(T,N))
fnOutput = HszDk5dwPs(BA,x,z,false,8,0)
theta    = fnOutput[1]
stdDKj   = fnOutput[9]

println("theta, std and t-stat")
println(round([theta stdDKj theta./stdDKj],4))
#------------------------------------------------------------------------------
