#------------------------------------------------------------------------------
#  ExampleAPTsT.jl
#
#  For testing the AssetPricingTest3bPs (alpha test) and 
#                  AssetPricingTest3dPs (cross-sectional test) functions
#
#  The data is from the webpage of Kenneth R. French,
#  http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html
#
#
#
#  Paul.Soderlind@unisg.ch   May 2015, To Julia Oct/Dec 2015
#------------------------------------------------------------------------------

include("D:/PsCode/PsJulia/PaulS/excise.jl")
#include("OlsPs.jl")
include("D:/PsCode/PsJulia/PaulS/NewEst3Ps.jl")
include("D:/PsCode/PsJulia/PaulS/HDirProdPs.jl")
include("D:/PsCode/PsJulia/PaulS/AssetPricingTest3bPs.jl")
include("D:/PsCode/PsJulia/PaulS/AssetPricingTest3dPs.jl")
#------------------------------------------------------------------------------

vvRe = [1;6;13;19;25]             #which of the 25 FF portfolios to use as test assets

xx   = readdlm("Data/FFmFactorsPs.csv",',',header=true)      
x    = xx[1]
Rme  = x[:,2]
RSMB = x[:,3]                    #small minus big firms
RHML = x[:,4]                    #high minus low book-to-market ratio
Rf   = x[:,5]                    #interest rate

f    = Rme + 0.0
g    = [RSMB RHML]
h    = [g f]
    

x = readdlm("Data/FF25Ps.csv",',')  #no header line: x is matrix     
R  = x[:,2:end]                  #returns for 25 FF portfolios
Re = R .- Rf                      #excess returns for the 25 FF portfolios
Re = Re[:,vvRe]
(T,n) = size(Re)
K     = size(f,2)
L     = size(h,2)
M     = L - K
#------------------------------------------------------------------------------

println("\n------------AssetPricingTest3bPs-----------------")

ExtraTest  = Any[collect(1:n),collect(n+n*K+1:n+n*K+n)]
ExtraTestq = Any[zeros(n,1),zeros(n,1)]

fnOutput = AssetPricingTest3bPs(Re,Re,f,h,0,0,ExtraTest,ExtraTestq)
Joint_b = fnOutput[1]
ExtraWaldStat_b = fnOutput[8]

println("\n 3x1, WaldStat for alpha, delta, alpha-delta ")
println(round(Joint_b,2))
println("2x1, WaldStat for alpha, delta but via ExtraTest")
println(round(ExtraWaldStat_b,2))
#----------------------------------------------------------


println("\n ---------------AssetPricingTest3dPs: test of ERe - beta*lambda--------")

fnOutput = AssetPricingTest3dPs(Re,Re,f,h,0,0,ExtraTest,ExtraTestq,"GLS")
Joint_d = fnOutput[1]
ExtraWaldStat_d = fnOutput[8]
println("\n 3x1, WaldStat for pricing error 1, pricing error 2, difference of pricing errors ")
println(round(Joint_d,2))
println("2x1, WaldStat for lambda and psi")
println(round(ExtraWaldStat_d,2))


println("\n ---------------AssetPricingTest3dPs: test of ERe - beta*lambda, II--------")
println("set up to replicate TS approach")

vvR = Any[1:n,1:n]
q   = Any[zeros(n,1),zeros(n,1)]
BG  = Any[[zeros(K,n) eye(K)],[zeros(L,n) eye(L)]]
fnOutput = AssetPricingTest3dPs([Re f],[Re h],f,h,0,q,ExtraTest,ExtraTestq,BG,vvR)
Joint_d = fnOutput[1]
ExtraWaldStat_d = fnOutput[8]
println("3x1, WaldStat for pricing error 1, pricing error 2, difference of pricing errors ")
println(round(Joint_d,2))
println("2x1, WaldStat for lambda and psi")
println(round(ExtraWaldStat_d,2))
#------------------------------------------------------------------------------

