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

include("excise.jl")
include("NewEst3Ps.jl")
include("HDirProdPs.jl")
include("AssetPricingTest3bPs.jl")
include("AssetPricingTest3dPs.jl")
#------------------------------------------------------------------------------

vvRe = [1;6;13;19;25]             #which of the 25 FF portfolios to use as test assets

xx   = readdlm("FFmFactorsPs.csv",',',header=true)
x    = xx[1]
Rme  = x[:,2]
RSMB = x[:,3]                    #small minus big firms
RHML = x[:,4]                    #high minus low book-to-market ratio
Rf   = x[:,5]                    #interest rate
xx = nothing
x  = nothing


f    = Rme + 0.0
g    = [RSMB RHML]
h    = [f g]


x = readdlm("FF25Ps.csv",',')       #no header line: x is matrix
R  = x[:,2:end]                     #returns for 25 FF portfolios
Re = R .- Rf                        #excess returns for the 25 FF portfolios
Re = Re[:,vvRe]
(T,n) = size(Re)
K     = size(f,2)
L     = size(h,2)
M     = L - K
x  = nothing
#------------------------------------------------------------------------------

println("\n------------AssetPricingTest3bPs-----------------")

fnOutput = AssetPricingTest3bPs(Re,Re,f,h)
Joint_b = fnOutput[1]
Indiv_b = fnOutput[2]

println("\n 3x1, WaldStat for alpha, delta, alpha-delta ")
println(round(Joint_b,2))
println("\n nx3, Individual t-stats for alpha, delta, alpha-delta ")
display(round(Indiv_b',2))
#----------------------------------------------------------


println("\n ---------------AssetPricingTest3dPs: test of ERe - beta*lambda--------")

fnOutput = AssetPricingTest3dPs(Re,Re,f,h,0,0,"GLS")
Joint_d = fnOutput[1]
println("\n 3x1, WaldStat for pricing error 1, pricing error 2, difference of pricing errors ")
println(round(Joint_d,2))


println("\n ---------------AssetPricingTest3dPs: test of ERe - beta*lambda, II--------")
println("set up to replicate TS approach")

vvR = Any[1:n,1:n]
q   = Any[zeros(n,1),zeros(n,1)]
BG  = Any[[zeros(K,n) eye(K)],[zeros(L,n) eye(L)]]
fnOutput = AssetPricingTest3dPs([Re f],[Re h],f,h,0,q,BG,vvR)
Joint_d = fnOutput[1]
println("3x1, WaldStat for pricing error 1, pricing error 2, difference of pricing errors ")
println(round(Joint_d,2))
#------------------------------------------------------------------------------

