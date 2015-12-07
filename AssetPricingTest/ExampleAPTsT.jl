#------------------------------------------------------------------------------
#  ExampleAPTsT.jl
#
#  For testing the AssetPricingTest3bPs function
#
#  The data is from the webpage of Kenneth R. French,
#  http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html
#
#
#
#  Paul.Soderlind@unisg.ch   May 2015, To Julia Oct 2015
#------------------------------------------------------------------------------


include("excise.jl")
include("NewEst3Ps.jl")
include("HDirProdPs.jl")
include("AssetPricingTest3bPs.jl")
#------------------------------------------------------------------------------

vvRe = [1;6;13;19;25]             #which of the 25 FF portfolios to use as test assets

xx   = readdlm("FFmFactorsPs.csv",',',header=true)      
x    = xx[1]
Rme  = x[:,2]
RSMB = x[:,3]                    #small minus big firms
RHML = x[:,4]                    #high minus low book-to-market ratio
Rf   = x[:,5]                    #interest rate

f    = Rme + 0.0
g    = [RSMB RHML]
h    = [g f]
    

x = readdlm("FF25Ps.csv",',')  #no header line: x is matrix     
R  = x[:,2:end]                  #returns for 25 FF portfolios
Re = R .- Rf                      #excess returns for the 25 FF portfolios
Re = Re[:,vvRe]
(T,n) = size(Re)
K     = size(f,2)
L     = size(h,2)
M     = L - K
#------------------------------------------------------------------------------

println("------------AssetPricingTest3bPs-----------------")

ExtraTest  = Any[collect(1:n),collect(n+n*K+1:n+n*K+n)]
ExtraTestq = Any[zeros(n,1),zeros(n,1)]

fnOutput = AssetPricingTest3bPs(Re,Re,f,h,0,0,ExtraTest,ExtraTestq)
Joint_b = fnOutput[1]
ExtraWaldStat_b = fnOutput[8]

println("AssetPricingTest3bPs")
println("3x1, WaldStat for alpha, delta, alpha-delta ")
println(round(Joint_b,2))
println("2x1, WaldStat for alpha, delta but via ExtraTest")
println(round(ExtraWaldStat_b,2))
#----------------------------------------------------------
