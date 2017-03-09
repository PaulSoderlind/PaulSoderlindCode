#------------------------------------------------------------------------------
#  EstimationEx.jl
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
#
#  Paul.Soderlind@unisg.ch   May 2013, to Julia Jan 2016
#------------------------------------------------------------------------------

using Optim, PyPlot

if VERSION >= v"0.5.9"
  using SpecialFunctions
end

include("NormalHistLoss.jl")


CaseQ    = 3               #try 1,2, or 3
#------------------------------------------------------------------------------

                                                #DATA
CatBounds = [-2.00;-1.00; 0.00; 1.00; 2.00; 3.00; 4.00; 5.00; 6.00]   #catogories, -Inf<x<=Bound(1),Bond(1)<x<=Bound(2),..,Bound(n)<x<=Inf
CatMid    = [-2.50;-1.50;-0.50; 0.50; 1.50; 2.50; 3.50; 4.50; 5.50; 6.50]   #mid points, first and last interval are artificially closed

if CaseQ == 1          #only one interval has non-zero prob
  Prob = [ 0.00; 0.00; 0.00; 0.00; 1.00; 0.00; 0.00; 0.00; 0.00; 0.00]
elseif CaseQ == 2      #two intervals have non-zero prob
  Prob = [ 0.00; 0.00; 0.00; 0.00; 0.30; 0.70; 0.00; 0.00; 0.00; 0.00]
elseif CaseQ == 3      #3 or more intervals have non-zero prob
  Prob = [ 0.00; 0.00; 0.00; 0.10; 0.30; 0.40; 0.15; 0.05; 0.00; 0.00]   #probabilities
end
#------------------------------------------------------------------------------

if sum(Prob) != 1
  error("probabilities do not sum to one")
end

MeanCrude = sum(Prob .* CatMid)                    #crude mean and variance
VarCrude  = sum(Prob .* (CatMid-MeanCrude).^2)

BinWidth    = mean(diff(CatBounds))
VarSheppard = VarCrude - BinWidth^2/12       #variance, Sheppard's correction
#------------------------------------------------------------------------------

NActiveCat = sum(Prob .> 0)      #no. intervals with non-zero probabilities

if NActiveCat == 1        #if only one active intervals: triangular distribution
   parM = [MeanCrude sqrt(BinWidth^2/24) NActiveCat]
elseif NActiveCat == 2        #if only two active intervals: N(MeanCrude,VarSheppard)
   parM = [MeanCrude sqrt(VarSheppard) NActiveCat]
elseif NActiveCat > 2         #if three or more active intervals: N(estimate,estimate)
  par0 = [MeanCrude;sqrt(VarSheppard)]
  Sol = optimize(par->NormalHistLoss(par,Prob,CatBounds),par0)
  par1 = Optim.minimizer(Sol)
  parM = [par1' NActiveCat]
else
  parM = fill(NaN,(1,3))
end
println("\nmean, std, no. active intervals: ",round.(parM,3))
#------------------------------------------------------------------------------

function NormPdfPs(y,mu=0,s2=1)
#NormPdfPs    Returns pdf value of normally distributed (univariate) variables
  s = sqrt(s2)
  z    = (y .- mu)./s
  pdfy = exp.(-0.5*z.^2)./(sqrt(2*pi)*s)   #pdf of y ~ N(mu,s^2)
  return pdfy
end
#------------------------------------------------------------------------------

if parM[3] >= 3                           #fitted N(mu,s^2)
  y    = linspace(-2,6,101)
  pdfy = NormPdfPs(y,parM[1],parM[2]^2)
  figure()
    ha = bar(CatMid,Prob/BinWidth,color="yellow",align="center",width=BinWidth)
    xlim(-2,6)
    ylim(0,maximum(Prob/BinWidth)+0.03)
    title("Histogram and fitted \$N(\\mu,\\sigma^2)\$")
    hb = plot(y,pdfy)
end
#-------------------------------------------------------------------------------

