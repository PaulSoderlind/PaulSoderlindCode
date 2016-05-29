#------------------------------------------------------------------------------
#  GarchEx.jl
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
#  Paul.Soderlind@unisg.ch   Oct 2015
#------------------------------------------------------------------------------

include("jlFiles/garch11LL.jl")

using Optim, PyPlot
#------------------------------------------------------------------------------

xx   = readdlm("Data/FFdSizePs.csv",',',header=true)
x    = xx[1]
ymd  = x[:,1]     #[YearMonthDay]
y    = x[:,2]     #returns for the size portfolio we want to study

yx = [y[2:end] y[1:end-1] ones(size(y,1)-1,1)]     #y(t),y(t-1),1
y  = yx[:,1]
x  = yx[:,2:3]
#------------------------------------------------------------------------------

#mean equation, y = x'*b
#GARCH(1,1) equation: s2(t) = alfa0 + alfa1*u(t-1)^2 + beta1*s2(t-1)


par0 = [0;mean(y,1);(var(y,1)/5);0.1;0.6]
(loglik,s2,yhat) = garch11LL(par0,yx)            #just testing the log lik

x1 = optimize(par->garch11LLLoss(par,yx),par0)   #do MLE by optimization with optimize, minimize -sum(LL)
parHata = x1.minimum
parHata[end-2:end] = abs(parHata[end-2:end])     #since the likelihood function uses abs(these values)
println("\nParameter estimates")
println(round(parHata,4))

(loglik,s2,ER) = garch11LL(parHata,yx)
VaR95          = -(ER - 1.64*sqrt(s2))

figure()
  plot(VaR95)
  title("VaR (95%)")

CovRatio = mean((-y) .>= VaR95)                          #coverage ratio for VaR
println("\nCoverage ratio for VaR(95%): ",round(CovRatio,3))
#------------------------------------------------------------------------------
