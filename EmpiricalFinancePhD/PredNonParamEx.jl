#------------------------------------------------------------------------------
#  PredNonParamEx.jl
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
#  Paul.Soderlind@unisg.ch   22 May 2013
#------------------------------------------------------------------------------

include("jlFiles/excise.jl")
include("jlFiles/lagnPs.jl")
include("jlFiles/KernRegNormalFn.jl")
include("jlFiles/KernRegUniformFn.jl")
include("jlFiles/OlsFn.jl")


CrossValCalcIt = 1       #1: do cross validation calculations
#------------------------------------------------------------------------------

xx  = readdlm("Data/FFdSizePs.csv",',',header=true)      
R   = xx[1]
ymd = R[:,1]             #[YearMonthDay]
R   = R[:,11]            #returns for the size portfolio we want to study

xx  = excise([R lagnPs(R)])             #return and lagged return: dependent variable and regressor
y   = xx[:,1]
x   = xx[:,2]
T   = size(x,1)
#------------------------------------------------------------------------------

xGrid = collect(-10:0.25:10)
h     = 1.5
bHat  = KernRegNormalFn(y,x,xGrid,h,1:T)

xGridU = collect(-9.5:0.25:9.5)
bHatU  = KernRegUniformFn(y,x,xGridU,1)

using PyPlot
figure()
  ha = plot(xGrid,bHat,xGridU,bHatU)
  xlim(minimum(xGrid),maximum(xGrid))
  ylim(-5,5)
  title("Return vs lagged return, kernel regression")
  xlabel("Lagged return")
  ylabel("Return")
  legend(["Normal kernel","Uniform kernel"])
#------------------------------------------------------------------------------

(b,res,) = OlsFn(y,[x.^2 x ones(T,1)])             #rule of thumb bandwidth

sigma = std(res)
gamm  = b[1]
xSort = sort(x)
x_10  = xSort[round(Int64,floor(T*0.1))]          #crude 10th and 90th percentiles
x_90  = xSort[round(Int64,floor(T*0.9))]

h_crude0 = 0.6*sigma^(2/5)*abs(gamm)^(-2/5)*(x_90-x_10)^(1/5)*T^(-1/5)
h_crude  = 0.67
println("initial choice of h: two versions")
println("\n",round([h_crude0 h_crude],3))
#------------------------------------------------------------------------------

if CrossValCalcIt == 1

  println("Cross-validation calculations take some time")

  hM = h_crude*[0.75 1 1.5 2 3 4 5 10]'

  Nh   = length(hM)
  EPEM = fill(NaN,(T,Nh))
  for t = 1:T
    v_No_t = setdiff(1:T,t)
    for j = 1:Nh
      b_t       = KernRegNormalFn(y,x,x[t],hM[j],v_No_t)   #fitted b[x(t)]
      tst       = (y[t] - b_t)^2
      EPEM[t,j] = tst[1]
    end       
  end
  EPE = mean(EPEM,1)'
  println("h and EPE")
  println(round([hM EPE],4))
  

  figure()
    ha = plot(hM,EPE/minimum(EPE))
    xlim(minimum(hM)-0.05,maximum(hM)+0.05)
    ylim(0.99,1.025)
    title("Cross validation simulations, kernel regression")
    xlabel("Bandwidth")
    ylabel("EPE")

end
#------------------------------------------------------------------------------
