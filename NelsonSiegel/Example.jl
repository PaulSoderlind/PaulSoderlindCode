#------------------------------------------------------------------------------
#  Example.jl
#
#
#
#
#
#
#
#
#
#  Paul.Soderlind@unisg.ch, April 2002, to Julia Nov 2015
#------------------------------------------------------------------------------

using Dates, Optim

include("jlFiles/BondPricePs.jl")
include("jlFiles/BondNSxPs.jl")
include("jlFiles/printmat.jl")


ytmLoss  = 0           #0: mimimize squared price errors: 1: minimize squared ytm errors
  weightLoss = 1       #1: 1/mat as weights for (fitted price-price)^2 in loss fn (only if ytmLoss==0)
#------------------------------------------------------------------------------

                                 #Swedish bond data for 29 Dec 1993

#coupons in %/yr
c  = [ 0;        0;         0;         0;         11.5;      10.75;
       11;       13;        10.25;     6;         9 ]

#time to maturity in year
tm = [ 0.00274;  0.21096;   0.46027;   0.88219;   1.67397;   3.06849;
       5.06301;  7.46027;   9.34795;   11.11507;  15.30685 ]

#interest rates.  Bills: simple rates in %/yr; bonds: yield to maturity in %/yr
y  = [ 7.75;     6.835;    6.655;     6.41;      6.215;     6.195;
       6.41;     6.755;    7.01;      7.21;      7.325 ]

n = length(y)       #number of bonds

Base.require_one_based_indexing(tm,y)
#----------------------------------------------------------------------------

#transform the data

c  = c/100 # -> column vector, coupons and yields as 0.05 rather than 5
y  = y/100

vvc    = c .== 0                #if bill, change from simple to effective rate
y[vvc] = (1 .+ tm[vvc].*y[vvc]).^(1.0./tm[vvc]) .- 1

P = fill(NaN,length(y))         #calculate bond prices
for i = 1:n
  local ti
  ti   = mod(tm[i],1):tm[i]
  P[i] = BondPricePs(y[i],c[i],ti)
end
#----------------------------------------------------------------------------

#estimating parameters in extended Nelson&Siegel model, restricted so b0 + b1 = log(1+y[1])

parX0 = [0.1045;-0.03;-0.0562;1.2;0;0.5]       #starting guess

if ytmLoss == 1                                #loss(ytm)
  NSXbR = BondNSxEstPs(parX0,y,tm,c,log(1+y[1]),1)
else
  if weightLoss == 1                           #loss(P/maturity)
    NSXbR = BondNSxEstPs(parX0,P,tm,c,log(1+y[1]),0,1.0./tm)
  else                                         #loss(P)
    NSXbR = BondNSxEstPs(parX0,P,tm,c,log(1+y[1]))
  end
end

println("\nParameter estimates: ")
printTable(NSXbR,[""],["b0","b1","b2","tau","b3","tau2"])
#----------------------------------------------------------------------------

ytmx = fill(NaN,n)            #model implied ytm
for i = 1:n
  local ti,d,Qx
  ti      = mod(tm[i],1):tm[i]
  d       = BondNSxPs(ti,NSXbR...)[3]    #... expands into NSXbR[1],NSXbR[2],...
  Qx      = sum(d.*c[i]) + d[end]     #model implied bond price
  #println(ti)
  ytmx[i] = BondYieldToMatPs(Qx,c[i],ti,1,1,0.05,1e-7)[1]
end

println("\nActual and model ytm, %: ")
printTable([y ytmx]*100,["ytm actual","ytm model"],string.(tm),width=15,cell00="maturity")
#----------------------------------------------------------------------------

#calculate model implied rates (spot, forward, yield to maturity) to plot

tmFig      = [1e-8;0.1:0.1:16]           #maturities to plot
(shx,fhx,) = BondNSxPs(tmFig,NSXbR...)
shx        = exp.(shx) .- 1     #effective interest rate
fhx        = exp.(fhx) .- 1

println("\nmodel spot and forward rates, %: ")
printTable([shx fhx]*100,["spot","forward"],string.(tmFig),width=15,cell00="maturity")
#----------------------------------------------------------------------------

#Comment out this if you do not have PyPlot installed.

using PyPlot
#PyPlot.svg(true)           #for ipynb notebooks
close("all")
                                       #plotting
figure()
  plot(tm,c,"+",tm,y,"s",tmFig,shx,"b--",tmFig,fhx,"r-",tm,ytmx,".")
  title("Swedish Interest Rates 29 Dec 1993")
  xlabel("Years to Maturity")
  legend(["Coupon rate","Yield to maturity","Estimated spot rate",
         "Estimated forward rate","Estimated yield to maturity"],loc=1)
  ylim(0.04,0.14)
  #display(gcf())            #uncomment in VsCode
#----------------------------------------------------------------------------
