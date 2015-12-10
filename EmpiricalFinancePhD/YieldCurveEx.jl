#------------------------------------------------------------------------------
#  YieldCurveEx.jl
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


include("jlFiles/VasicekABFn.jl")
include("jlFiles/VasicekTsCsFn.jl")
include("jlFiles/lagnPs.jl")
using Optim, PyPlot
#-----------------------------------------------------------------------------

xx  = readdlm("Data/USCMRatesPs.csv",',',header=true)      
y   = xx[1]

YearMonth = 1970 + 1/24 + (1/12)*collect(0:size(y,1)-1)

y[:,1] = ((1-0.25*y[:,1]/100).^(-1/0.25) - 1)*100    #discount basis->effective
y[:,2] = ((1-0.50*y[:,2]/100).^(-1/0.5) - 1)*100

m      = [0.25 0.5 1 3 5 7 10]                             #time to maturity (in years)
mMonth = round(Int64,m*12)                            #maturites, in months (integers)
y      = log(1+y/100)                    #continuously compounded interest rates
y      = y/12                            #interest rates per month (the period length of data)

(T,n) = size(y)

vvo = 1                               #indices of yields that are observed without errors
vvu = 2:n                             #indices of yields with observation errors
#------------------------------------------------------------------------------

yo = y[:,vvo]                         #observable yields
yu = y[:,vvu]                         #unobservable yields
nMo = mMonth[vvo]                     #maturity of yo
nMu = mMonth[vvu]                     #maturity of yu

#notice: important that the elements of the parameter vector have similar magnitude
#lambda/100,mu*1200,p = 1 - 2/(1+exp(par(3))) so 5.6 -> p = 0.993,sigma*1200,omega*1200
par0 = [-2.3;10;5.5;0.5;0.8]

                         #just testing the VasicekABFn function
(ao,bo,xt,au,bu,yuHat) = VasicekABFn(1,0.5,0.9,0.02,nMo,nMu,yo)


(MinusLL,yuHat,u,xt) = VasicekTsCsFn(par0,yo,yu,nMo,nMu)

                             #do MLE by optimization with optimize, minimize -sum(LL)
x1 = optimize(par->VasicekTsCsLossFn(par,yo,yu,nMo,nMu),par0)        
par1 = x1.minimum
println("\n par0 and par1")
println(round([par0 par1],3))
                                                #fitted yields at parameter estimates
(MinusLL,yuHat,) = VasicekTsCsFn(par1,yo,yu,nMo,nMu)
yhatTsCs = [yo yuHat]                      #also yo
                                            #a and b at estimated parameters, average yield curve
(ao,bo,xt,au,bu,yuHatAvg) = VasicekABFn(par1[1]*100,par1[2]/1200,1 - 2/(1+exp(par1[3])),par1[4]/1200,
                                        nMo,nMu,mean(yo))
yuHatAvg = [mean(yo) yuHatAvg]             #also yo
#------------------------------------------------------------------------------

figure()
  yy = y[:,[3,6]] - yhatTsCs[:,[3,6]]
  ha = plot(YearMonth,12*yy)
  xlim(1970,2014)
  ylim(-0.04,0.04)
  title("Pricing errors, Vasicek (TS and CS)")
  legend(["1 year";"7 years"])

figure()
  ha = plot(m',12*[yuHatAvg' mean(y,1)'])
  xlim(0,11)
  ylim(0.05,0.08)
  title("Avg yield curve, Vasicek (TS and CS)")
  legend(["at average x(t)";"data"])
  xlabel("Maturity (years)")
#------------------------------------------------------------------------------


bOls = [ones(T,1) y[:,1]]\y                 #LS of yields on short yield
alfa = bOls[1,:]
beta = bOls[2,:]
                         #rescaling Vasicek au and bu to be comaparable with OLS on yo
bu_b = bu/bo             #y = au + bu*xt, but xt = (yo-ao)/bo
au_b = au - bu/bo*ao     #y = au + bu*(yo-ao)/bo = y = au - bu*ao/bo + bu/bo*yo
bu_b = [1;bu_b]'         #also yo
au_b = [0;au_b]'

figure()
  ha = plot(m',[alfa' au_b']*1200)
  xlim(0,maximum(m)+1)
  title("Intercept times 1200, Vasicek and OLS")
  xlabel("maturity (years)")
  legend(["Vasicek";"OLS"])

figure()
  ha = plot(m',[beta' bu_b'])
  xlim(0,maximum(m)+1)
  title("Slope on short rate, Vasicek and OLS")
  xlabel("maturity (years)")
  legend(["Vasicek";"OLS"])
#------------------------------------------------------------------------------