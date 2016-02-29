#-------------------------------------------------------------------------
#  Some finance examples in Julia
#
#
#
#  Paul Söderlind (Paul.Soderlind at unisg.ch), Dec 2015
#-------------------------------------------------------------------------

                                  #do Pkg.add("package name") if needed
using StatsBase, Gadfly, Roots, MathProgBase, Ipopt

set_default_plot_size(20cm, 13cm)


println("\n","------------------------- CAPM regression ----------------","\n")

xx  = readdlm("MyData.csv",',',header=true)
x   = xx[1]
ym  = x[:,1]                 #yearmonth
Rme = x[:,2]                 #market excess return
Rf  = x[:,3]                 #interest rate
R   = x[:,4]                 #return small growth stocks
Re  = R - Rf                 #excess returns
T   = size(Rme,1)

x    = [ones(T) Rme]             #regressors
y    = Re                        #just to get standard OLS notation
b    = x\y                       #OLS
u    = y - x*b                   #residuals
covb = inv(x'x)*var(u)           #cov(b)
stdb = sqrt(diag(covb))          #std(b)
R2   = 1 - var(u)/var(y)

println("OLS coefficients and std")
println(round([b stdb],3))
println("R2: ",round(R2,3))
println("no. of observations: ",T)


println("\n","---------------------- Autocorrelations ------------------","\n")

plags = 1:5
xCorr = autocor(Re,plags)
println("\n lag autocorrr and its t-stat of excess returns")
println(round([plags xCorr sqrt(T)*xCorr],3))


println("\n","------------------------ Value at Risk -------------------","\n")

lambda = 0.95                         #weight on old volatility

Rm     = Rme + Rf                     #equity market return in %
sigma2 = var(Rm)
mu     = mean(Rm)

vol2 = fill(sigma2,T)
for t = 2:T
  vol2[t] = lambda*vol2[t-1] + (1-lambda)*(Rm[t-1]-mu)^2    #RiskMetrics approach
end

VaR95 = -(mu - 1.64*sqrt(vol2))            #VaR at 95% level

YearFrac = floor(ym/100) + (mod(ym,100)-1)/12    #eg 1990.92 for Dec 1990

plot1 = plot(x=YearFrac,y=VaR95,Geom.line,Theme(default_color=colorant"blue"),
             Scale.x_continuous(minvalue=1978,maxvalue=2012),
             Scale.y_continuous(minvalue=0,maxvalue=11),
             Guide.xticks(ticks=[1980,1990,2000,2010]),
             Guide.yticks(ticks=[0,5,10]),
             Guide.title("1-month Value at Risk (95%) for US equity market"),
             Guide.xlabel(""),
             Guide.ylabel("%"))
display(plot1)                              #show plot in browser


println("\n","------------------ Mean-variance frontier ----------------","\n")

Rf    = 0.04                       #some inputs, here riskfree
mu    = [0.125; 0.105; 0.07]       #mean returns
stdd  = [0.129; 0.08; 0.1]         #std
Corr  = [ 1   0.33   0.45;         #correlation matrix
          0.33  1    0.05;
          0.45  0.05  1];
Sigma  = (stdd*stdd').*Corr        #covariance matrix
N      = length(mu)                #number of assets

mu_p = linspace(0,0.15,151)        #average returns to calculate min(Std(R_p)) at

ettor   = ones(N)                  #ettor means ones in Swedish: this a vector of those
Sigma_1 = inv(Sigma)               #mean-variance arithmetic
A       = (ettor'Sigma_1*mu)[1]          # [1] transforms from 1x1 array to scalar
B       = (mu'Sigma_1*mu)[1]
C       = (ettor'Sigma_1*ettor)[1]
D       = B*C - A^2
K       = length(mu_p)

Var_p = fill(NaN,K)                        #to store results in
for i = 1:K                                #loop over mean returns in Std x  mean figure
  g        = ( B*(Sigma_1*ettor) - A*(Sigma_1*mu) )/D
  h        = ( C*(Sigma_1*mu)  - A*(Sigma_1*ettor) )/D
  w_pi     = g + h.*mu_p[i]                        #portfolio weights at mean mu_p[i]
  Var_p[i] = (w_pi'Sigma*w_pi)[1]                  #variance of portfolio with mean mu_p[i]
end

SR_T    = sqrt((mu-Rf)'Sigma_1*(mu-Rf))          #Sharpe ratio of Tangency portfolio
CML_std = sqrt(((mu_p-Rf)./SR_T).^2)             #std according to CLM

plot1 = plot(layer(x=sqrt(Var_p),y=mu_p,Geom.line(preserve_order=true),Theme(default_color=colorant"red")),
             layer(x=CML_std,y=mu_p,Geom.line(preserve_order=true),Theme(default_color=colorant"blue")),
             Guide.title("Mean-variance frontier"),
             Guide.xlabel("Std"),
             Guide.ylabel("mean"))
display(plot1)                                   #show plot in browser


println("\n","---------- Mean-variance frontier w.o. short sales--------","\n")

include("MyFunctions.jl")     #load in file with some functions

Var_no_ss, = MeanVarNoSSPs(mu,Sigma,mu_p)   #see MyFunctions.jl file

plot1 = plot(layer(x=sqrt(Var_p),y=mu_p,Geom.line(preserve_order=true),Theme(default_color=colorant"red")),
             layer(x=CML_std,y=mu_p,Geom.line(preserve_order=true),Theme(default_color=colorant"blue")),
             layer(x=sqrt(Var_no_ss),y=mu_p,Geom.line(preserve_order=true),Theme(default_color=colorant"green")),
             Guide.title("Mean-variance frontier"),
             Guide.xlabel("Std"),
             Guide.ylabel("mean"),
             Guide.manual_color_key(" ",
                  ["Risky only","Risky&riskfree","Risky only (no ss)"],
                  ["red", "blue","green"]))
display(plot1)


println("\n","--------------------- Black-Scholes ----------------------","\n")

function Stdn_cdfPs(a)
#Stdn_cdfPs    Calculates Pr(X<=a) for X which is N(0,1)
  cdf = 0.5*erfc(-a/sqrt(2))
  return cdf
end

function OptionBlackSPs(S,X,T,r,sigma)
#Calculates Black-Scholes european call option price
  d1 = ( log(S./X) + (r+1/2*sigma.^2)*T ) ./ (sigma*sqrt(T))
  d2 = d1 - sigma*sqrt(T)
  c  = S.*Stdn_cdfPs(d1) - X.*exp(-r*T).*Stdn_cdfPs(d2)
  return c
end

sigma = 0.4
c1 = OptionBlackSPs(10,10,0.5,0.1,sigma)
println("\n","call price according to Black-Scholes: ",round(c1,3))

X = linspace(7,13,51)
c = OptionBlackSPs(10,X,0.5,0.1,sigma)
plot1 = plot(x=X,y=c,Geom.line,Theme(default_color=colorant"red"),
             Guide.title("Black-Scholes call option price"),
             Guide.xlabel("strike price"),
             Guide.ylabel("option price"))
display(plot1)

                                #simple (crude) way to solve for implied vol
iv = Roots.fzero(sigma->OptionBlackSPs(10,10,0.5,0.1,sigma)-c1,[-1;1])
println("Implied volatility (should be the same as above): ",round(iv,3))

#  LIFFE Bunds option data, trade date April 6, 1994
X = [                        #strike prices; Mx1 vector
      92.00;  94.00;  94.50;  95.00;  95.50;  96.00;  96.50;  97.00;
      97.50;  98.00;  98.50;  99.00;  99.50;  100.0;  100.5;  101.0;
     101.5;  102.0;  102.5;  103.0;  103.5 ];
C = [                        #call prices; Mx1 vector
      5.13;    3.25;    2.83;    2.40;    2.00;    1.64;    1.31;    1.02;
      0.770;   0.570;   0.400;   0.280;   0.190;   0.130;  0.0800;  0.0500;
      0.0400;  0.0300;  0.0200;  0.0100;  0.0100 ];
F = 97.05                #Forward price
m = 48/365               #time to expiry in years
r = 0.0                  #Interest rate: LIFFE=>no discounting
N = length(X)

iv = fill(NaN,N)
for i = 1:N
  iv[i] = Roots.fzero(sigma->OptionBlackSPs(exp(-m*r)*F,X[i],m,r,sigma)-C[i],[-1;1])
end

println(round([X iv],4))

plot1 = plot(x=X,y=iv,Geom.line,Theme(default_color=colorant"red"),
             Guide.title("Implied volatility, Bunds options April 6, 1994"),
             Guide.xlabel("strike price"),
             Guide.ylabel(" "))
display(plot1)
