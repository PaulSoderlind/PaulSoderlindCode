#------------------------------------------------------------------------------
#  OlsEx.jl
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

include("jlFiles/NWFn.jl")
include("jlFiles/OlsFn.jl")
include("jlFiles/Ols2Fn.jl")
include("jlFiles/OlsDiagnosticsFn.jl")
include("jlFiles/excise.jl")
#include("jlFiles/lagnPs.jl")

using StatsBase, Distributions
#------------------------------------------------------------------------------

xx   = readdlm("Data/FFmFactorsPs.csv",',',header=true)      
x    = xx[1]
ym   = x[:,1]                                      #[yearmonth]
x    = x[:,2:end]/100
Rme  = x[:,1]
RSMB = x[:,2]                #small minus big firms
RHML = x[:,3]                #high minus low book-to-market ratio

Y = Rme
X = [ones(size(Rme,1),1) RSMB RHML]

(T,K) = size(X)
S_xx = 0.0
S_xy = 0.0
for t = 1:T
  x_t = X[t,:]'            #x_t is 2x1
  y_t = Y[t,:]'
  S_xx = S_xx + x_t*x_t'/T   #2x2
  S_xy = S_xy + x_t*y_t/T    #2x1
end
b1 = inv(S_xx)*S_xy          #OLS coeffs, version 1

b2 = inv(X'X)*X'Y            #OLS coeffs, version 2

b3 = X\Y                     #OLS coeffs, version 3

println("\n b1, b2 and b3")
println(round([b1 b2 b3],3))

b = X\Y
u = Y - X*b              #residuals
g = X.*repmat(u,1,K)     #moment conditions
println("\n avg moment conditions")
println(round(mean(g,1),3))

S = NWFn(g,1)            #Newey-West covariance matrix
D = -X'X/T
V = inv(D'inv(S)*D)     #Cov(sqrt(T)*b)

println("\n b and std(b)")
println(round(b3,3))
println(round(sqrt(diag(V/T)),3))

(b4,res,yhat_,CovbLS_,R2_,T_,CovbNW4) = Ols2Fn(Y,X,1)
println("\n OLS with NW standard errors")
println(round(b4,3))
println(round(sqrt(diag(CovbNW4)),3))

R = [0 1 0;               #testing if b(2)=0 and b(3)=0
     0 0 1]
a = [0;0]
Gamma = R*V*R'
test_stat = (R*b-a)'inv(Gamma/T)*(R*b-a)
println("\n test-statictic and 10% critical value of chi-square(2)")
println([round(test_stat,3) 4.61])

(AutoCorr,DW,BoxPierce,White,Regr) = OlsDiagnosticsFn(Y,X,u,2)     #diagnostics
println("\n","diagnostics with std (and df): AutoCorr,DW,BoxPierce,White,Regr")
println("AutoCorr")
println(round(AutoCorr,3))
println("DW: ",round(DW[1],3))
println("BoxPierce: ", round(BoxPierce,3))
println("White: ",round(White,3))
println("Test of all slopes: ",round(Regr,3))
#------------------------------------------------------------------------------
