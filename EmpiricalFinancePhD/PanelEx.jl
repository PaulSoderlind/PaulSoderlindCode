#------------------------------------------------------------------------------
#  PanelEx.m
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

include("jlFiles/OlsFn.jl")
include("jlFiles/HszDkFn.jl")
include("jlFiles/HDirProdFn.jl")
#------------------------------------------------------------------------------

ER1 = readdlm("Data/PPM_ER1.csv",',')                   #load data from csv files
ER2 = readdlm("Data/PPM_ER2.csv",',')
ER  = [ER1;ER2]
#ER = randn(2354,2637)           #uncomment this line (and comment the previous 3 lines)
                                 #if you do not have ER1.csv and ER2.csv

Factors   = readdlm("Data/PPM_Factors.csv",',')         #no header line: x is matrix
Investors = readdlm("Data/PPM_N_Changes.csv",',')
N_Changes = Investors[:,1]

(T,N) = size(ER)
D = N_Changes .> 50                      #logical dummies: [active]
#------------------------------------------------------------------------------


alphaM = fill(NaN,N)                                #individual alphas
for i = 1:N
   b, = OlsFn(ER[:,i],[Factors ones(T,1)])
   alphaM[i] = b[end]
end
println("\nAverage annualised alphas for the two groups")
println(round([mean(alphaM[~D]) mean(alphaM[D])]*252,3))

PortfER      = fill(NaN,(T,2))     #create portfolios as average across individuals
PortfER[:,1] = mean(ER[:,~D],2)   #Tx1, portfolio return = average individual return
PortfER[:,2] = mean(ER[:,D],2)

println("\n\nLS, group by group")
Avg = mean(PortfER,1)*252          #average excess return on portfolios
Std = std(PortfER,1)*sqrt(252)
SR  = Avg./Std
(b,res,yhat,Covb) = OlsFn(PortfER,[ones(T,1) Factors])
println("\nStats for the portfolios: Avg, Std, SR, alpha")
display(round([Avg' Std' SR' b[1:1,:]'*252],3))

R       = [1 0 0 0 -1 0 0 0]                       #testing if alpha(1) = alpha(2)
a_diff  = R*vec(b)                                 #b(:) = vec(b)
tstatLS = a_diff/sqrt(R*Covb*R')
println("\ndiff of annual alpha (inactive - 51+), tstat (LS)")
println(round([a_diff*252 tstatLS],3))
#------------------------------------------------------------------------------

println("\n\npanel")
(theta,CovDK,CovLS) = HszDkFn(ER,[ones(T,1) Factors],[~D D]+0.0)

R       = [1 0 0 0 -1 0 0 0]                #testing if alpha(1) = alpha(2)
a_diff  = R*theta
tstatLS = a_diff/sqrt(R*CovLS*R')
tstatDK = a_diff/sqrt(R*CovDK*R')
println("\nLS of panel, diff of annual alpha (inactive - 51+), tstat (LS), tstat (DK)")
println(round([a_diff*252 tstatLS tstatDK],3))
#------------------------------------------------------------------------------
