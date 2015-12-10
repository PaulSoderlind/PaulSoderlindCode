#------------------------------------------------------------------------------
#  OlsBootstrapEx.m
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
#  Paul.Soderlind@unisg.ch   12 July 2013
#------------------------------------------------------------------------------

include("jlFiles/lagnPs.jl")
include("jlFiles/excise.jl")
include("jlFiles/OlsFn.jl")
#------------------------------------------------------------------------------

#xx   = readdlm("Data/FFmFactorsPs.csv",',',header=true)      
#x    = xx[1]
#Rme  = x[:,2]
#RSMB = x[:,3]                #small minus big firms
#RHML = x[:,4]                #high minus low book-to-market ratio
#Rf   = x[:,5]                    #interest rate
#
#
#x = readdlm("Data/FF25Ps.csv",',')  #no header line: x is matrix     
#R  = x[:,2:end]                    #returns for 25 FF portfolios
#Re = R - repmat(Rf,1,size(R,2))   #excess returns for the 25 FF portfolios
#y = Re(:,[1,25])                 #use only portfolio 1 (small growth) and 25 (large value)
#x = [ones(size(Re,1),1) Rme RSMB RHML]
#------------------------------------------------------------------------------

xx  = readdlm("Data/BondPremiaPs.csv",',',header=true)      
z   = xx[1]
rx  = z[:,2:5]
f   = z[:,6:end]

x = [ones(size(f,1),1) lagnPs(f,12)]
K = size(x,2)

yx = excise([rx[:,4] x])
y  = yx[:,1]
x  = yx[:,2:end]
#------------------------------------------------------------------------------

(bLS,res,yhat,Covb,) = OlsFn(y,x)              #OLS estimate and classical std errors
StdbLS = sqrt(diag(Covb))
println("\n","LS coeffs and std")
println(round([bLS';StdbLS'],3))

T = size(y,1)                 #no. obs and no. test assets
n = size(y,2)               
K = size(x,2)

BlockSize = 10                  #size of blocks
NSim      = 2000                 #no. of simulations

nBlocks = round(Int64,ceil(T/BlockSize))             #number of blocks, rounded up
bBoot   = fill(NaN,(NSim,K*n))                       #vec(b), [beq1 beq2..beqn]
for i = 1:NSim                                       #loop over simulations
  t_i       = rand(1:T,nBlocks,1)                    #nBlocks x 1, random starting row of blocks
  t_i       = broadcast(+,t_i,collect(0:BlockSize-1)')  #nBlocks x BlockSize, each row is a block
  #println(t_i)                                      #uncomment to see which rows that are picked out
  t_i        = reshape(t_i',length(t_i))             #column vector of the blocks
  vv_i       = t_i .> T                             
  t_i[vv_i]  = t_i[vv_i] - T                         #wrap around if index > T
  epsilon    = res[round(Int64,t_i),:]        
  epsilon    = epsilon[1:T,:]                        #get exact sample length
  y_i        = x*bLS + epsilon
  b_i,       = OlsFn(y_i,x)                          #,skips the remaining outputs
  bBoot[i,:] = b_i'                     
end

println("\n","Average bootstrap estimates and bootstrapped std")
println(round([mean(bBoot,1); std(bBoot,1)],3))
#------------------------------------------------------------------------------
