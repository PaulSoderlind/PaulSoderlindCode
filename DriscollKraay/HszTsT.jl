#------------------------------------------------------------------------------
#  HszTsT.jl
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
#  Paul.Soderlind@unisg.ch   Nov 2015
#------------------------------------------------------------------------------



include("HDirProdPs.jl")
include("excise.jl")
include("HszDk5cPs.jl")
#------------------------------------------------------------------------------
 
yh = readdlm("yh.csv",',')             #load from csv files instead
x  = readdlm("x.csv",',')
z  = readdlm("z.csv",',')

T = size(yh,1)
N = size(yh,2)


nGroups = size(z,2)                  #no. of groups
vvM = z .== 1                        #booleans for membership in group
#------------------------------------------------------------------------------


yj = fill(NaN,(T,nGroups))                      #EW portfolios of each group
for i = 1:nGroups
  vvi     = vvM[:,i]
  yj[:,i] = mean(yh[:,vvi],2)
end

println("b2 and Std(b2), LS group by group")
for i = 1:nGroups
  b2  = x\yj[:,i] 
  println(round(b2,3))
end
#------------------------------------------------------------------------------

z4 = fill(NaN,(T,N,size(z,2)))                  #create TxNxK z matrix 
for t = 1:T
  z4[t,:,:] = z
end

fnOutput = HszDk5cPs(yh,x,z4,true,2,true)

println("\n","b4b StdWhite, StdDKlags")
println(round([fnOutput[1] fnOutput[3] fnOutput[9]],5))
println("\n","R2: ",round(fnOutput[6],5))
#------------------------------------------------------------------------------