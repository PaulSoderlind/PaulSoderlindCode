#------------------------------------------------------------------------------
#  QuadFormEx.m
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
#  Paul.Soderlind@unisg.ch   05 July 2013
#------------------------------------------------------------------------------


xx   = readdlm("Data/FFmFactorsPs.csv",',',header=true)      
x    = xx[1]
x    = x[:,2:end]/100
Rme  = x[:,1]
RSMB = x[:,2]                #small minus big firms
RHML = x[:,3]                #high minus low book-to-market ratio
#------------------------------------------------------------------------------

v = [Rme RSMB RHML]
v = v - repmat(mean(v,1),size(v,1),1)

S   = cov(v)
S_1 = inv(S)

T = size(v,1)

tic()
Qa = fill(NaN,(T,1))
for t = 1:T                   #do by a loop
  v_t = v[t,:]'               #v_t is Kx1
  Qa[t] = (v_t'S_1*v_t)[1]    #1x1 matrix to scalar... 
end
toc()

tic()
Qb = sum(v*S_1.*v,2)      #do it by the sum(,2)
toc()

println("max(abs(Qa-Qb)): ",maximum(abs(Qa-Qb)))
