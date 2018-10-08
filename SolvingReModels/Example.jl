#--------------------------------------------------------------------------
#  Example.jl
#
#  Test program for the functions in ComItAlg.jl, DiscAlg, and SimpRulT.jl.
#  The program also uses Var1SimPs to calculate impulse response functions.
#
#  The program uses a very simplified version of the model in
#  Fuhrer, J. C. (1997), "Inflation/Output Variance Trade-Offs and Optimal
#  Monetary Policy," JMCB, 29, 214-234.
#
#
#
#  This example file and the procedures were created by translating MatLab code
#  to Julia.
#
#
#  Paul Soderlind (Paul.Soderlind@unisg.ch), 20 June 2000, to Julia Jan 2016
#----------------------------------------------------------------------------

using PyPlot                        #comment out this is PyPlot is not installed

using LinearAlgebra
include("Fuhrer1.jl")
include("Var1SimPs.jl")
include("SimpRulT.jl")
include("ComitAlg.jl")
include("DiscAlg.jl")
#-----------------------------------------------------------------------

                               #simplified Fuhrer model, parameter values
a1    = 0.85                   #IS equation
ap    = -0.41
Varey = 0.84^2

w     = 0.65                   #contracting equation
gamm  = 0.002
Varep = 0.19^2

Dbig  = 40                     #arbitrage condition

qy    = 1                      #Loss function
qpi   = 1
qf    = 0.5
bet  = 0.99
                                #Parameters -> system matrices
(A,B,K,Q,U,R) = Fuhrer1(a1,ap,Varey,w,gamm,Varep,Dbig,qy,qpi,qf,bet)
n1 = 3
n2 = 2
k  = 1

Tbig   = 12                       #periods to simulate in VAR
Shock0 = [sqrt(Varep);zeros(2)]   #Impulse response wrt ep
#-----------------------------------------------------------------------
                                  #SIMPLE RULE, u(t)=-Fx(t)

Fy  = 0.5     #Simple decision rule: i(t) = = Fy*y(t) + Fpi*[p(t)-p(t-1)]
Fpi = 1.2
F   = [0  -Fy  -Fpi*4*(1-w)  0  -Fpi*4*w ]

x10 = zeros(n1)                      #initial state vector
SigmaXX = diagm(0=>[Varep,Varey,0])  #covariance matrix of shocks
(M_Simp,C_Simp,J0) = SimpRulT(A,B,Q,R,U,bet,n1,n2,F,SigmaXX,x10,1.0)


x1 = Var1SimPs(M_Simp,Shock0,Tbig)     #VAR of x1(t)
x2 = x1*C_Simp'                        #x2(t)
x  = [x1 x2]                           #x(t) = [x1(t),x2(t)]
uu = -x*F'                             #u(t)
ypii_Simp = [x[:,2]  (4*(1-w)*x[:,3]+4*w*x[:,5]) uu]


println("\nSimple rule: impulse response to a one std of price shock (y,pi,i)")
display(round.([1:Tbig ypii_Simp],digits=3))

#comment out this is PyPlot is not installed
figure()
  plot(1:Tbig,ypii_Simp)
  title("Simple rule: impulse response to a one std of price shock")
  legend(["Output","Inflation","Short interest rate"])
#-----------------------------------------------------------------------
                                  #COMMITMENT

(M_Commit,C_Commit) = ComItAlg(A,B,Q,R,U,bet,n1,n2,1.0)


k      = Var1SimPs(M_Commit,[Shock0;zeros(2)],Tbig) #VAR of x1(t),p2(t)
lambda = k*C_Commit'                          #x2(t),u(t),p1(t)
uu     = lambda[:,n2+1:n2+1]                  #u(t)
x      = [k[:,1:n1] lambda[:,1:n2]]           #x(t) = [x1(t),x2(t)]
ypii_Commit = [x[:,2] (4*(1-w)*x[:,3]+4*w*x[:,5]) uu]

println("\nCommitment: impulse response to a one std of price shock (y,pi,i)")
display(round.([1:Tbig ypii_Commit],digits=3))

#comment out this is PyPlot is not installed
figure()
  plot(1:Tbig,ypii_Commit)
  title("Commitment: impulse response to a one std of price shock")
  legend(["Output","Inflation","Short interest rate"])
#-----------------------------------------------------------------------
                                 #DISCRETION

println("\nPlease wait, discretionary case takes some time to solve")
(M_Disc,C_Disc,V,F) = DiscAlg(A,B,Q,R,U,bet,n1,n2,Matrix(1.0I,n1,n1),zeros(n2,n1),
                               [1e-1;1e-5],0,1,0,1,0,1e+4)

x1 = Var1SimPs(M_Disc,Shock0,Tbig)    #VAR of x1(t)
x2 = x1*C_Disc'                        #x2(t)
x  = [x1 x2]                           #x(t) = [x1(t),x2(t)]
uu = -x1*F'                            #u(t)
ypii_Disc = [x[:,2] (4*(1-w)*x[:,3]+4*w*x[:,5]) uu]

println("\nDiscretion: impulse response to a one std of price shock (y,pi,i)")
display(round.([1:Tbig ypii_Disc],digits=3))

#comment out this is PyPlot is not installed

figure()
  plot(1:Tbig,ypii_Disc)
  title("Discretion: impulse response to a one std of price shock")
  legend(["Output","Inflation","Short interest rate"])
#-----------------------------------------------------------------------

