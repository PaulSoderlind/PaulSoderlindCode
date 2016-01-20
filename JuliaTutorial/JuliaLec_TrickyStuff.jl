#  This file highlights some tricky aspects of Julia
#
#
#
#
#
#  Paul Söderlind (Paul.Soderlind at unisg.ch), December 2015
#-------------------------------------------------------------------------

println("\n","--------------------------------------------------------------","\n")
println("\n","------------------------ Handle NaNs -------------------------","\n")

xx     = readdlm("JuliaLecData.csv",',',header=true)     #loading into matrix
x      = xx[1]
x[1,4] = NaN                            #set to NaN
z      = x[:,3:4]

println("First five lines of z:")
println(z[1:5,:])

if any(isnan(z))       #check if any NaNs
  println("\nz has some missing values. Cannot do calculations on those.")
end
println("If data has NaNs, then most calculations give a NaN.")
println("Standard deviations of Tx2 matrix z: ",round(std(z,1),3))

println("\n","Cut out row 1, since the NaN is there")
println("Standard deviations of (T-1)x2 matrix z[2:end,:]: ",round(std(z[2:end,:],1),3))

vv = !any(isnan(z),2)             #a more automatic approach: find rows with no NaN
z = z[vv,:]                       #only keep rows with no NaNs. clumsy but works
println("\n","Standard deviations of automatically trimmed matrix: ",round(std(z,1),3))
#-------------------------------------------------------------------------

println("\n","--------------------------------------------------------------","\n")
println("\n","--- A = B means that A and B are one and the same (forever)---","\n")
A = [2 2]
B = A
C = sum(B)
D = A/2
println("(the arrays) A,B,C,D before: ",A,B,C,D)
A[2] = 3
println("A,B,C,D after after changing A[2]: ",A,B,C,D)
println("\nNotice that when A is changed, then it carries
         over to B since A and B are one and the same")
println("Actually, if you instead changed B, then it carries over to A")
println("\nIn contrast, C and D are not changed when A is: they are not the same as A")

println("\n","------------- A = B + 0 means A and B are not one and the same --------","\n")
A = [2 2]
B = A + 0
println("A,B before: ",A,B)
A[2] = 3
println("A,B after after changing A[2]: ",A,B)
#-------------------------------------------------------------------------

println("\n","--------------------------------------------------------------","\n")
println("\n","----------------- 1x1 arrays are not scalars  ----------------")
A = [1 2]
b = [3]                                  #A + b does not work
println("(1)--- cannot do A + b, if A is Txn array and b is 1x1 array  ---","\n")
println("instead, use A .+ b: ",A .+ b)  #works since b is expanded ('broadcasted')
                                         #to have same dimension as A

println("\n","(2)-- cannot do A[2] = b, if A is Txn array and b is 1x1 array --","\n")
A = [1 2]
println("old A and b: ",A," ",b)
b = [3]
A[2] = b[1]                 #A[1] = b does not work, A[2] = b[1] does
println("instead use A[2] = b[1]: "," new A ",A)
#-------------------------------------------------------------------------

println("\n","--------------------------------------------------------------","\n")
println("\n","-- variables CREATED in a loop are not visible outside the loop  ----")

for i = 1:5
  Tor = cos(i)
end
println("\nTrying to print Tor would give an error message")

println("\n","-- variables CHANGED in a loop are visible outside the loop  ----")
Tor = []
for i = 1:5
  Tor = cos(i)
end
println(round(Tor,4))
#-------------------------------------------------------------------------





