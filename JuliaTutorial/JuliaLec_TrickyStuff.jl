#  This file highlights some tricky aspects of Julia (from the perspective of
#  a matlab user)
#
#
#
#
#
#  Paul Söderlind (Paul.Soderlind at unisg.ch), December 2015
#-------------------------------------------------------------------------



println("\n\n------- A = B means that A and B are the same forever-------","\n")
A = [2 2]
B = A
C = sum(B)
D = A + 0
println("(the arrays) old A,B,C,D: ",A," ",B," ",C," ",D)
A[2] = 3
println("new A,B,C,D after after changing element A[2]: ",A," ",B," ",C," ",D)
println("\nNotice that when A is changed, then it carries over to B since
A and B are one and the same. Actually, if you instead changed B, then it
would carry over to A as well. In contrast, C and D are not changed when A is:
they are not the same as A.")


println("\n\n------------------- 1x1 arrays are not scalars  ------------------")
A = [1 2]
b = [3]
println("\nA and b: ",A," ",b)
println("You cannot do A + b, if A is a Txn array and b is a 1x1 array.
Instead, use A .+ b: ",A .+ b)
println("This works since b is expanded ('broadcasted') to have the same dimension as A")

A[2] = b[1]
println("\nYou cannot do A[2] = b, if A is a Txn array and b is a 1x1 array.
Instead use A[2] = b[1]: "," new A ",A)

c = ones(2)
d = [10;11]
println("\nc and d: ",c," ",d)
A[2] = (c'd)[1]
println("You cannot do A[2] = c'd, if A is a Txn array and c and d are vectors.
Instead use A[2] = (c'd)[1]: "," new A ",A)


println("\n\n------------------Creating variables in loop------------------")

for i = 1:5
  Tor = cos(i)
end
println("\nvariables CREATED in a loop are not visible outside the loop
Trying to print Tor after the loop would give an error message (try it)")

println("\n","In contrast, variables CHANGED in a loop are visible outside the loop")
Oden = []
for i = 1:5
  Oden = cos(i)
end
println("Oden ",round(Oden,4))


println("\n\n----------------------Adding rows to an array----------------------")

A =  [1 11]
B =  [3 13]
println("\nA and B: ",A," ",B)
println("\nTo append B at the end of A, you have to use [A;B],
doing A[2,:] = B does not work ")
println([A;B])


println("\n\n---------------------------Cell arrays------------------------")

println("\nCreate Any[x1,x2,...]")
A = Any[rand(3,2),"A nice dog",27]
println("\nThe array A: ")
println(A)
println("\nAlternatively, you can use A= cell(3) and fill as A[3] = 27")
println("\nElement 3,2 of A[1]: ",A[1][3,2])



