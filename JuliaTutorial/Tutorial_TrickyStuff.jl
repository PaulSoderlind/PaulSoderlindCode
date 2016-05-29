#  This file highlights some tricky aspects of Julia (from the perspective of
#  a matlab user)
#
#
#
#
#
#  Paul Söderlind (Paul.Soderlind at unisg.ch), December 2015
#-------------------------------------------------------------------------



print_with_color(:red,"\n\n---for arrays, A = B means that A and B stay the same---\n")
A = [2 2]
B = A
C = sum(B)
D = A + 0
println("old A,B,C,D: ",A," ",B," ",C," ",D)
A[2] = 3
println("after changing element A[2]: ")
println("new A,B,C,D: ",A," ",B," ",C," ",D)
println("\nNotice that when A is changed, then it carries over to B since
A and B are one and the same. Actually, if you instead changed B, then it
would carry over to A as well. In contrast, C and D are not changed when A is:
they are not the same as A.")


print_with_color(:red,"\n\n----changing array in function can have effect outside function----\n")
function f1(A)
  A[1:end] = A[1:end]*2
  return A
end
function f2(A)
  A = A*2
  return A
end
function f3(A)
  #B = deepcopy(A)                    #works too
  #B[1:end] = B[1:end]*2
  #return B
  A = A + 0
  A[1:end] = A[1:end]*2
  return A
end
x0 = [1;2]
x1 = deepcopy(x0)
x2 = deepcopy(x0)
x3 = deepcopy(x0)
y1 = f1(x1)
y2 = f2(x2)
y3 = f3(x3)
println("original x: ",x0,", new x after calling f1(): ",x1,", f1(x1): ",y1)
println("original x: ",x0,", new x after calling f2(): ",x2,", f2(x2): ",y2)
println("original x: ",x0,", new x after calling f3(): ",x3,", f3(x3): ",y3)
println("\nNotice that when individual ELEMENTS of an array are changed inside a
function, then this carries over to the array used in the function call. This is true
also when we change all individual elements (as in f1()). It is not true when
we work on the entire array (as in f2()) or change its shape. The solution
to the problem with f1() is to do as in f3(): work on a copy of the input array.")


print_with_color(:red,"\n\n------------------- 1x1 arrays are not scalars  ------------------\n")
A = [1 2]
b = [3]
println("A and b: ",A," ",b,"\n")
println("You cannot do A + b if A is a Txn array and b is a 1x1 array.
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


print_with_color(:red,"\n\n------------------Creating variables in loop------------------\n")
for i = 1:5
  Tor = cos(i)
end
println("Variables CREATED in a loop are not visible outside the loop
Trying to print Tor after the loop would give an error message")

println("\n","In contrast, variables CHANGED in a loop are visible outside the loop")
Oden = Float64[]
for i = 1:5
  Oden = cos(i)
end
println("Oden ",round(Oden,4))


print_with_color(:red,"\n\n----------------------Adding rows to an array----------------------\n")
A =  [1 11]
B =  [3 13]
println("A and B: ",A," ",B)
println("\nTo append B at the end of A, you have to use [A;B]. A[2,:] = B does not work ")
println([A;B])


print_with_color(:red,"\n\n---------------------------Cell arrays------------------------\n")
println("Create Any[x1,x2,...]")
A = Any[[11 12;21 22],"A nice dog",27]
println("\nThe array A: ")
println(A)
println("\nAlternatively, you can use A = cell(3) and fill as A[3] = 27")
println("\nElement 2,2 of A[1]: ",A[1][2,2])

