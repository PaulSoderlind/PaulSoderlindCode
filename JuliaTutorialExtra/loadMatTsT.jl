using MAT
#see https://github.com/simonster/MAT.jl for more examples

fh = matopen("loadMatTsT_Data.mat")         #open the mat file

println("\nVariables in mat file: ",names(fh))

A = read(fh,"A")           #read variable A
println("\nA is: ")
println(A)

xx = read(fh)              #read whole mat file into xx
Ab = xx["A"]               #using a dictionary instead
println("\nAb is: ")
println(Ab)

close(fh)                  #close the mat file
