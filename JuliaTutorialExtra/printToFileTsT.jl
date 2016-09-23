import Formatting
include("printmat.jl")            #a function for prettier matrix printing


fh = open("NewTxtFile.txt", "w")        #open the file, "w" for writing

println(fh,"testing")
println(fh,"more testing")

x = rand(1:7,(4,3))
println(fh,printmat(x,10,3,"f",true))   #to pretty print the matrix

close(fh)                               #close the file

println("newfile.txt has been created in the current folder")
