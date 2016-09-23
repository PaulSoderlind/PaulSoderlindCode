import Formatting
include("printmat.jl")            #just function for prettier matrix printing


z = [1 NaN;2 12;3 13]                 #a matrix with NaNs
println("z: ")
printmat(z)

if any(isnan(z))                      #check if any NaNs
  println("z has some NaNs")
end

println("\nIf data has NaNs, then most calculations give a NaN. For, instance,
the column means of Tx2 matrix z: ")
printmat(mean(z,1))

vv = vec(!any(isnan(z),2))             #an automatic approach
z2 = z[vv,:]                           #keep only rows with no NaNs
println("z2: a new matrix where all rows with any NaNs have been pruned:")
printmat(z2)
