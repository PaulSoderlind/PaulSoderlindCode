
z = [1 NaN;2 12;3 13]                 #a matrix with NaNs
println("z ",z)

if any(isnan(z))                      #check if any NaNs
  println("\nz has some missing values. Cannot do calculations on those.")
end
println("If data has NaNs, then most calculations give a NaN.
Column means of Tx2 matrix z: ",mean(z,1))

println("Column means of (T-1)x2 matrix z[2:end,:]: ",mean(z[2:end,:],1))

vv = !any(isnan(z),2)             #a more automatic approach
z = z[vv,:]                       #keep only rows with no NaNs
println("Column means of automatically trimmed matrix: ",mean(z,1))
