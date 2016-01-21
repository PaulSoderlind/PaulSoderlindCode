
include("readdlmFixPs.jl")

xx = readdlm("loadCsvTsT_Data.csv",',',header=true)     #works, but gives x[1,4] = "  "
x = xx[1]
println("\nHeaders: ",xx[2])
println("\nx")
println(x)

y = readdlmFixPs(x)
println("\nafter fix")
println(y)

