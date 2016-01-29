

println("\n---------------------loading csv, fixing missing values--------------")
include("readdlmFixPs.jl")

xx = readdlm("Data/loadCsvTsT_Data.csv",',',header=true)     #works, but gives x[1,4] = "  "
x = xx[1]
println("\nHeaders: ",xx[2])
println("\nx")
println(x)

y = readdlmFixPs(x)
println("\nafter fix")
println(y)

println("\n-----------------------------------------------------------------")
println("----------------------loading matlab mat file----------------------")

using MAT
#see https://github.com/simonster/MAT.jl for more examples

println("\n------------Approach 1: with matopen-------------")

fh = matopen("Data/loadMatTsT_Data.mat")         #open the mat file

println("\nVariables in mat file: ",names(fh))

A = read(fh,"A")           #read variable A
println("\nA is: ")
println(A)

xx = read(fh)              #read whole mat file into xx
Ab = xx["A"]               #using a dictionary instead
println("\nAb is: ")
println(Ab)

close(fh)                  #close the mat file
xx = nothing

println("\n----------Approach 2: with matread-------------")

xx = matread("Data/loadMatTsT_Data.mat")     #read whole mat file into xx
Ab = xx["A"]
println("\nAb is: ")
println(Ab)

println("\n-----------------------------------------------------------------")
println("-------------------------loading xlsfile---------------------------")

using ExcelReaders
#see https://github.com/davidanthoff/ExcelReaders.jl for more examples
#Notice: you need python's xlrd libarary for this to work.

println("\n------------Approach 1: readxl-------------------")

data1 = readxl("Data/readXlsTsT_Data.xlsx","Data!B2:C11")       #reading in using range
println("\nRaw input")
println(data1)
x1 = convert(Array{Float64},data1)
println("\nNumeric part after conversion")
println(x1)

println("\n------------Approach 2: readxlsheet--------------")

data2  = readxlsheet("Data/readXlsTsT_Data.xlsx","Data",skipstartrows=1)  #reading in all columns
println("\nRaw input")
println(data2)
x2 = convert(Array{Float64},data2[:,2:end])
println("\nNumeric part after conversion")
println(x2)

x2[x2 .== -999.99] = NaN                     #converting to NaNs
println("\nNumeric part after changing -999.99 to NaN")
println(x2)

T  = size(data2,1)                           #convertind^g datetimes to dates
dN = Array(Date,T)
for t = 1:T
  dN[t] = Date(data2[t,1])
end
println("\ndate part after convering days, together with numeric part")
println([dN x2])
#------------------------------------------------------------------------------
