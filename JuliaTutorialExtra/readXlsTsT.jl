
using ExcelReaders
#see https://github.com/davidanthoff/ExcelReaders.jl for more examples
#Notice: you need python's xlrd libarary for this to work.


data1 = readxl("readXlsTsT_Data.xlsx","Data!B2:C11")       #reading in using range
println("\nRaw input")
println(data1)
x1 = convert(Array{Float64},data1)
println("\nNumeric part after conversion")
println(x1)
println("------------------------------------------------------------")

data2  = readxlsheet("readXlsTsT_Data.xlsx","Data",skipstartrows=1)  #reading in all columns
println("\nRaw input")
println(data2)
x2 = convert(Array{Float64},data2[:,2:end])
println("\nNumeric part after conversion")
println(x2)
println("------------------------------------------------------------")

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
