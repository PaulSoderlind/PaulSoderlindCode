################################################################################
println("\n---------------------------------------------------------------------")
println("-------------------saving and loading csv----------------------------")

A = reshape(1:40,5,8)
writedlm("NewCsvFile.csv",A,',')
println("NewCsvFile.csv has been created in the current folder")

A2 = readdlm("NewCsvFile.csv",',')
println("\nA and A2")
println(A)
println(A2)

################################################################################
println("\n\n---------------------------------------------------------------------")
println("-----------------------saving and loading jld file-------------------")

using JLD       #see https://github.com/JuliaLang/JLD.jl for more examples

B   = reshape(1:40,5,8)
B27 = 1
save("NewJldFile.jld","B",B,"B27",B27)
println("NewJldFile.jld has been created in the current folder")

B2 = load("NewJldFile.jld","B")
println("\nB from jld file is")
println(B2)
                                                 #alternative approach
xx = load("NewJldFile.jld")                      #load all variables into dictionary
println("\nAll variables in the JLD file:")
println(keys(xx))
B3 = xx["B"]
println("\nB from jld file is (2nd way of loading)")
println(B3)


################################################################################
println("\n\n---------------------------------------------------------------------")
println("------------------saving and loading matlab mat file-------------------")
using MAT
#see https://github.com/simonster/MAT.jl for more examples

C   = reshape(1:40,5,8)
C27 = 1
fh = matopen("NewMatFile.mat","w")               #create a mat file
write(fh,"C",C)
write(fh,"C27",C27)
close(fh)
println("\nNewMatFile.mat has been created in the current folder")

fh = matopen("NewMatFile.mat")                   #open the mat file
println("\nVariables in mat file: ",names(fh))
C2 = read(fh,"C")                                #read variable A
close(fh)                                        #close the mat file
println("\nC from mat file is ")
println(C2)
                                                 #alternative approach
xx = matread("NewMatFile.mat")                   #read whole mat file into xx
C3 = xx["C"]                                     #xx is a dictionary
println("\nC from mat file is (2nd way of loading) ")
println(C3)


################################################################################
println("\n\n---------------------------------------------------------------------")
println("--------------loading csv with headers, fixing missing values--------")


xx = readdlm("Data/loadCsvTsT_Data.csv",',',header=true)     #xx is Any[] array
x = xx[1]
println("\nx")
println(x)

include("readdlmFixPs.jl")                       #function for fixing missing values
y = readdlmFixPs(x)                              #fixing missing value at x[1,4]
println("\nafter fix")
println(y)


################################################################################
println("\n\n---------------------------------------------------------------------")
println("-----------------------loading xlsfile-------------------------------")

using ExcelReaders
#see https://github.com/davidanthoff/ExcelReaders.jl for more examples
#Notice: you need python's xlrd libarary for this to work.

println("\n------------Approach 1: readxl-------------------")

data1 = readxl("Data/readXlsTsT_Data.xlsx","Data!B2:C11")       #reading using range
println("\nRaw input")
println(data1)
x1 = convert(Array{Float64},data1)
println("\nNumeric part after conversion")
println(x1)

println("\n------------Approach 2: readxlsheet--------------")

data2  = readxlsheet("Data/readXlsTsT_Data.xlsx","Data",skipstartrows=1)  #reading all columns
println("\nRaw input")
println(data2)
x2 = convert(Array{Float64},data2[:,2:end])
println("\nNumeric part after conversion to floating point numbers")
println(x2)

x2[x2 .== -999.99] = NaN                     #converting to NaNs
println("\nNumeric part after changing -999.99 to NaN")
println(x2)

T  = size(data2,1)                           #converting datetimes to dates
dN = Array{Date}(T)
for t = 1:T
  dN[t] = Date(data2[t,1])
end
println("\nBoth dates (after conversion) and numeric part")
println([dN x2])
