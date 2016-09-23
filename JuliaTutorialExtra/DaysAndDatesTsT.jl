#--------------------------building a calendar---------------------------------

#dNb      = Array{Date}(1)                   #clumsy way to build monthly calendar
#dNb[1]   = Date(2014,1,1)
#for i = 1:11
#  push!(dNb,dNb[end] + Dates.Month(1))
#end
                                             #better way to build monthly calendar
dNb = collect(Date(2014,1,1):Dates.Month(1):Date(2014,12,1))

println("\nMy calendar and day of week:")
for i = 1:length(dNb)
  println(dNb[i]," ", Dates.dayofweek(dNb[i]))
end

#--------------------------convert from datetime to date-----------------------

#Background: importing xls sheets gives DateTime

xlsDate = DateTime(2014,1,1)
jlDate  = Date(xlsDate)

println("\nDateTime from from xls and converted to Date:")
println(xlsDate," ",jlDate)

#if xlsDate is a vector, loop over the elements


#--------------------------convert from matlab's datenum to date-----------------------

#Background: in matlab datenum(2016,3,31) gives 736420.0

dNml = 736420.0
dN = round(Int,dNml) - 366            #convert to integer, subtract 366
dN = Date(Dates.rata2datetime(dN))
println("\nmatlab datenum and Julia Date:")
println(dNml," ",dN)

#if dNml is a vector, loop over the elements
#------------------------------------------------------------------------------
