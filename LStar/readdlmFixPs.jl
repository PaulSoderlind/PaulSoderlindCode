function readdlmFixPs(x,missVal=NaN)
#readdlmFixPs    Change from " " to missVal (default to NaN) in arrays 
#                imported by readdlm and then convert to float


  (T,N) = size(x)

  for i = 1:T                      #clumsy, but works 
    for j = 1:N
      if !isa(x[i,j],Number)
        x[i,j] = missVal
      end 
    end
  end
  
  x = convert(Array{Float64},x)   #convert to float

  return x
end   



#--------------------------OLD APPROACH----------------------------------------
                                    
#  x2 = fill(NaN,(T,N))             #x2=x (eats memory...), but with NaNs and Float64
#  for i = 1:T                      #clumsy, but works 
#    for j = 1:N                   
#      if isa(x[i,j],Number)
#        x2[i,j] = x[i,j]
#      end  
#    end
#  end
#  return x2
#------------------------------------------------------------------------------