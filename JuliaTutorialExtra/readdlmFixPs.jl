function readdlmFixPs(x,missVal=NaN)
#readdlmFixPs    Change from " " to missVal (defaults to NaN) in arrays
#                imported by readdlm and then convert to float
#
#
#
#
#
#  Usage:    readdlmFixPs!(x[,missVal])
#
#  Input:    x         TxK heterogeneous matrix with some missing values "  "
#            missVal   scalar, to replace missing values in X with
#
#  Output:   y         TxK matrix
#
#
#
#
#
#
#  Paul.Soderlind@unisg.ch
#------------------------------------------------------------------------------

  (T,N) = size(x,1,2)

  y = copy(x)

  for i = 1:T                      #clumsy, but works
    for j = 1:N
      if !isa(x[i,j],Number)
        y[i,j] = missVal
      end
    end
  end

  y = convert(Array{Float64},y)   #convert to float

  return y

end
#------------------------------------------------------------------------------



function readdlmFixPs!(x,missVal=NaN)
#readdlmFixPs    Change from " " to missVal (defaults to NaN) in arrays
#                imported by readdlm and then convert to float
#
#
#
#
#  Usage:    readdlmFixPs!(x[,missVal])
#
#  Input:    x         TxK heterogeneous matrix with some missing values "  "
#            missVal   scalar, to replace missing values in X with
#
#  Output:   changes x
#
#
#
#
#  Paul.Soderlind@unisg.ch
#------------------------------------------------------------------------------

  (T,N) = size(x,1,2)

  for i = 1:T                      #clumsy, but works
    for j = 1:N
      if !isa(x[i,j],Number)
        x[i,j] = missVal
      end
    end
  end

  x = convert(Array{Float64},x)   #convert to float

end
#------------------------------------------------------------------------------
