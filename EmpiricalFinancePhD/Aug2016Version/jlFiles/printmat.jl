#------------------------------------------------------------------------------
function printmat(x,width=10,prec=3,Numfmt="f",NoPrinting=false)
#printmat   Prints all elements of matrix with a predefined formatting.
#           Uses the package Formatting (import Formatting)
#
#
#
#  Input:      x           numerical matrix to print
#              width       (optional) scalar, width of printed cells. [10]
#              prec        (optional) scalar, precision of printed cells. []
#              Numfmt      (optional) string, format of numerical cells, ["f"]
#              NoPrinting  (optional) bool, true: no printing, just return formatted string
#
#  Output:     str         (if NoPrinting) string, (otherwise nothing)
#
#
#  Uses:  import Formatting
#
#
#  Paul.Soderlind@unisg.ch, May 2016
#------------------------------------------------------------------------------

  if isa(x,String)                      #strings need special treatment
    str = string(x,"\n")
    if NoPrinting
      return str
    else
      println(str)
      return nothing
    end
  end

  numdims = ndims(x)

  if numdims > 3
    warn("$numdims > 3")
    return nothing
  end

  eltype_x = eltype(x)
                                                        #floats or intergers
  if any([eltype_x <: AbstractFloat eltype_x <: Signed eltype_x <: Unsigned])
    if eltype_x <: Integer
      prec = 0
    end
    fmt2a = string("%",string(width),".",string(prec),Numfmt)  #fmt numerical data
  else                                 #assume strings
    fmt2a = string("%",string(width),"s")
  end

  (m,n,p) = size(x,1,2,3)

  str = ""
  for k = 1:p                          #loop over dim 3
    s = ""
    if p > 1
      s = s * "x[:,:,$k]\n"
    end
    for i = 1:m
      for j = 1:n
        s = s * Formatting.sprintf1(fmt2a,x[i,j,k])
      end
      s = s * @sprintf("%s\n","")
    end
    str = str * s
  end

  if NoPrinting                              #no printing, just return str
    return str
  else                                       #print, return nothing
    println(str)
    return nothing
  end

end
#------------------------------------------------------------------------------
