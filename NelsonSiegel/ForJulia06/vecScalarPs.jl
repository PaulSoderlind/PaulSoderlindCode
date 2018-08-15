#------------------------------------------------------------------------------
function vecPs(x::Number)
#  vecPs
    y = [x]
    return y
end
function vecPs(x::AbstractArray)
#  Notice:  changing y (after fn) will affect x (before fn)
    y = vec(x)
    return y
end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function ScalarOf1x1Ps(x)
  if length(x) == 1
    y = x[1]
  else
    error("x has more than one element")
  end
  return y
end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function isScalarPs(x)
  y = typeof(x) <: Number
  return y
end
#------------------------------------------------------------------------------

