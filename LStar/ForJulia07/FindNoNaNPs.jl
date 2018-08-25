#------------------------------------------------------------------------------
"""
    FindNoNaNPs(Keepdim,x...)

Find rows (if Keepdim==1) which have no NaNs missing in other dimensions (eg. in no columns).

# Input
- `Keepdim::Int`: 1 if check rows, 2 if check columns, etc
- `z::Array`: one or several numerical arrays

# Output
- `vvb::BitArray(T)`: vector, element t is true if row (if Keepdim==1) t has no NaN or missing

# Notice
- Set Keepdim=2 if we should instead look for NaNs/missings along rows (and other dimensions).


Paul.Soderlind@unisg.ch

"""
function FindNoNaNPs(Keepdim,x...)

  xNum  = length(x)
  T     = size(x[1],Keepdim)                 #length of output

  xDims = maximum([ndims(x[i]) for i=1:xNum])
  dims  = setdiff(1:xDims,Keepdim)           #dimensions to check

  vvM = falses(T,xNum)
  for i = 1:xNum                             #loop over inputs
    vvM[:,i] = any(isnan.(x[i]) .| ismissing.(x[i]),dims=dims)
  end

  vvb = vec(.!any(vvM,dims=2))      #rows witout NaN/missing in any of the x matrices

  return vvb

end
#------------------------------------------------------------------------------
