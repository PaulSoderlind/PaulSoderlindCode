"""
    FindNoNaNPs(Keepdim,x...)

Find rows (if Keepdim==1) which have no NaNs in other dimensions (eg. in no columns).
Returns a BitArray(1) vector. Set Keepdim=2 if we should instead look for NaNs along
rows (and other dimensions).


Paul.Soderlind@unisg.ch
"""
function FindNoNaNPs(Keepdim,x...)

  xNum  = length(x)
  T     = size(x[1],Keepdim)                 #length of output

  #xDims = ndims(x[1])                        #not good enough
  xDims = maximum([ndims(x[i]) for i=1:xNum])
  dims  = setdiff(1:xDims,Keepdim)           #dimensions to check

  vvM = BitArray(T,xNum)
  for i = 1:xNum                             #loop over inputs
    vvM[:,i] = any(isnan.(x[i]),dims)
  end

  vvb = broadcast(!,vec(any(vvM,2)))         #.! in 0.6

  return vvb

end
#------------------------------------------------------------------------------
