#------------------------------------------------------------------------------
"""
    FindNNPs(x...;Keepdim=1)

Find rows (if Keepdim==1) which have no NaNs missing in other dimensions (eg. in no columns).

# Input
- `z::Array`: one or several numerical arrays
- `Keepdim::Int`: (keyword) 1 if check rows, 2 if check columns, etc

# Output
- `vvb::BitArray(T)`: vector, element t is true if row (if Keepdim==1) t has no NaN or missing

# Notice
- Set Keepdim=2 if we should instead look for NaNs/missings along rows (and other dimensions).
- For heterogenous arrays like `x=[x1,x1]`, use `FindNNPs(x...)`

Paul.Soderlind@unisg.ch

"""
function FindNNPs(x...;Keepdim=1)

  Base.require_one_based_indexing(x)

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


#------------------------------------------------------------------------------
"""
    FindNNDictPs(D...)

Find rows that have no missing values in any of the arrays in the dictionaries D...,
assuming all arrays have the same number of rows.

# Input
- `D::Dict`:    One (D) or several (D1,D2) dictionaries  with N different Tx? arrays.

"""
function FindNNDictPs(D...)

  N    = length(D)

  TM = hcat( [hcat([size(x,1) for x in values(D[i])]...) for i=1:N]... )  #size(x,1)
  all(TM .== TM[1]) ? T = TM[1] : error("different number of rows")

  vv2 = falses(T,N)
  for i = 1:N
    vv       = hcat([FindNNPs(x) for x in values(D[i])]...)
    vv2[:,i] = vec(all(vv,dims=2))
  end

  vvb = vec(all(vv2,dims=2))

  return vvb

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
"""
    DictSelectRowsPs(D,vv)

Create new dictionary where we pick rows of D[i] to vv ()

"""
function DictSelectRowsPs(D,vv)
  Dnew = Dict()
  for (k,v) in D
    ndims(v) == 1 ? Dnew[k] = v[vv] : Dnew[k] = v[vv,:]
  end
  return Dnew
end
#------------------------------------------------------------------------------


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
