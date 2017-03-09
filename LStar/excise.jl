#------------------------------------------------------------------------------
function excise(x)
#notice: if !any(vv), then this does NOT create a copy, so changing the output will
#will change the input

  vv = vec(any(isnan.(x),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)               #use .~vv in 0.6?
    x = x[vvb,:]                        #use view(x,vvb,:) in the future?
  end

  return x

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function exciseCopy(x)

  vv = vec(any(isnan.(x),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)
    x1 = x[vvb,:]
  else
    x1 = deepcopy(x)
  end

  return x1

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise1mPs(x)

  vv = vec(any(isnan.(x),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)
    x = x[vvb,:]
  end

  return x,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise1mCopyPs(x)

  vv = vec(any(isnan.(x),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)
    x1 = x[vvb,:]
  else
    x1 = deepcopy(x)
  end

  return x1,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise2mPs(x,y)

  vv = vec(any(isnan.([x y]),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)
    x = x[vvb,:]
    y = y[vvb,:]
  end

  return x,y,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise2mCopyPs(x,y)

  vv = vec(any(isnan.([x y]),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)
    x1 = x[vvb,:]
    y1 = y[vvb,:]
  else
    x1 = deepcopy(x)
    y1 = deepcopy(y)
  end

  return x1,y1,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise3mPs(x,y,z)

  vv = vec(any(isnan.([x y z]),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)
    x = x[vvb,:]
    y = y[vvb,:]
    z = z[vvb,:]
  end

  return x,y,z,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise3mCopyPs(x,y,z)

  vv = vec(any(isnan.([x y z]),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)
    x1 = x[vvb,:]
    y1 = y[vvb,:]
    z1 = z[vvb,:]
  else
    x1 = deepcopy(x)
    y1 = deepcopy(y)
    z1 = deepcopy(z)
  end

  return x1,y1,z1,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise4mPs(x,y,z,w)

  vv = vec(any(isnan.([x y z w]),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)
    x = x[vvb,:]
    y = y[vvb,:]
    z = z[vvb,:]
    w = w[vvb,:]
  end

  return x,y,z,w,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise4mCopyPs(x,y,z,w)

  vv = vec(any(isnan.([x y z w]),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)
    x1 = x[vvb,:]
    y1 = y[vvb,:]
    z1 = z[vvb,:]
    w1 = w[vvb,:]
  else
    x1 = deepcopy(x)
    y1 = deepcopy(y)
    z1 = deepcopy(z)
    w1 = deepcopy(w)
  end

  return x1,y1,z1,w1,vv

end
#------------------------------------------------------------------------------
