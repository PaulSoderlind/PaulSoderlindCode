#------------------------------------------------------------------------------
function excise(x)

  vv = vec(any(isnan(x),2))

  if any(vv)              #only keep rows with no NaNs
    x = x[!vv,:]
  end

  return x

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function exciseCopy(x)

  vv = vec(any(isnan(x),2))

  if any(vv)              #only keep rows with no NaNs
    x1 = x[!vv,:]
  else
    x1 = deepcopy(x)
  end

  return x1

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise1mPs(x)

  vv = vec(any(isnan(x),2))

  if any(vv)              #only keep rows with no NaNs
    x = x[!vv,:]
  end

  return x,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise1mCopyPs(x)

  vv = vec(any(isnan(x),2))

  if any(vv)              #only keep rows with no NaNs
    x1 = x[!vv,:]
  else
    x1 = deepcopy(x)
  end

  return x1,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise2mPs(x,y)

  vv = vec(any(isnan([x y]),2))

  if any(vv)              #only keep rows with no NaNs
    x = x[!vv,:]
    y = y[!vv,:]
  end

  return x,y,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise2mCopyPs(x,y)

  vv = vec(any(isnan([x y]),2))

  if any(vv)              #only keep rows with no NaNs
    x1 = x[!vv,:]
    y1 = y[!vv,:]
  else
    x1 = deepcopy(x)
    y1 = deepcopy(y)
  end

  return x1,y1,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise3mPs(x,y,z)

  vv = vec(any(isnan([x y z]),2))

  if any(vv)              #only keep rows with no NaNs
    x = x[!vv,:]
    y = y[!vv,:]
    z = z[!vv,:]
  end

  return x,y,z,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise3mCopyPs(x,y,z)

  vv = vec(any(isnan([x y z]),2))

  if any(vv)              #only keep rows with no NaNs
    x1 = x[!vv,:]
    y1 = y[!vv,:]
    z1 = z[!vv,:]
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

  vv = vec(any(isnan([x y z w]),2))

  if any(vv)              #only keep rows with no NaNs
    x = x[!vv,:]
    y = y[!vv,:]
    z = z[!vv,:]
    w = w[!vv,:]
  end

  return x,y,z,w,vv

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function excise4mCopyPs(x,y,z,w)

  vv = vec(any(isnan([x y z w]),2))

  if any(vv)              #only keep rows with no NaNs
    x1 = x[!vv,:]
    y1 = y[!vv,:]
    z1 = z[!vv,:]
    w1 = w[!vv,:]
  else
    x1 = deepcopy(x)
    y1 = deepcopy(y)
    z1 = deepcopy(z)
    w1 = deepcopy(w)
  end

  return x1,y1,z1,w1,vv

end
#------------------------------------------------------------------------------
