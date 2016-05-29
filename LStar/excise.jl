function excise(x)

  vv = vec(any(isnan(x),2))

  if any(vv)              #only keep rows with no NaNs
    x = x[!vv,:]
  end

  return x

end
#------------------------------------------------------------------------------

function excise2mPs(x,y)

  vv = vec(any(isnan([x y]),2))

  if any(vv)              #only keep rows with no NaNs
    x = x[!vv,:]
    y = y[!vv,:]
  end

  return x,y

end
#------------------------------------------------------------------------------

function excise3mPs(x,y,z)

  vv = vec(any(isnan([x y z]),2))

  if any(vv)              #only keep rows with no NaNs
    x = x[!vv,:]
    y = y[!vv,:]
    z = z[!vv,:]
  end

  return x,y,z

end
#------------------------------------------------------------------------------

function excise4mPs(x,y,z,w)

  vv = vec(any(isnan([x y z w]),2))

  if any(vv)              #only keep rows with no NaNs
    x = x[!vv,:]
    y = y[!vv,:]
    z = z[!vv,:]
    w = w[!vv,:]
  end

  return x,y,z,w

end
#------------------------------------------------------------------------------
