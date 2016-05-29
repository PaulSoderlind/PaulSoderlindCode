function excise(x)

  vv = vec(any(isnan(x),2))

  if any(vv)              #only keep rows with no NaNs
    x = x[!vv,:]
  end

  return x

end
#------------------------------------------------------------------------------
