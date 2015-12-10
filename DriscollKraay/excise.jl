function excise(x)

  #vv = find(!any(isnan(x),2))       
  vv = !any(isnan(x),2)
  z  = x[vv,:]           #only keep rows with no NaNs. clumsy but works

  return z

end