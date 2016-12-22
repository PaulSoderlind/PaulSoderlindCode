function KernRegNormalFn(y,x,xGrid,h,vv)

  Ngrid = length(xGrid)                  #number of grid points

  bHat = fill(NaN,Ngrid)                 #y(t) = b[x(t)] + e(t)

  for i = 1:Ngrid                        #loop over elements in xGrid
    zi        = (x - xGrid[i])/h
    w         = exp(-zi.^2/2)./(h*sqrt(2*pi))       #gaussian kernel, with "std" = h
    bHat[i]   = sum(w[vv].*y[vv])/sum(w[vv])        #sum over observations (data)
  end

  return bHat

end