function KernRegUniformFn(y,x,xGrid,h)

  Ngrid = length(xGrid)                 #number of grid points
  
  bHat  = fill(NaN,(Ngrid,1))            #y(t) = b[x(t)] + e(t)
  for i = 1:Ngrid                        #loop over elements in xGrid
    zi = (x - xGrid[i])/h
    w  = (abs(zi) .< 0.5) + 0.0         # + 0.0 transforms 'True' into 1
    w  = w/h                            #uniform kernel over xGrid +/- h/2 
    if sum(w)  > 0.0                    #avoids dividing by 0
      bHat[i,:] = sum(w.*y)/sum(w)      #sum over observations (data)
    end
  end

  return bHat

end  