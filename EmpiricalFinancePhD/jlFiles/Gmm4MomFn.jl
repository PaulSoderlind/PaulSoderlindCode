function Gmm4MomFn(par,x)

  mu = par[1]
  s2 = par[2]
  
  g = [(x-mu) ((x-mu).^2-s2) ((x-mu).^3) ((x-mu).^4-3*s2^2)]    #Tx4
  
  gbar = mean(g,1)'                         #4x1

  return g,gbar

end  
#------------------------------------------------------------------------------

function Gmm4MomLossFn(par,x,W=1)

  (g,gbar) = Gmm4MomFn(par,x)

  Loss = gbar'W*gbar                #to be minimized
  Loss = Loss[1]

  return Loss  

end    
