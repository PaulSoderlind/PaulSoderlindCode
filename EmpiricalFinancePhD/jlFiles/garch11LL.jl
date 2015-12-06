function garch11LL(parm::Vector,yx)

  y = yx[:,1]               #split up yx
  x = yx[:,2:end]
  
  T = size(x,1)
  k = size(x,2)
  b     = parm[1:k]         #mean equation, y = x'*b
  omega = abs(parm[k+1])    #GARCH(1,1) equation: s2(t) = omega + alpha*u(t-1)^2 + beta1*s2(t-1)
  alpha = abs(parm[k+2])
  beta1 = abs(parm[k+3])    #do not use label 'beta' since that is an ml function
  
  yhat = x*b
  u    = y - yhat
  s2_0 = var(u)                                 #var(u,1) gives a matrix, var(u) a scalar
                                                
  s2    = fill(NaN,(T,1))
  s2[1] = omega + alpha*s2_0 + beta1*s2_0        #simple, but slow apparoach,
  for t = 2:T                                    #using filter() is quicker
    s2[t] = omega + alpha*u[t-1]^2 + beta1*s2[t-1]
  end
  
  LL    = -(1/2)*log(2*pi) - (1/2)*log(s2) - (1/2)*(u.^2)./s2
  LL[1] = 0               #effectively skip the first observation
  
                          #output scalar: overall LL function value
  LL = -sum(LL)           #to minimize -sum(LL), notice
  
  
  return LL,s2,yhat

end
#------------------------------------------------------------------------------

function garch11LLLoss(parm::Vector,yx)

  LL, = garch11LL(parm::Vector,yx)

  return LL
  
end
#------------------------------------------------------------------------------
   