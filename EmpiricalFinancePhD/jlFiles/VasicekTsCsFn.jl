function VasicekTsCsFn(par,yo,yu,nMo,nMu)
#VasicekTsCsFn    Loss function for estimating the Vasicek model using both
#                 time series and cross sectional information
#
#------------------------------------------------------------------------------

  J = length(nMu)
  T = size(yo,1)
  
  lambda  = par[1]*100
  mu      = par[2]/1200
  p       = 1 - 2/(1+exp(par[3]))          #par(3) = log((1+p)/(1-p))
  s2      = (par[4]/1200)^2
  omega_i = abs(par[5]/1200)
  
  if length(nMo) != 1
    error("yo must be a single yield")
  end
  
  (ao,bo,xt,au,bu,yuHat) = VasicekABFn(lambda,mu,p,sqrt(s2),nMo,nMu,yo)
  
  Et1xt    = (1-p)*mu + p*lagnPs(xt)    #E(t-1)x(t)
  Et1xt[1] = mu
  Et1yo    = ao + Et1xt*bo              #E(t-1)yo(t)
  v        = yo - Et1yo                 #Tx1, forecast error of yo
  S        = bo's2*bo                   #variance of forecast error of yo
  S_1      = inv(S*1000)*1000           #1000/1000 improves the precision a bit
  LLo_t    = -0.5*log(pi) - 0.5*log(det(S)) - 0.5*sum(v*S_1.*v,2)   #Tx1
  
  u        = yu -  yuHat                 #TxL, fitted errors of yu
  Omega    = eye(J)*omega_i^2
  Omega_1  = inv(Omega)                  #covariance matrix of u
  LLu_t    = -0.5*J*log(pi) - 0.5*log(det(Omega)) - 0.5*sum(u*Omega_1.*u,2)  #Tx1

  LL_t    = LLo_t + LLu_t                #Tx1, log likelihood(t), sum of TS and CS
  MinusLL = -sum(LL_t[2:end])            #to be minimized

  return MinusLL,yuHat,u,xt

end
#------------------------------------------------------------------------------

function VasicekTsCsLossFn(par,yo,yu,nMo,nMu)

  MinusLL, = VasicekTsCsFn(par,yo,yu,nMo,nMu)

  return MinusLL

end  
#------------------------------------------------------------------------------
