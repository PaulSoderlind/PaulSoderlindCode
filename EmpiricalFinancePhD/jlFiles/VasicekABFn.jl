function VasicekABFn(lambda,mu,rho,sigma,nMo,nMu,yo)

  nMax = maximum([nMo;nMu])    #longest maturity to calculate (a,b) for
  Nvec = (1:nMax)'
  
  A = fill(NaN,(1,nMax))              #recursive solution of AR(1) model
  B = fill(NaN,(1,nMax))
  B[1] = 1 + 0*rho
  A[1] = 0 + 0*(1-rho)*mu - (lambda+0)^2*sigma^2/2
  for n = 2:nMax
    B[n] = 1 + B[n-1]*rho
    A[n] = A[n-1] + B[n-1]*(1-rho)*mu - (lambda+B[n-1])^2*sigma^2/2
  end
  
  a  = A./Nvec
  b  = B./Nvec
  
  ao = a[nMo]                           #required maturities
  bo = b[nMo]
  
  xt    = inv(bo)*(yo-ao)               #value of state variable x(t)
  T     = size(yo,1)
  yoHat = ao + xt*bo                    #Tx1, should be the same as yo
  
  au    = a[nMu]
  bu    = b[nMu]
  yuHat = repmat(au',T,1) + xt*bu'     #Txn, fitted values of yu

  return ao,bo,xt,au,bu,yuHat

end
