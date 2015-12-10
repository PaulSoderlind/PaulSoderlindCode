function HszDkFn(y,x,z)
#HszDkFn   LS and Driscoll-Kray standard errors for panel, assuming x(t,i) = x(t) * z(i)
#
#
#  Paul.Soderlind@unisg.ch   Oct 2015
#------------------------------------------------------------------------------

  T = size(y,1)
  N = size(y,2)
  K = size(x,2)*size(z,2)
    
  Sxx = 0.0
  Sxy = 0.0
  for t = 1:T                          #OLS by looping over t
    y_t  = y[t,:]'                     #dependent variable, Nx1
    x0_t = repmat(x[t,:],N,1)          #factors, NxK
    x_t  = HDirProdFn(z,x0_t)          #effective regressors, z is NxKz, x_t is NxK
    Sxx  = Sxx + x_t'x_t/(T*N)         #building up Sxx and Sxy
    Sxy  = Sxy + x_t'y_t/(T*N)
  end
  
  theta = Sxx\Sxy
  
  s2     = 0.0
  omegaj = zeros(K,K)
  for t = 1:T                          #Covariance matrix by looping over t
    y_t  = y[t,:]'                     #create y_t and x_t (again)
    x0_t = repmat(x[t,:],N,1)
    x_t  = HDirProdFn(z,x0_t)
    e_t  = y_t - x_t*theta             #residuals in t
    h_t  = (x_t'e_t)'/N                #moment conditions in t (divided by N)
    omegaj = omegaj + h_t'h_t          #building up covariance matrix
    s2     = s2 + sum(e_t.^2)/N^2
  end
  Shat = omegaj/T^2                     #estimate of S
  s2   = s2/T^2
    
  zx_1  = inv(Sxx)
  CovDK = zx_1 * Shat * zx_1'                     #covariance matrix, DK
  stdDK = sqrt( diag(CovDK) )                     #standard errors, DK
  
  CovLS = zx_1 * s2                               #covariance matrix, LS iid
  stdLS = sqrt(diag(CovLS))                       #standard errors, LS iid
  
  return theta,CovDK,CovLS

end
