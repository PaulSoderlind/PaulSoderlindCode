#----------------------------------------------------------------------------
"""
    OlsPs(y,x,ExciseIt=false,UnExciseIt=false,SkipCovbIt=false)

Do OLS of y on x, for one dependent variable or SURE with same regressors

# Usage
`(b,res,yhat,Covb,R2a,T) = OlsPs(y,x,ExciseIt,UnExciseIt,SkipCovbIt)`

# Input
- `y::Array`:         Tx1 or Txn matrix of the dependent variables
- `x::Array`:         Txk matrix of regressors (including deterministic ones)
- `ExciseIt::Bool`:   (optional) bool, true to doing excise([y,x])
- `UnExciseIt::Bool`: (optional) bool, true to put residuals and fitted values in Txn matrix
- `SkipCovbIt::Bool`: (optional) bool, true to skip calculation of Covb (sensitive inversion)

# Output
- `b::Array`:     kxn matrix, regression coefficients
- `res::Array`:   Tx1 or Txn matrix, residuals y - yhat
- `yhat::Array`:  Tx1 or Txn matrix, fitted values x*b
- `Covb::Array`:  matrix, covariance matrix of vec(b) = [beq1;beq2;...]
- `R2a::Array`:   1xn vector, R2 values
- `Ty::Int`:      scalar, number of obs (after excise)

# Requires
- FindNoNaNPs

"""
function OlsPs(y,x,ExciseIt=false,UnExciseIt=false,SkipCovbIt=false)

  (T,n) = (size(y,1),size(y,2))

  if ExciseIt
     vvNoNaNR = FindNNPs(y,x)                #find rows with no NaNs
    (y,x)     = (y[vvNoNaNR,:],x[vvNoNaNR,:])
  end

  (Ty,Tx,k) = (size(y,1),size(x,1),size(x,2))

  if Tx != Ty
    error("y and x must have same number of observations")
  end
  if any(isnan.([y x]))
    error("NaN in x or y")
  end
  if ndims(x) == 1               #T vector to (T,1) matrix, could also do reshape(x,T,1)
    x = x[:,:]
  end

  if Ty >= k
    b = x\y
    #b = (x'x)\(x'y)
    #ndims(x) == 1 ? yhat2  = x*b[1] : yhat2  = x*b  #if? then: else. to handle x is vector
    yhat2  = x*b
    res2   = y - yhat2
    Covres = cov(res2)*(Ty-1)/Ty
    if SkipCovbIt
      Covb = NaN
    else
      Covb   = kron(Covres,inv(x'x))
    end
    R2a = 1 .- var(res2,dims=1)./var(y,dims=1)
  else
    (b,yhat2,res2,Covb,R2a) = (NaN,NaN,NaN,NaN,NaN)
  end

  if UnExciseIt && ExciseIt   #put fitted value and residual in Tyxn matrix
    (yhat,res)       = [fill(NaN,(T,n)) for i=1:2]
    yhat[vvNoNaNR,:] = yhat2
    res[vvNoNaNR,:]  = res2
  else
    yhat = yhat2
    res  = res2
  end

  if n == 1                #if only one regression
    R2a = R2a[1]
  end

  return b,res,yhat,Covb,R2a,Ty

end
#------------------------------------------------------------------------------
