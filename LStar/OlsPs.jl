function OlsPs(y,x,ExciseIt=false,UnExciseIt=false,SkipCovbIt=false)
#OlsPs    LS of y on x, for one dependent variable or SURE with same regressors
#
#
#
#  Usage:    (b,res,yhat,Covb,R2,T) = OlsPs(y,x[,ExciseIt[,UnExciseIt[,SkipCovbIt]]]])
#
#  Input:    y            Tx1 or Txn matrix of the dependent variables
#            x            Txk matrix of regressors (including deterministic ones)
#            ExciseIt     (optional) scalar, if 1, excise([y,x])
#            UnExciseIt   (optional) scalar, 1 if putting residuals and fitted values in Txn matrix
#            SkipCovbIt   (optional) scalar, 1 if skipping calculation of Covb (sensitive inversion)
#
#  Output:   b            kxn matrix, regression coefficients
#            res          Tx1 or Txn matrix, residuals y - yhat
#            yhat         Tx1 or Txn matrix, fitted values x*b
#            Covb         matrix, covariance matrix of vec(b) =[beq1;beq2;...]
#            R2           1xn vector, R2 values
#            T            scalar, number of obs
#
#
#
#
#  Calls on:  excisePs
#
#
#
#  Paul.Soderlind@unisg.ch, Jan 2001, Jun 2005, to Julia Oct 2015
#----------------------------------------------------------------------------

  T = size(y,1)
  n = size(y,2)

  if ExciseIt
    (yx,vvNaN,vvNaNR,vvNoNaNR) = excisePs([y x])
    y = yx[:,1:n]
    x = yx[:,n+1:end]
  end
  yx = nothing

  Ty = size(y,1)
  Tx = size(x,1)
  k  = size(x,2)

  if Tx != Ty
    error("y and x must have same number of observations")
  end
  if any(isnan([y x]))
    error("NaN in x or y")
  end

  if Ty >= k
    b      = x\y
    yhat2  = x*b
    res2   = y - yhat2
    Covres = cov(res2)*(Ty-1)/Ty
    if SkipCovbIt
      Covb = NaN
    else
      Covb   = kron(Covres,inv(x'x))
    end
    R2 = 1 - var(res2,1)./var(y,1)
  else
    b      = NaN
    yhat2  = NaN
    res2   = NaN
    Covb   = NaN
    R2     = NaN
  end

  if UnExciseIt        #put fitted value and residual in Tyxn matrix
    yhat             = fill(NaN,(T,n))
    res              = fill(NaN,(T,n))
    yhat[vvNoNaNR,:] = yhat2
    res[vvNoNaNR,:]  = res2
  else
    yhat = yhat2
    res  = res2
  end

  return b,res,yhat,Covb,R2,T

end
#------------------------------------------------------------------------------
