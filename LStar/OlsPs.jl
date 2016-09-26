function OlsPs(y,x,ExciseIt=false,UnExciseIt=false,SkipCovbIt=false)
#OlsPs    LS of y on x, for one dependent variable or SURE with same regressors
#
#
#
#  Usage:    (b,res,yhat,Covb,R2a,T) = OlsPs(y,x[,ExciseIt[,UnExciseIt[,SkipCovbIt]]]])
#
#  Input:    y            Tx1 or Txn matrix of the dependent variables
#            x            Txk matrix of regressors (including deterministic ones)
#            ExciseIt     (optional) bool, true to doing excise([y,x])
#            UnExciseIt   (optional) bool, true to put residuals and fitted values in Txn matrix
#            SkipCovbIt   (optional) bool, true to skip calculation of Covb (sensitive inversion)
#
#  Output:   b            kxn matrix, regression coefficients
#            res          Tx1 or Txn matrix, residuals y - yhat
#            yhat         Tx1 or Txn matrix, fitted values x*b
#            Covb         matrix, covariance matrix of vec(b) =[beq1;beq2;...]
#            R2a          1xn vector, R2 values
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

  (T,n) = size(y,1,2)

  if ExciseIt
    (yx,_,_,vvNoNaNR) = excisePs([y x])      #does not affect calling scope
    y = yx[:,1:n]
    x = yx[:,n+1:end]
    (yx,_) = (nothing,nothing)
  end

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
    R2a = 1 - var(res2,1)./var(y,1)
  else
    (b,yhat2,res2,Covb,R2a) = (NaN,NaN,NaN,NaN,NaN)
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

  return b,res,yhat,Covb,R2a,T

end
#------------------------------------------------------------------------------
