function HszDk5dwPs(y,x,z,yhatQ=false,m=0,ScaleByNtQ=0,vvzx=[],wM=[])
#HszDk5dwPs LS and Driscoll-Kray standard errors for unbalanced panel, assuming
#           x are the same across individuals, while z are time-varying individual
#           characteristics. The effective regressors are kron(z,x).
#
#
#
#
#  Usage: fnOutput = HszDk5dwPs(y,x,z[,yhatQ[,m[,ScaleByNtQ[,vvzx[,wM]]]]])
#
#  Input:    y             TxN matrix with the dependent variable, y(t,i) is for period t, individual i
#            x             TxK matrix with K factors that are common for all investors, should
#                          at least contain Tx1 matrix ones(T)
#            z             TxNxL matrix with L (time-varying) individual characteristics, should
#                          at least contain TxNx1 array of ones, ones(T,N,1)
#            yhatQ         (optional) scalar, 1: generate and report fitted values
#            m             (optional), scalar, number of lags in covariance estimation [0]
#            ScaleByNtQ    (optional), scalar, 1: scales all moment conditions by N(t), not by N;
#                          2: scale all moment conditions by wM; [0]
#                          can be used to replicate a portfolio (Calendar time) approach
#            vvzx          (optional), vector of indices of colums in HDirProdPs(z_t,x0_t) to keep
#            wM            (optional) TxN matrix, weights on observation (t,i)
#
#  Output:   fnOutput      named tuple with:
#              theta         (K*L)x1 vector, LS estimates of regression coeefficients on kron(z,x)
#              stdDK         (K*L)x1 vector, Driscoll-Kraay standard errors
#              stdW          (K*L)x1 vector, White's standard errors
#              CovDK         (K*L)x(K*L) matrix, Driscoll-Kraay covariance matrix
#              yhat          TxN matrix with fitted values
#              R2a           scalar, (pseudo-) R2
#              CovW          covariance matrix, White's
#              CovDKj        covariance matrix, DK with lags
#              stdDKj        standard errors, DK with lags
#
#
#  Notice:   (a) the effective regressors are kron(z,x). For instance, with z = [1,z1] and
#                x = [1,x1,x2,x3], we have [1,x1,x2,x3,z1,z1*x1,z1*x2,z1*x3].
#
#
#
#  Uses:     FindNoNaNPs and HDirProdPs
#
#
#
#
#  Paul.Soderlind@unisg.ch   May 2010, to Julia Oct 2015
#------------------------------------------------------------------------------

  if isempty(vvzx)
    vvzx = 1:(size(x,2)*size(z,3))
  end

  (T,N) = (size(y,1),size(y,2))
  K     = size(x,2)
  L     = size(z,3)
  KL    = sum(vvzx .>0)

  msg = " "
  xx  = zeros(KL,KL)                   #Sum[x(t)*x(t)',t=1:T]
  xy  = zeros(KL)                      #Sum[x(t)*y(t),t=1:T]

  Nb  = zeros(Int,T)                   #effective number of obs, after pruning NaNs
  for t = 1:T                             #loop over time
    x_t  = HDirProdPs(reshape(z[t,:,:],N,L),x[t:t,:])  #effective regressors, z_t is NxL, x_t is 1xK
    x_t  = x_t[:,vvzx]
    vv   = FindNoNaNPs(1,y[t,:],x_t)
    N_t  = sum(vv)
    if ScaleByNtQ == 1
      Nb[t] = N_t
      w_t   = ones(N_t)
    elseif ScaleByNtQ == 2
      Nb[t] = N
      w_t   = wM[t,vv]
    else
      Nb[t] = N
      w_t   = ones(N_t)
    end
    if N_t > 0                           #don't accumulate [] to xx and xy (generates [])
      (y_t,x_t) = (y[t,vv],x_t[vv,:])    #pruning NaNs
      y_t = y_t .* w_t                   #put weights on observation i (in t), y_t .* w_t
      lmul!(Diagonal(w_t),x_t)           #x_t = w_t .* x_t
      xx .= xx + x_t'x_t/Nb[t]
      xy .= xy + x_t'y_t/Nb[t]
    end
    msg = IterationPrintPs(t,T,msg,100)
  end

  Tb    = sum(Nb .> 0)                    #number of effective time periods
  xx    = xx/Tb
  xy    = xy/Tb
  theta = xx\xy                           #ols estimates, solves xx*theta-xy=0 in MSE sense

  if yhatQ
    yhat = fill(NaN,(T,N))
  else
    yhat = Float64[]
  end
  omega0DK = zeros(KL,KL)               #DK, lag 0
  omega0W  = zeros(KL,KL)               #White's
  omegajDK = zeros(KL,KL,m)             #DK, lags 1 to m
  h_tLag   = zeros(m,KL)                #lag1;lag2;...,lagm
  for t = 1:T                           #loop over time
    x_t    = HDirProdPs(reshape(z[t,:,:],N,L),x[t:t,:])
    x_t    = x_t[:,vvzx]
    yhat_t = x_t*theta
    r_t    = y[t,:] - yhat_t
    vv     = FindNoNaNPs(1,r_t,x_t)
    N_t    = sum(vv)
    if ScaleByNtQ == 1
      w_t = ones(N_t)
    elseif ScaleByNtQ == 2
      w_t = wM[t,vv]
    else
      w_t = ones(N_t)
    end
    if N_t > 0                       #don't accumulate [] to omega0DK
      (r_t,x_t) = (r_t[vv],x_t[vv,:])
      lmul!(Diagonal(r_t .* w_t.^2),x_t)  #x_t is henceforth the moment condition
      h_t       = sum(x_t,dims=1)/Nb[t]     #for (i,t), x_t .* (r_t .* w_t.^2)
      omega0DK .= omega0DK + h_t'h_t
      omega0W  .= omega0W  + x_t'x_t/Nb[t]^2
      for j = 1:m
        omegajDK[:,:,j] .= omegajDK[:,:,j] + h_t'h_tLag[j:j,:]    #h(t)*h(t-j)'
      end
      h_tLag = [h_t;h_tLag[1:end-1,:]]  #update only if !isempty(rx_t), effectively disregarding t if no data
    end
    if yhatQ
      yhat[t,:] = yhat_t
    end
    msg = IterationPrintPs(t,T,msg,100)
  end

  Shatj = copy(omega0DK)
  for j = 1:m
    Shatj .= Shatj + (1-j/(m+1))*(omegajDK[:,:,j]+omegajDK[:,:,j]')
  end

  zx_1   = inv(xx)
  CovDK  = zx_1 * (omega0DK/Tb^2) * zx_1'           #covariance matrix, DK
  CovW   = zx_1 * (omega0W/Tb^2)  * zx_1'           #covariance matrix, White's
  CovDKj = zx_1 * (Shatj/Tb^2)    * zx_1'           #covariance matrix, DK with lags
  stdDK  = sqrt.(diag(CovDK))                       #standard errors, DK
  stdW   = sqrt.(diag(CovW))                        #standard errors, White's
  stdDKj = sqrt.(diag(CovDKj))                      #standard errors, DK with lags

  if yhatQ
    vv  = FindNoNaNPs(1,vec(yhat),vec(y))
    R2a = cor(vec(yhat)[vv],vec(y)[vv])^2
  else
    R2a  = Float64[]
  end

  #fnOutput = [theta,stdDK,stdW,CovDK,yhat,R2a,CovW,CovDKj,stdDKj]
  fnOutput = (theta=theta,stdDK=stdDK,stdW=stdW,CovDK=CovDK,yhat=yhat,R2a=R2a,
              CovW=CovW,CovDKj=CovDKj,stdDKj=stdDKj)

  return fnOutput

end
#------------------------------------------------------------------------------
