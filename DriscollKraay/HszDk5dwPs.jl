function HszDk5dwPs(y,x,z,yhatQ=false,m=0,ScaleByNtQ=0,vvzx=[],wM=[])
#HszDk5dPs  LS and Driscoll-Kray standard errors for unbalanced panel, assuming
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
#                          at least contain Tx1 matrix ones(T,1)
#            z             TxNxL matrix with L (time-varying) individual characteristics, should
#                          at least contain TxNx1 array of ones, ones(T,N,1)
#            yhatQ         (optional) scalar, 1: generate and report fitted values
#            m             (optional), scalar, number of lags in covariance estimation
#            ScaleByNtQ    (optional), scalar, 1: scales all moment conditions by N(t), not by N;
#                          2: scale all moment conditions by wM; [0]
#                          can be used to replicate a portfolio (Calendar time) approach
#            vvzx          (optional), vector of indices of colums in HDirProdPs(z_t,x0_t) to keep
#            wM            (optional) TxN matrix, weights on observation (t,i)
#
#  Output:   fnOutput      9x1 heterogeneous array with the following contents:
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
#  Uses:     excisePs and HDirProdPs
#
#
#
#  Paul.Soderlind@unisg.ch   May 2010, to Julia Oct 2015
#------------------------------------------------------------------------------

  if isempty(vvzx)
    vvzx = 1:(size(x,2)*size(z,3))
  end

  (T,N) = size(y,1,2)
  L     = size(z,3)
  #KL    = size(x,2)*L
  KL    = sum(vvzx .>0)


  xx = 0.0                            #Sum[x(t)*x(t)',t=1:T]
  xy = 0.0                            #Sum[x(t)*y(t),t=1:T]
  Nb = Array{Int}(T)                  #effective number of obs, after pruning NaNs
  for t = 1:T                            #loop over time
    y_t   = y[t,:]'                        #dependent variable, Nx1
    x0_t  = repmat(x[t,:],N,1)             #factors, NxK
    z_t   = reshape(z[t,:,:],N,L)          #NxL, better than squeeze (cf 2-d arrays)
    x_t   = HDirProdPs(z_t,x0_t)           #effective regressors, z_t is NxL, x_t is NxK
    x_t   = x_t[:,vvzx]
    (yx_t,_,_,vvNoNaNRow) = excisePs([y_t x_t])   #pruning NaNs
    N_t   = size(yx_t,1)
    if ScaleByNtQ == 1
      Nb[t] = N_t
      w_t   = ones(N_t,1)
    elseif ScaleByNtQ == 2
      Nb[t] = N
      w_t   = wM[t,vvNoNaNRow]'
    else
      Nb[t] = N
      w_t   = ones(N_t,1)
    end
    if !isempty(yx_t)                    #don't accumulate [] to xx and xy (generates [])
      yx_t = yx_t .* repmat(w_t,1,1+KL)  #put weights on observation i (in t)
      y_t  = yx_t[:,1]
      x_t  = yx_t[:,2:end]
      xx   = xx + x_t'x_t/Nb[t]
      xy   = xy + x_t'y_t/Nb[t]
    end
  end

  Tb  = sum(Nb .> 0)                    #number of effective time periods
  xx  = xx/Tb
  xy  = xy/Tb
  theta  = xx\xy                      #ols estimates, solves xx*theta = xy

  if yhatQ
    yhat = fill(NaN,(T,N))
  else
    yhat = Float64[]
  end
  omega0DK = zeros(KL,KL)               #DK, lag 0
  omega0W  = zeros(KL,KL)               #White's
  omegajDK = zeros(KL,KL,m)             #DK, lags 1 to m
  h_tLag   = zeros(m,KL)                #lag1;lag2;...,lagm
  for t = 1:T                            #loop over time
    y_t    = y[t,:]'
    x0_t   = repmat(x[t,:],N,1)
    z_t    = reshape(z[t,:,:],N,L)
    x_t    = HDirProdPs(z_t,x0_t)
    x_t    = x_t[:,vvzx]
    yhat_t = x_t*theta
    r_t    = y_t - yhat_t
    (rx_t,_,_,vvNoNaNRow) = excisePs([r_t x_t])
    N_t = size(rx_t,1)
    if ScaleByNtQ == 1
      w_t = ones(N_t,1)
    elseif ScaleByNtQ == 2
      w_t = wM[t,vvNoNaNRow]'
    else
      w_t = ones(N_t,1)
    end
    if !isempty(rx_t)                    #don't accumulate [] to omega0DK
      rx_t     = rx_t .* repmat(w_t,1,1+KL);  #put weights on observation i (in t)
      r_t      = rx_t[:,1]
      x_t      = rx_t[:,2:end]
      hi_t     = x_t.*repmat(r_t,1,KL)    #moment condition for (i,t)
      h_t      = sum(hi_t,1)/Nb[t]
      omega0DK = omega0DK + h_t'h_t
      omega0W  = omega0W + hi_t'hi_t/Nb[t]^2
      for j = 1:m
        omegajDK[:,:,j] = omegajDK[:,:,j] + h_t'h_tLag[j,:]    #h(t)*h(t-j)'
      end
      h_tLag = [h_t;h_tLag[1:end-1,:]]  #update only if !isempty(rx_t), effectively disregarding t if no data
    end
    if yhatQ
      yhat[t,:] = yhat_t'
    end
  end
  Shat  = omega0DK/Tb^2                   #estimate of S, DK
  Shatw = omega0W/Tb^2                    #estimate of S, White's
  Shatj = omega0DK/Tb^2
  for j = 1:m
    Shatj = Shatj + (1-j/(m+1))*(omegajDK[:,:,j]+omegajDK[:,:,j]')/Tb^2
  end


  zx_1  = inv(xx)
  CovDK = zx_1 * Shat * zx_1'                      #covariance matrix, DK
  stdDK = sqrt( diag(CovDK) )                      #standard errors, DK
  CovW  = zx_1 * Shatw * zx_1'                     #covariance matrix, White's
  stdW  = sqrt( diag(CovW) )                       #standard errors, White's
  CovDKj = zx_1 * Shatj * zx_1'                    #covariance matrix, DK with lags
  stdDKj = sqrt( diag(CovDKj) )                    #standard errors, DK with lags


  if yhatQ
    yy,   = excisePs([vec(yhat) vec(y)])
    R2a  = cor(yy[:,1],yy[:,2])^2
  else
    R2a  = Float64[]
  end

  fnOutput = Any[theta,stdDK,stdW,CovDK,yhat,R2a,CovW,CovDKj,stdDKj]

  return fnOutput

end

#------------------------------------------------------------------------------
