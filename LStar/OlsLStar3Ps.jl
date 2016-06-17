function OlsLStar3Ps(y,x0,w,ExciseIt,z,gM,cM,gcKeep=[],xwzHat=[])
#OlsLStar3Ps  LSTAR LS of y on (x0,w), y = (1-G)*b1*x0 + G*b2*x0 + d*w,
#             G = 1./(1+exp(-g*(z-c)))
#
#
#
#
#  Usage:    fnOutput = OlsLStar3Ps(y,x0,w,ExciseIt,z,gM,cM[,gcKeep[,xwzHat]])
#
#
#  Input:    y            Tx1, dependent variable
#            x0           Txk, regressors that have regime shifts (including deterministic ones)
#            w            Txkw, regressors that do not have regime shifts, can be empty, Array{Float64}(T,0)
#            ExciseIt     bool, true: excise([y,x0,w,z])
#            z            Tx1, argument of G(z) function
#            gM           Ngx1, different values of g to try in a loop
#            cM           Ncx1, different values of c to try in a loop
#            gcKeep       (optional) 2x1, [g,c], if an element is a NaN, then this parameter is estimated
#            xwzHat       (optional) Tz x(k+kw+1) x nPred, [x0,w,z] values at which to calculate predicted values
#
#  Output:   fnOutput     heterogeneous (Any[]) array with
#              [1]  sse          scalar, loss fn value at point estimate
#              [2]  theta        (2+2k+kw)x1, (1+2k+kw)x1 or (0+2k+kw)x1, parameter estimates:
#                                [g;c;b;d],[c;b;d],[g;b;d] or [b;d],
#                                b is a (k+k)x1 vector for x0, the first k elements are for z=-Inf
#                                and the second k elements for z=Inf
#                                d are coefficiets for w (no regimes)
#              [3]  Stdtheta       ""              , standard errots of theta
#              [4]  Covtheta     cov(theta)
#              [5]  b            coeffs from traditional LS (conditional on g and c)
#              [6]  Stdb_ols     standard errors according to traditional LS (conditional on g and c)
#              [7]  R2a          scalar, coefficient of determination
#              [8]  Gquant       19x2, [linspace(0.01,0.99,99)',quantiles of G(z)]
#              [9]  gc           1x2, [g,c] prespecified/estimated
#              [10] sseM         NgxNc, sum of squared errors for different values of g and c
#              [11] yHat         Tzx1, predicted values at xwzHat
#              [12] slopeDiff    matrix, [b tstat]
#              [13] yHatLH       Tz x 3 x nPred matrix of pedicted values, [yHat from b1,yHat from b2,G]
#
#  Calls on: excise, OlsPs, NewEst3Ps, NumJac3Ps
#
#
#  Notice:  (a) z is NOT standardized inside this function
#           (b) the optimization for Octave does not report an exit flag, so the results need to
#               checked
#
#
#
#  Paul.Soderlind@unisg.ch                Jan 2013, to Julia Nov 2015
#----------------------------------------------------------------------------

  Ng = length(gM)
  Nc = length(cM)

  k  = size(x0,2)
  kw = size(w,2)

  if ExciseIt
    (y,x0,w,z) = excise4mPs(y,x0,w,z)
  end

  sseM = fill(NaN,(Ng,Nc))             #calculate sse in loop over g and c values
  for i = 1:Ng
    for j = 1:Nc
      sse_ij    = OlsLStar3LossPs([gM[i] cM[j]],y,x0,w,z,[gM[i] cM[j]])
      sseM[i,j] = sse_ij
    end    #j
  end   #i
  (minLoss,vvMin) = findmin(sseM)            #minimum loss, for which i
  (i,j) = ind2sub(size(sseM),vvMin)
  (g,c,b,par0,) = OlsLStar3Par(gcKeep,[NaN;NaN;NaN],[gM[i];cM[j]])

  if !isempty(par0)
    Sol = optimize(par->OlsLStar3LossPs(par,y,x0,w,z,gcKeep),par0,x_tol=1e-6)
    parX = Optim.minimizer(Sol)
    if !Optim.converged(Sol)
      warn("no convergence")
      return
    end
  else
    parX = Float64[]
  end

  #fnOutput = Any[sse,theta,Stdtheta,Covtheta,b,Stdb_ols,R2a,Gquant,gc]
  fnOutput  = OlsLStar3LossAllPs(parX,y,x0,w,z,gcKeep,1)
  theta     = fnOutput[2]
  Covtheta  = fnOutput[4]

  bDiff      = fill(NaN,k)                             #calculate and test b2-b1=0
  tstatbDiff = fill(NaN,k)
  vv    = length(theta) - kw - 2*k + (1:2*k)   #pick out slopes from theta and their Cov matrix
  b2    = theta[vv]
  Covb2 = Covtheta[vv,vv]
  for j = 1:k
    R             = zeros(length(b2),1)
    vv            = [j j+k]      #location of b1 and b2 for same x0(:,j)
    R[vv]         = [-1;1]
    bDiff[j]      = (R'b2)[1]
    tstatbDiff[j] = (R'b2/sqrt(R'Covb2*R))[1]
  end
  slopeDiff = [bDiff tstatbDiff]

  if isempty(xwzHat) || all(isnan(xwzHat))
    yHat   = Float64[]
    yHatLH = Float64[]
  else                                  #predicted values
    nPred = size(xwzHat,3)
    yHatLH = fill(NaN,(size(xwzHat,1),3,nPred))
    for i = 1:nPred
      (yHat,yHat2,G) = OlsLStar3PredPs(xwzHat[:,:,i],k,kw,theta,gcKeep)
      yHatLH[:,:,i]  = [(yHat - yHat2) yHat2 G]
    end
    yHat = OlsLStar3PredPs(xwzHat[:,:,1],k,kw,theta,gcKeep)
    if nPred == 1
      yHatLH = reshape(yHatLH,size(xwzHat,1),3)     #better than squeeze(yHatLH,3)
    end
  end

  #fnOutput = Any[sse,theta,Stdtheta,Covtheta,b,Stdb_ols,R2a,Gquant,gc]
  push!(fnOutput,sseM,yHat,slopeDiff,yHatLH)
  return fnOutput

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function OlsLStar3PredPs(xwzHat,k,kw,theta,gcKeep)

  (g,c,b,) = OlsLStar3Par(gcKeep,theta,[NaN;NaN])

  x0 = xwzHat[:,1:k]             #[k,kw,1]
  w  = xwzHat[:,1+k:k+kw]
  z  = xwzHat[:,end]

  G  = 1./(1+exp(-g*(z-c)))
  x1 = x0.*repmat(1-G,1,k)
  x2 = x0.*repmat(G,1,k)
  x  = [x1 x2 w]

  yHat = x*b

  yHat2 = [zeros(size(x1)) x2 zeros(size(w))]*b        #constribution of x2 only

  return yHat,yHat2,G

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function OlsLStar3LossPs(par,y,x0,w,z,gcKeep=Float64[])        #just sse from OlsLStar3LossAllPs

  fnOutput = OlsLStar3LossAllPs(par,y,x0,w,z,gcKeep,0)
  sse = 1.0 + fnOutput[1]

  return sse

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function OlsLStar3LossAllPs(par,y,x0,w,z,gcKeep=Float64[],DetailsIt=0)

  (g,c,) = OlsLStar3Par(gcKeep,par,[NaN;NaN])
  (T,k) = size(x0)

  G  = 1./(1+exp(-g*(z-c)))
  x1 = x0.*repmat(1-G,1,k)
  x2 = x0.*repmat(G,1,k)
  x  = [x1 x2 w]

  (b,res,yhat,Covb,R2a,) = OlsPs(y,x)
  Stdb_ols              = sqrt(diag(Covb))
  sse                   = sum(res.^2)
  theta                 = [vec(par);b]

  if DetailsIt == 1
    (mbar,m) = OlsLStar3MomCondAllPs(theta,y,x0,w,z,gcKeep)
    S0       = NewEst3Ps(m,0)                                    #ACov(sqrt(T)*mbar)
    D0       = NumJac3Ps(x->OlsLStar3MomCondPs(x,y,x0,w,z,gcKeep),theta,Float64[],3)  #gradient of mbar
    Covtheta = inv(D0)*S0*inv(D0)'/T                        #Cov(theta)
    Stdtheta = sqrt(diag(Covtheta))
    Gquant   = collect(linspace(0.01,0.99,99))
    Gquant   = [Gquant quantile(vec(excise(G)),Gquant)]
    gc       = [g c]
  else
    Covtheta = Float64[]
    Stdtheta = Float64[]
    Gquant   = Float64[]
    gc       = Float64[]
  end

  fnOutput = Any[sse,theta,Stdtheta,Covtheta,b,Stdb_ols,R2a,Gquant,gc]
  return fnOutput

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

function OlsLStar3MomCondPs(theta,y,x0,w,z,gcKeep)   #just mbar from  OlsLStar3MomCondAllPs

  mbar, = OlsLStar3MomCondAllPs(theta,y,x0,w,z,gcKeep)

  return mbar

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function OlsLStar3MomCondAllPs(theta,y,x0,w,z,gcKeep)   #moment conditions

  k  = size(x0,2)
  kw = size(w,2)

  #(res,G,g,c,b,x) = OlsLStar3RegFuncPs(theta,y,x0,w,z,gcKeep)
  (g,c,b,_,EstType) = OlsLStar3Par(gcKeep,theta,[NaN;NaN])

  G   = 1./(1+exp(-g*(z-c)))
  x1  = x0.*repmat(1-G,1,k)
  x2  = x0.*repmat(G,1,k)
  x   = [x1 x2 w]
  res = y - x*b

  b1 = b[1:k]
  b2 = b[k+1:2*k]

  b1_b2x0 = x0*(b2-b1)
  dF_dg   = (1-G).*G.*(z-c).*b1_b2x0
  dF_dc   = (1-G).*G.*(-g).*b1_b2x0

  mg = res.*dF_dg
  mc = res.*dF_dc
  mb = repmat(res,1,2*k+kw).*x

  #(xx_,c_,b_,par0_,EstType) = OlsLStar3Par(gcKeep,theta,[NaN NaN])
  if EstType == 1
    m = -[mg mc mb]
  elseif EstType == 2
    m = -[   mc mb]
  elseif EstType == 3
    m = -[mg    mb]
  elseif EstType == 4
    m = -[      mb]
  else
    error("invalid case")
  end

  mbar = mean(m,1)'

  return mbar,m

end
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
function OlsLStar3Par(gcKeep,theta,gcM)

  if isempty(gcKeep) || all(isnan(gcKeep))       #gcKeep=[] or [NaN,NaN], estimate both g and c
    g = abs(theta[1])
    c = theta[2]
    b = theta[3:end]
    EstType = 1
    par0    = [gcM[1];gcM[2]]
  elseif !isnan(gcKeep[1]) && isnan(gcKeep[2])   #gcKeep=[1.5,NaN], don't estimate g, estimate c
    g = gcKeep[1]
    c = theta[1]
    b = theta[2:end]
    EstType = 2
    par0    = gcM[2]
  elseif isnan(gcKeep[1]) && !isnan(gcKeep[2])   #gcKeep=[NaN,-0.75], estimate g, don't estimate c
    g = theta[1]
    c = gcKeep[2]
    b = theta[2:end]
    EstType = 3
    par0    = gcM[1]
  elseif all(!isnan(gcKeep))                    #gcKeep=[1.5,-0.75], don't estimate g or c
    g = gcKeep[1]
    c = gcKeep[2]
    b = theta
    EstType = 4
    par0    = Float64[]
  else
    error("invalid case")
  end

  return g,c,b,par0,EstType

end
#------------------------------------------------------------------------------
