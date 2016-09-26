function AssetPricingTest3bPs(Re1,Re2,f,h,bandwidth=0,q=0,PointEstOnlyQ=false)
#AssetPricingTest3bPs  Testing difference of two linear factor models by testing the
#                      difference of the intercepts, or some other combination of
#                      parameters
#
#
#
#
#  Usage:  fnOutput = AssetPricingTest3bPs(Re1,Re2,f,h[,bandwidth[,q[,PointEstOnlyQ]]])
#
#  Input:   Re1             Txn1, test assets
#           Re2             Txn2, test assets     (in many applications, Re2 = Re1)
#           f               TxK, 1st set of factors
#           h               TxL, 2nd set of factors
#           bandwidth       (optional) scalar of Newey-West band width,  [0]
#           q               (optional) heterogenous array of vectors,
#                           expected values of R*theta: Any[alpha,delta], [0]
#           PointEstOnlyQ   (optional), bool: if true, then only the point estimates are calculated [false]
#
#  Output:  fnOutput      8x1 heterogeneous array with the following content:
#             [1] Joint         3x1, in each row WaldStat for alpha, delta, alpha-delta (if n1=n2)
#             [2] Indiv         3xmax(n1,n2), t-stat of alfa,delta and alpha-delta (if n1=n2)
#             [3] bf            n1x(1+K), [alpha beta] in Re1 = alpha + beta*f + errors
#             [4] bh            n2x(1+L), [delta gamma] in Re2 = delta+ gamma*h + errors
#             [5] EReFit        1x2 cell array, each cell is a nx2 matrix with [actual,fitted ERe]
#             [6] theta         [n1(1+K)+n2(1+L)]x1, estimated parameters, [alpha;vec(beta);delta;vec(gamma)]
#             [7] V             matrix, Cov[sqrt(T)*theta]
#
#
#
#  Notice: (a) the order of parameters in theta are  [alpha;vec(beta);delta;vec(gamma)] from
#          Re1 = alpha + beta*f + errors
#          Re2 = delta + gamma*h + errors,
#          where Re1 is an n1x1 vector, Re2 is n2x1, f is is Kx1 and h is Lx1
#
#          (b) theta2 = reshape(theta,n,1+K+1+L)
#          creates a matrix with [alpha beta delta gamma]
#
#
#
#
#  Uses:    excise4mPs, HDirProdPs, NewEst3Ps
#
#
#
#
#  Paul.Soderlind@unisg.ch, April 2015, To Julia Oct 2015
#------------------------------------------------------------------------------

  n1 = size(Re1,2)
  n2 = size(Re2,2)
  if q == 0                                    #scalar 0
    q = Any[zeros(n1,1),zeros(n2,1)]           #default values for options
  end
  #----------------------------------------------------------

  K  = size(f,2)
  L  = size(h,2)

  (Re1,Re2,f,h,) = excise4mPs(Re1,Re2,f,h)
  T              = size(Re1,1)                      #no. obs

  Ef   = mean(f,1)'
  Eh   = mean(h,1)'
  ERe1 = mean(Re1,1)'
  ERe2 = mean(Re2,1)'

  f1 = [ones(T,1) f]
  K1 = 1 + K
  h1 = [ones(T,1) h]
  L1 = 1 + L
                                         #bf = f1\Re1; epsf = Re1 - bf*f1
  bf         = f1\Re1
  epsf       = Re1 - f1*bf
  bf         = bf'                       #n1xK1
  alfaM      = bf[:,1]
  betaM      = bf[:,2:end]               #n1xK, Re1(t) = alfa * beta*f(t) + e(t)
  gf         = HDirProdPs(f1,epsf)       #moment conditions
  ERe1Hatf   = betaM*Ef                  #n1xK1 * K1x1

  bh         = h1\Re2
  epsh       = Re2 - h1*bh
  bh         = bh'                       #n2xL1
  deltaM     = bh[:,1]
  gammaM     = bh[:,2:end]               #n2xL
  gh         = HDirProdPs(h1,epsh)
  ERe2Hath   = gammaM*Eh                 #n2xL1 * L1x1

  ERe1Fit = [ERe1 ERe1Hatf]              #actual and fitted ERe
  ERe2Fit = [ERe2 ERe2Hath]
  EReFit  = Any[ERe1Fit,ERe2Fit]

  theta = [vec(bf);vec(bh)]              #parameters (combined)

  if PointEstOnlyQ
    fnOutput = Any[NaN,NaN,bf,bh,EReFit,theta,NaN,NaN]
    return fnOutput
  end
  #----------------------------------------------------------

  g     = [gf gh]                        #moment conditions (combined)
  S0    = NewEst3Ps(g,bandwidth)         #Cov[sqrt(T)*gbar]
  g  = nothing                           #can save some memory
  gf = nothing
  gh = nothing

  Sff  = f1'f1/T
  Df_1 = -kron(inv(Sff),eye(n1))         #upper left corner of inverse of Jacobian
  Shh  = h1'h1/T
  Dh_1 = -kron(inv(Shh),eye(n2))         #lower right corner of inverse of Jacobian
  D_1  = cat([1,2],Df_1,Dh_1)            #inverse of Jacobian

  V    = D_1*S0*D_1'                     #Cov[sqrt(T)*theta]
  #----------------------------------------------------------

  RM = Any[[eye(n1)      zeros(n1,n1*K) zeros(n1,n2) zeros(n1,n2*L)],        #alpha
           [zeros(n2,n1) zeros(n2,n1*K) eye(n2)      zeros(n2,n2*L)]]        #delta

  if n1 == n2                     #alpha - delta, if n1=n2
    RM = Any[RM[1],RM[2],RM[1] - RM[2]]
    q  = Any[ q[1], q[2], q[1] -  q[2]]
  end

  Indiv = fill(NaN,(length(RM),max(n1,n2)))     #t-stat
  Joint = fill(NaN,length(RM))                  #WaldStat
  for i = 1:length(RM)                          #loop over (alpha,delta,alpha-delta)
    R          = RM[i]
    q_i        = q[i]
    Rstat      = R*theta - q_i                  #eg. alpha - delta - expected value
    Lambda     = R*V*R'                         #eg. Cov[sqrt(T)*(alpha-delta)]
    stat_std   = sqrt(diag(Lambda/T))
    tstat      = Rstat./stat_std
    Joint[i]   = (T * Rstat'inv(Lambda)*Rstat)[1]     #[1] to make it a scalar
    Indiv[i,1:length(tstat)] = tstat'
  end
  #----------------------------------------------------------

             #    1      2    3 4    5      6   7
  fnOutput = Any[Joint,Indiv,bf,bh,EReFit,theta,V]

  return fnOutput

end
#------------------------------------------------------------------------------
