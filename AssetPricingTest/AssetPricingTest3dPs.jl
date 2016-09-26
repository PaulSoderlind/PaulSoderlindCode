function AssetPricingTest3dPs(Re1,Re2,f,h,bandwidth=0,q=0,CrossRegType="beta",vvR=0)
#AssetPricingTest3dPs   Testing difference of two linear factor models by testing the
#                       difference of cross-sectional pricing errors. Allows Re1 ~= Re2
#
#
#
#
#  Usage: fnOutput = AssetPricingTest3dPs(Re1,Re2,f,h,[,bandwidth[,q[,CrossRegType[,vvR]]]])
#
#  Input:   Re1            Txn1, test assets
#           Re2            Txn2, test assets     (in many applications, Re2 = Re1)
#           f              TxK, 1st set of factors, not necessarily excess returns
#           h              TxL, 2nd set of factors
#           bandwidth      (optional) scalar of Newey-West band width,  [0]
#           q              (optional) vector, expected values of R*theta: [alpha delta]
#           CrossRegType   (optional) ('beta' or 'GLS') or {Kxn1;Lxn2} cell arrray, ['beta']
#           vvR            (optional) heterogeneous array of vectors, indices of
#                          moment conditions to test with R*g-q
#
#
#  Output:  fnOutput      8x1 heterogeneous array with the following contents:
#             [1] Joint         3x1, in each row WaldStat for alpha, delta, alpha-delta (if n1=n2)
#             [2] Indiv         3xmax(n1,n2), t-stat of alfa,delta and alpha-delta (if n1=n2)
#             [3] bf
#             [4] bh
#             [5] EReFit        1x2 cell array, each cell is a nx2 matrix with [actual,fitted ERe]
#             [6] theta         [n1(1+K)+n2(1+L)]x1, estimated parameters, [alpha;vec(beta);delta;vec(gamma)]
#                              Notice the order of parameters: [for const; for factor 1;...]
#             [7] V             matrix, Cov[sqrt(T)*theta]
#
#
#
#  Uses:    excise, HDirProdPs, NewEst3Ps
#
#
#
#
#  Paul.Soderlind@unisg.ch, June 2015, to Julia Dec 2015
#------------------------------------------------------------------------------

  n1 = size(Re1,2)
  n2 = size(Re2,2)
  if q == 0                                  #scalar 0
    q = Any[zeros(n1,1),zeros(n2,1)]         #default values for options
  end
  if vvR == 0
    vvR = Any[1:n1,1:n2]
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

  bf         = f1\Re1
  epsf       = Re1 - f1*bf
  bf         = bf'                       #n1xK1
  alfaM      = bf[:,1]
  betaM      = bf[:,2:end]               #n1xK, Re1(t) = alfa * beta*f(t) + e(t)
  gf         = HDirProdPs(f1,epsf)       #moment conditions

  bh         = h1\Re2
  epsh       = Re2 - h1*bh
  bh         = bh'                       #n2xL1
  deltaM     = bh[:,1]
  gammaM     = bh[:,2:end]               #n2xL
  gh         = HDirProdPs(h1,epsh)
  #----------------------------------------------------------

  if CrossRegType == "beta"
    B = betaM'                               #traditional cross-sectional regression
    G = gammaM'
  elseif CrossRegType == "GLS"
    B = betaM'* inv(cov(epsf))               #GLS
    G = gammaM'*inv(cov(epsh))
  else
    B = CrossRegType[1]                      #predefined matrix
    G = CrossRegType[2]
  end
                                             #cross-sectional estimate of price of factor risk
  lambda  = (B*betaM)\(B*ERe1)               #solve (B*betaM)*lambda = B*ERe
  psi     = (G*gammaM)\(G*ERe2)
  theta   = [vec(bf);lambda;vec(bh);psi]     #parameters (combined)

  EReHatf = betaM*lambda                     #nxK1 * K1x1
  EReHath = gammaM*psi                       #nxL1 * L1x1
  EReFit  = Any[[ERe1 EReHatf],[ERe2 EReHath]]  #actual and fitted ERe

  gfb     = Re1 .- EReHatf'                  #moment condition ERe - beta*lambda
  ghb     = Re2 .- EReHath'                  #moment condition ERe - gamma*psi

  g       = [gf gfb gh ghb]                  #moment conditions (combined)
  S       = NewEst3Ps(g,bandwidth)           #Cov[sqrt(T)*gbar]
  gf   = nothing                               #can save some memory
  gfb  = nothing
  gh   = nothing
  ghb  = nothing

  gbar = mean(g,1)'                          #moment conditions, sample average
  gNum = size(g,2)                           #no. moment conditions


  Sff  = f1'f1/T                             #Jacobian
  Df11 = kron(Sff,eye(n1))                   #upper left corner
  Shh  = h1'h1/T
  Dh11 = kron(Shh,eye(n2))
  Df12 = kron([0.0 vec(lambda)'],eye(n1))    #lower left
  Dh12 = kron([0.0 vec(psi)'],   eye(n2))
  Df   = - [ Df11  zeros(n1*(1+K),K);
             Df12  betaM           ]
  Dh   = - [ Dh11  zeros(n2*(1+L),L);
             Dh12  gammaM          ]
  D    = cat([1,2],Df,Dh)                    #block diagonal cat

  A    = cat([1,2],eye(n1*K1),B,eye(n2*L1),G)


  AD_1 = inv(A*D)
  V    = AD_1*A*S*A'AD_1'                    #Cov[sqrt(T)*theta]
  Psia = eye(gNum) - D*AD_1*A
  Psi  = Psia*S*Psia'                        #Cov[sqrt(T) gbar], rank gNum-k
  #----------------------------------------------------------

  RM1 = [zeros(n1,n1*K1) eye(n1)      zeros(n1,n2*L1) zeros(n1,n2)]    #Re1 - beta*lambda
  RM2 = [zeros(n2,n1*K1) zeros(n2,n1) zeros(n2,n2*L1) eye(n2)     ]    #Re2 - gamma*psi
  RM  = Any[RM1[vvR[1],:],RM2[vvR[2],:]]                               #pick out rows

  if size(RM[1],1) == size(RM[2],1)
    RM = Any[RM[1],RM[2],RM[1] - RM[2]]
    q  = Any[ q[1], q[2], q[1] -  q[2]]
  end

  Indiv = fill(NaN,(length(RM),max(n1,n2)))     #t-stat
  Joint = fill(NaN,length(RM))                  #WaldStat
  for i = 1:length(RM)                          #loop over (alpha,delta,alpha-delta)
    R          = RM[i]
    q_i        = q[i]
    Rstat      = R*gbar - vec(q_i)              #eg. gbar[2:3]
    Lambda     = R*Psi*R'                       #eg. Cov[sqrt(T)*gbar[2:3]
    stat_std   = sqrt(diag(Lambda/T))
    tstat      = Rstat./stat_std
    Joint_i    =  T * Rstat'pinv(Lambda)*Rstat  #reduced rank
    Joint[i]   = Joint_i[1]                     #[1] to make it a scalar
    Indiv[i,1:length(tstat)] = tstat'
  end
  #----------------------------------------------------------

              #    1     2   3   4   5      6   7
  fnOutput = Any[Joint,Indiv,bf,bh,EReFit,theta,V]

  return fnOutput

end
#------------------------------------------------------------------------------
