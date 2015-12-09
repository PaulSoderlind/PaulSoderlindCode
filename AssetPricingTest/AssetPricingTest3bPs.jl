function AssetPricingTest3bPs(Re1,Re2,f,h,bandwidth=0,q=0,ExtraTest=[],ExtraTestq=[])
#AssetPricingTest3bPs  Testing difference of two linear factor models by testing the
#                      difference of the intercepts, or some other combination of
#                      parameters
#
#
#
#
#  Usage:  fnOutput = AssetPricingTest3bPs(Re1,Re2,f,h[,bandwidth[,ExtraTest[,ExtraTestq]]]])
#
#  Input:   Re1           Txn1, test assets
#           Re2           Txn2, test assets     (in many applications, Re2 = Re1)
#           f             TxK, 1st set of factors 
#           h             TxL, 2nd set of factors 
#           bandwidth     (optional) scalar of Newey-West band width,  [0]
#           options       (optional) structure, struct('q',0,'ExtraTest', [],'ExtraTestq',0)
#             q             cell array of vectors, expected values of R*theta: {alpha,delta}
#             ExtraTest     cell array of vectors, indices of coeffs in theta 
#                           to test whether = 0 or = ExtraTestq, eg. {{[1,4,5:9];[2:3]}} 
#                           to make 2 tests
#             ExtraTestq    cell array of vectors, matching ExtraTest, expected values of
#                           extra tests, eg {{zeros(7,1);ones(2,1}}
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
#             [8] ExtraWaldStat Wald test stats of the coefficients in ExtraTest
#
#
#
# 
#  Notice:  ExtraTest = {[1,4,5:9]}; ExtraTestq = {[1,2,10,11,12,13,14]}
#           tests the joint hypothesis that theta([1,4,5:9]) = [1,2,10,11,12,13,14] 
#
#
#  Uses:    excise, HDirProdPs, NewEst3Ps
#
#
#
#
#  Paul.Soderlind@unisg.ch, April 2015, To Julia Oct 2015
#------------------------------------------------------------------------------ 
  
  n1         = size(Re1,2)   
  n2         = size(Re2,2) 
  if q == 0                                    #scalar 0   
    q = Any[zeros(n1,1),zeros(n2,1)]           #default values for options
  end
  #----------------------------------------------------------
  
  K  = size(f,2)            
  L  = size(h,2)
  
  yx  = excise([Re1 Re2 f h])          #cut rows with some NaNs
  Re1 = yx[:,1:n1]
  Re2 = yx[:,n1+1:n1+n2]
  f   = yx[:,n1+n2+1:n1+n2+K]            #factors, set f
  h   = yx[:,n1+n2+K+1:n1+n2+K+L]        #factors, set h
  T   = size(Re1,1)                      #no. obs
  yx = nothing                              #can save some memory
  
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
  #(bf,epsf,) = OlsPs(Re1,f1)             #first set of factors
  bf         = bf'                        #n1xK1
  alfaM      = bf[:,1]
  betaM      = bf[:,2:end]                #n1xK, Re1(t) = alfa * beta*f(t) + e(t)
  gf         = HDirProdPs(f1,epsf)        #moment conditions
  ERe1Hatf   = betaM*Ef                   #n1xK1 * K1x1
  
  bh         = h1\Re2
  epsh       = Re2 - h1*bh
  #(bh,epsh,) = OlsPs(Re2,h1)            #second set of factors
  bh         = bh'                       #n2xL1
  deltaM     = bh[:,1]
  gammaM     = bh[:,2:end]               #n2xL
  gh         = HDirProdPs(h1,epsh)
  ERe2Hath   = gammaM*Eh                 #n2xL1 * L1x1
  
  ERe1Fit = [ERe1 ERe1Hatf]              #actual and fitted ERe
  ERe2Fit = [ERe2 ERe2Hath]
  EReFit  = Any[ERe1Fit,ERe2Fit]       
  #----------------------------------------------------------
  
  theta = [vec(bf);vec(bh)]              #parameters (combined)
  g     = [gf gh]                        #moment conditions (combined)
  S0    = NewEst3Ps(g,bandwidth)         #Cov[sqrt(T)*gbar]
  g  = nothing                               #can save some memory
  gf = nothing
  gh = nothing
  
  Sff  = f1'f1/T
  Df_1 = -kron(inv(Sff),eye(n1))          #upper left corner of inverse of Jacobian
  
  Shh  = h1'h1/T
  Dh_1 = -kron(inv(Shh),eye(n2))          #lower right corner of inverse of Jacobian
  
  D_1  = [ Df_1                zeros(n1*K1,n2*L1);  #inverse of Jacobian of g wrt theta'
           zeros(n2*L1,n1*K1)  Dh_1              ]
  
  V    = D_1*S0*D_1'                      #Cov[sqrt(T)*theta] 
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
    Rstat      = R*theta - vec(q_i)             #eg. alpha - delta - expected value
    Lambda     = R*V*R'                         #eg. Cov[sqrt(T)*(alpha-delta)] 
    stat_std   = sqrt(diag(Lambda/T))
    tstat      = Rstat./stat_std
    Joint_i    =  T * Rstat'inv(Lambda)*Rstat
    Joint[i]   = Joint_i[1]                     #[1] to make it a scalar
    Indiv[i,1:length(tstat)] = tstat'
  end  
  #----------------------------------------------------------
  
  ExtraWaldStat = fill(NaN,length(ExtraTest))
  for i = 1:length(ExtraTest)
    vvi = ExtraTest[i]                 #which variables to test whether = 0 or ExtraTestq
    R   = eye(n1*K1+n2*L1)
    R   = R[vvi,:]
    ExtraTestq_i     = ExtraTestq[i]
    Rstat            = R*theta - vec(ExtraTestq_i)     #eg. gamma - 0.25
    Lambda           = R*V*R'                          #eg. Cov[sqrt(T)*gamma] 
    ExtraWaldStat_i  = T * Rstat'inv(Lambda)*Rstat
    ExtraWaldStat[i] = ExtraWaldStat_i[1]
  end
  #---------------------------------------------------------
  
  fnOutput = Any[Joint,Indiv,bf,bh,EReFit,theta,V,ExtraWaldStat] 
  
  return fnOutput

end
#------------------------------------------------------------------------------
