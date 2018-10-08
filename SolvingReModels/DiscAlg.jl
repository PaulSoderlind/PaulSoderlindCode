function DiscAlg(A,B,Q,R,U,bet,n1,n2,Vt1,Ct1,
          ConvCrit=[1e-1;1e-5],Vweight=0,Fweight=1,CritLags=0,step=1,PrintIt=0,MaxIter=1e+4)
#DiscAlg    Solves the LQ problem under discretion, iterating backwards in time.
#
#
#
#  Usage:     (M,C,V,F) = DiscAlg( A,B,Q,R,U,bet,n1,n2,Vt1,Ct1,...
#                                  ConvCrit,Vweight,Fweight,CritLags,step,PrintIt,MaxIter )
#
#  Input:     A,B,Q,R,U,bet,n1,n2:   - see ComItAlg
#             Vt1        n1xn1 matrix: initial guess of value matrix
#             Ct1        n2xn1 matrix, initial guess of C in x2(t)=C*x1(t)
#             ConvCrit   1x2 convergence criteria for abs(V-Vt1)|abs(F-Ft1)
#             Vweight    scalar or n1xn1 matrix with weights for
#                        difference criterion for V
#             Fweight    scalar or (n1+n2)x1 vector with weights for
#                        difference criterion for F
#             CritLags   no. lags of CritVar compared with ConvCrit
#             step       scalar in (0,1): factor of updating of V and F as
#                        in Vt1 = step*V + (1-step)*Vt1
#             PrintIt    1: printing iteration number
#                        2: printing  iteration number and convergence criteria
#             MaxIter    scalar, maximum number of iterations (eg. 10000).
#
#  Output:    M        n1x1 matrix, x1(t+1) = M*x1(t) + e(t+1)
#             C        n2xn1 matrix, x2(t)  = C*x1(t)
#             V        n1xn1 matrix, value function is x1(t)'*V*x1(t)
#             F        kxn1 matrix, decision rule is u(t) = -F*x1(t), where
#                      k is number of elements in u(t)
#
#
#  Paul Söderlind, Paul.Soderlind@unisg.ch, Aug 2000, to Julia Jan 2016
#-----------------------------------------------------------------------

  Q = (Q + Q')/2                #to make symmetric
  R = (R + R')/2

  n        = n1 + n2
  Ft1      = fill(1000,(1,n1))
  ConvCrit = vec(collect(ConvCrit))'             #to row vector

  (M,C,F,V) = (Float64[],Float64[],Float64[],Float64[])    #to make visible outside loop
  Cdiff = 1000*ones(1+CritLags,2)
  iter  = 1
  while any(maximum(Cdiff,dims=1) .> ConvCrit) & (iter < MaxIter)   #iterations

    (M,C,F,V) = DiscAlg2(A,B,Q,R,U,bet,n1,n2,Ct1,Vt1)  #solve period t

    Vdiff = maximum(Vweight.*abs.(V-Vt1))        #changes t+1 -> t
    Fdiff = maximum(Fweight.*abs.(F-Ft1))
    Cdiff = [ Cdiff[2:end,:];                    #latest is last
              Vdiff Fdiff ]
    Vt1 = step*V + (1-step)*Vt1                        #"downdating"
    Ct1 = copy(C)
    Ft1 = step*F + (1-step)*Ft1

    if PrintIt == 1
      println(iter)
    elseif PrintIt == 2
      println([iter maximum(Cdiff)])
    end

    iter = iter + 1

  end                                 #end iterations

  if iter >= MaxIter
    @warn("Maximum number of iterations reached")
  end

  return M,C,V,F

end
#-----------------------------------------------------------------------


function DiscAlg2(A,B,Q,R,U,bet,n1,n2,Ct1,Vst1)
#DiscAlg2   Solves the LQ problem under discretion in t, given C(t+1) and V(t+1).
#
#
#
#
#  Usage:     [M,C,Fs,Vs] = DiscAlg2( A,B,Q,R,U,bet,n1,n2,Ct1,Vst1 )
#
#  Input:     A,B,Q,R,U,bet,n1,n2: see ComItAlg
#             Ct1            n2xn1  matrix, C(t+1) in x2(t+1) = C(t+1)x1(t+1)
#             Vst1           n1xn1  matrix, V(t+1) value function matrix in t+1
#
#  Output:    M              n1xn1 matrix,
#             C              n2xn1 matrix,
#             Fs             kxn1 matrix,
#             Vs             n1xn1 matrix,
#
#  Paul Söderlind, Paul.Soderlind@unisg.ch, Aug 2000, to Julia Jan 2016
#-----------------------------------------------------------------------

  n = n1 + n2

  v1  = 1:n1
  v2  = n1+1:n

                #partitioning of A,Q,B,U could be moved outside of the loop
  (A11,A12,A21,A22) = (A[v1,v1],A[v1,v2],A[v2,v1],A[v2,v2])
  (Q11,Q12,Q21,Q22) = (Q[v1,v1],Q[v1,v2],Q[v2,v1],Q[v2,v2])
  (B1,B2) = (B[v1,:],B[v2,:])
  (U1,U2) = (U[v1,:],U[v2,:])

  d1mat = inv(A22-Ct1*A12)
  D     = d1mat * (Ct1*A11-A21)
  G     = d1mat * (Ct1*B1-B2)

  As = A11 + A12*D
  Bs = B1  + A12*G

  Qs = Q11 + Q12*D + D'Q21 + D'Q22*D
  Us = Q12*G + D'Q22*G + U1 + D'U2
  Rs = R .+ G'Q22*G + G'U2 + U2'G              #R is scalar, G'Q22*G is 1x1 array

  Fs = inv(Rs + bet*Bs'Vst1*Bs) * (Us' + bet*Bs'Vst1*As)    #u(t) = -F*x1(t)
  Vs = Qs - Us*Fs - Fs'Us' + Fs'Rs*Fs + bet*(As-Bs*Fs)'Vst1*(As-Bs*Fs)

  C = D - G*Fs                   #x2(t)=C*x1(t)
  M = A11 + A12*C - B1*Fs        #x1(t+1) = M*x1(t) + e(t+1)

  return M,C,Fs,Vs

end
#-----------------------------------------------------------------------

