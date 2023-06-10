function SimpRulT(A,B,Q,R,U,bet,n1,n2,F,SigmaXX,x10,cutoff=1.0)
#SimpRulT   Solving RE model for given simple rule u(t)  = -F * x(t),
#           where   x(t)  = [x1(t)' x2(t)']'.
#
#
#
#  Usage:    (M,C,J0) = SimpRulT(A,B,Q,R,U,bet,n1,n2,F,SigmaXX,x10,cutoff)
#
#  Input:     A,B,Q,R,U,bet,n1,n2: see ComItAlg
#             F          kxn matrix, decision rule
#             SigmaXX    n1xn1 covariance matrix of shocks
#             x10        n1 vector, initial values for predetermined
#                        variables. Not n1x1 array.
#             cutoff     scalar, max modulus of backward looking variables
#
#  Output:    M          n1xn1 matrix, x1(t+1) = M*x1(t) + e(t+1)
#             C          n2xn1 matrix, x2(t) = C*x1(t)
#             J0         scalar, value of loss function
#
#
#
#  Paul Soderlind, Paul.Soderlind@unisg.ch, Aug 2000, to Julia Jan 2016
#-----------------------------------------------------------------------

  Base.require_one_based_indexing(A,B,Q,R,U,F,SigmaXX,x10)

  Q = (Q + Q')/2                  #to make symmetric
  R = (R + R')/2

  A = A - B*F                     #subst for u(t) in evolution
  n = n1 + n2

  GS = schur(Matrix(1.0I,n,n),A)
  logconA = abs.(GS.beta) .<= (abs.(GS.alpha)*cutoff) #selecting stable eigenvalues
  ordschur!(GS,logconA)                                 #reordering, stable first

  if sum(logconA) < n1
    @warn("Too few stable roots: no stable solution")
    return NaN,NaN,NaN
  elseif sum(logconA) > n1
    @warn("Too many stable roots: inifite number of stable solutions")
    return NaN,NaN,NaN
  end

  Stt = GS.S[1:n1,1:n1]
  Zkt = GS.Z[1:n1,1:n1]
  Zlt = GS.Z[n1+1:n,1:n1]
  Ttt = GS.T[1:n1,1:n1]

  if cond(Zkt) > 1e+14
    @warn("Zkt is singular: rank condition for solution not satisfied")
    return NaN,NaN,NaN
  end

  Zkt_1 = inv(Zkt)         #inverting
  Stt_1 = inv(Stt)

  M = real(Zkt*Stt_1*Ttt*Zkt_1)       #x1(t+1) = M*x1(t) + e(t+1)
  C = real(Zlt*Zkt_1)                 #x2(t) = C*x1(t)

  P = [ I;
        C;
        (-F*[I;C]) ]
  QUUR = [ Q   U;
           U'  R ]
  PQUURP = P'QUUR*P

  vecV = (I - kron(M',bet*M'))\vec(PQUURP)
  V    = reshape(vecV,n1,n1)                 #solves V = W + bet*M'V*M

  J0 =  x10'V*x10 + (bet/(1-bet)) * tr(V*SigmaXX)

  return M,C,J0

end
#----------------------------------------------------------------------
