function DiscAlg(A,B,Q,R,U,bet,n1,n2,Vt1,Ct1,
                ConvCrit,Vweight,Fweight,CritLags,step,PrintIt,MaxIter)
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
#             MaxIter    scalar, maximum number of iterations (eg. 10000). NEW (March 2003)
#
#  Output:    M        n1x1 matrix, x1(t+1) = M*x1(t) + e(t+1)
#             C        n2xn1 matrix, x2(t)  = C*x1(t)
#             V        n1xn1 matrix, value function is x1(t)'*V*x1(t)
#             F        kxn1 matrix, decision rule is u(t) = -F*x1(t), where
#                      k is number of elements in u(t)
#
#
#  Remark:    (For Octave users) If PrintIt is 1 or 2, then set
#             page_output_immediately = 1
#             This seems to make Octave work better.
#
#  Calls on:  DiscAlg2
#
#
#  Paul Söderlind, Paul.Soderlind@unisg.ch, Aug 2000, Mar 2003
#-----------------------------------------------------------------------

  Q = (Q + Q')/2                #to make symmetric
  R = (R + R')/2

  n        = n1 + n2
  Ft1      = 1000
  ConvCrit = collect(ConvCrit)'             #to row vector

  M = []
  C = []
  F = []
  V = []
  Cdiff = 1000*ones(1+CritLags,2)
  iter  = 1
  while any(maximum(Cdiff,1) .> ConvCrit) & (iter < MaxIter)   #iterations

    (M,C,F,V) = DiscAlg2(A,B,Q,R,U,bet,n1,n2,Ct1,Vt1)  #solve period t

    Vdiff = maximum(Vweight.*abs(V-Vt1))        #changes t+1 -> t
    Fdiff = maximum(Fweight.*abs(F-Ft1))
    Cdiff = [ Cdiff;
              Vdiff Fdiff ]
    Cdiff = Cdiff[2:size(Cdiff,1),:]                   #latest is last

    Vt1 = step*V + (1-step)*Vt1                        #"downdating"
    Ct1 = C + 0.0
    Ft1 = step*F + (1-step)*Ft1

    if PrintIt == 1
      println(iter)
    elseif PrintIt == 2
      println([iter maximum(Cdiff)])
    end

    iter = iter + 1

  end                                 #end iterations

  if iter >= MaxIter
    warn("Maximum number of iterations reached")
  end

  return M,C,V,F

end
#-----------------------------------------------------------------------

