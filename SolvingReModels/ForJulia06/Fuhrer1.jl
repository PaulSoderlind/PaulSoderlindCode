function Fuhrer1(a1,ap,Varey,w,gamm,Varep,Dbig,qy,qpi,qf,bet)
#Fuhrer1 Generate A,B,K,Q,U, and R from parameter values of simplified Fuhrer model.
#
#
#  Purpose:   Generate A,B,K,Q,U, and R from parameter values of simplified
#             Fuhrer model.
#
#  Usage:     (A,B,K,Q,U,R) =
#              Fuhrer1(a1,ap,Varey,w,gamm,Varep,Dbig,qy,qpi,qf,bet)
#
#  The period loss function is
#
#  qy*y^2 +  qpi*pi^2 + qf*i^2, which is rewritten here as
#
#  x'Qx + 2x'Uu + u'Ru.
#
#  x = [x1,x2], with x1(t) = [ep(t),y(t),Dx(t-1)] and x2(t) = [R[t],Dx(t)]
#
#  See lecture notes for further comments.
#
#
#
#
#
#
#  Paul Soderlind, Paul.Soderlind@hhs.se, 20 August 1997, to MatLab Aug 2000
#----------------------------------------------------------------------------

                             #model->general setup
  Ax = [ (  -1/(w*(1-w)^2)                        );
         (  -gamm*a1/(w*(1-w)) - gamm/(1-w)^2     );
         (  -1                                    );
         (  -gamm*ap/(w*(1-w))                    );
         (   2   ) ]'

  Ar = [ 0  0  0  ((1+Dbig)/Dbig)  (4/Dbig*(1-w) ) ]
  Ar =  Ar + 4/Dbig*w*Ax


  A = [ zeros(1,5);                    #x(t+1)=A*x(t)+B*u(t)+ e(t+)
        0  a1   0   ap  0;
        0  0    0    0  1;
        Ar;
        Ax ]

  B = [ 0;0;0;(-1/Dbig);0 ]

  K  = [ 0  1  0         0   0      0;   #[y,pi,f] = K*[x,u]
         0  0  4*(1-w)   0   4*w    0;
         0  0  0         0   0      1  ]

  Wbig = K'diagm([qy;qpi;qf])*K
  Q = Wbig[1:5,1:5]
  U = Wbig[1:5,6]
  R = Wbig[6,6]

  return A,B,K,Q,U,R

end
#-----------------------------------------------------------------------
