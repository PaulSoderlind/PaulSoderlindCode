function NormalHistLoss(par,Probs,Bounds)
# ------------------------------------------------------------------------------
#   NormalHistLoss
#
#   proc for creating loss function from observed histogram and fitted
#   probabilities from normal distribution.
#
#
#
#
#
#   n interval boundaries, z(1), z(2), and z(3) for which
#   the survey gives the following probabilities:
#
#   Pr(x < z(1)), Pr(z(1) <= x <z(2)), Pr(z(2) <= x< z(3)), and Pr(z(3) <= x),
#
#   so there are n+1 categories (and therefore probabilities)
#
#   Notice: both the lowest and the highest intervals are open
#
#
#   Paul Soderlind (Paul.Soderlind@hhs.se), 10 February 2000, to Julia Dec 2016
#------------------------------------------------------------------------------

  n = length(Bounds)                        #n bounds, n+1 categories

  mu = par[1]
  s2 = par[2]^2

  xnorm        = (Bounds-mu)./sqrt(s2)                 #(x-mu)/s
  TheoryProb_x = 0.5 + 0.5*erf(xnorm/sqrt(2))          #theoretical Pr( z<=x(i) )

  TheoryProb_Interval = [ TheoryProb_x[1]                        ;
                          TheoryProb_x[2:n] - TheoryProb_x[1:n-1];
                          1.0 - TheoryProb_x[n]                  ]

  loss = 10000*sum( (TheoryProb_Interval-Probs).^2 )  #loss function

  return loss

end
#----------------------------------------------------------------------------

