"""
    NormalHistLoss
Create loss function from observed histogram and fitted probabilities from normal distribution.

# Notice
n interval boundaries, z(1), z(2), and z(3) for which the survey gives the following probabilities:
Pr(x < z(1)), Pr(z(1) <= x <z(2)), Pr(z(2) <= x< z(3)), and Pr(z(3) <= x),
so there are n+1 categories (and therefore probabilities)

Notice: both the lowest and the highest intervals are open

Paul Soderlind (Paul.Soderlind@hhs.se), 10 February 2000, to Julia Dec 2016
"""
function NormalHistLoss(par,Probs,Bounds)

  n = length(Bounds)                        #n bounds, n+1 categories

  mu = par[1]
  s2 = par[2]^2

  xnorm        = (Bounds.-mu)./sqrt(s2)                 #(x-mu)/s
  TheoryProb_x = 0.5 .+ 0.5*erf.(xnorm/sqrt(2))          #theoretical Pr( z<=x(i) )

  TheoryProb_Interval = [ TheoryProb_x[1]                        ;
                          TheoryProb_x[2:n] - TheoryProb_x[1:n-1];
                          1.0 - TheoryProb_x[n]                  ]

  loss = 1.0 + 10000*sum( (TheoryProb_Interval-Probs).^2 )  #loss function

  return loss

end
#----------------------------------------------------------------------------

"""
    NormPdfPs
Pdf value of a normally distributed variable
"""
function NormPdfPs(x,mu=0,s2=1)
  s = sqrt(s2)
  z    = (x - mu)/s
  pdfx = exp(-0.5*z^2)/(sqrt(2*pi)*s)
  return pdfx
end
#----------------------------------------------------------------------------

"""
  TriangularPdfPs
Pdf value of a triangularly distributed variable over (a,b) with mode c.
"""
function TriangularPdfPs(x,a,b,c)
  if (x<a) || (x>b)     #outside (a,b)
    pdfx = 0
  elseif a <= x < c
    pdfx = 2*(x-a)/((b-a)*(c-a))
  elseif y == c
    pdfx = 2/(b-a)
  else
     pdfx = 2*(b-x)/((b-a)*(b-c))
  end
  return pdfx
end
#------------------------------------------------------------------------------
