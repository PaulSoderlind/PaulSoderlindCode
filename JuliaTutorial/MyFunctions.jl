function MeanVarNoSSPs(mu,Sigma,mu_p)
#MeanVarNoSSPs     Calculate mean variance frontier when no short sales are allowed
#
#
#
#
#  Usage:  (Var_p,w_p) = MeanVarNoSSPs(mu,Sigma,mu_p)
#
#  Input:     mu        Nx1, vector of mean returns
#             Sigma     NxN, covariance matrix of returns, can contain riskfree assets
#             mu_p      Kx1, mean returns to calculate results for
#
#
#  Output:    Var_p     Kx1, variance of mean-variance portfolio (risky only) with mean mu_p
#             w_p       KxN, portfolio weights of       ""
#
#
#
#  Notice:  (a) uses quadprog from MathProgBase
#           (b) with a riskfree asset, put zeros row j and column j in Sigma
#
#
#
#  Paul.Soderlind@unisg.ch   27 October 2005, to Julia Nov 2015
#------------------------------------------------------------------------------

  K = length(mu_p)                #MV with no short-sales, numerical minimization
  N = length(mu)

  vv = find((mu_p .>= minimum(mu)) & (mu_p .<= maximum(mu)))  #solve only if feasible

  lb   = zeros(N,1)              #w >= 0
  ub   = ones(N,1)               #w <= 1
  Aeq  = [ones(1,N);mu']         #1'w=1, mu'w = mu_p

  w_p   = fill(NaN,(K,N))
  Var_p = fill(NaN,K)                       #to store results in
  for i = vv
    beq  = [1;mu_p[i]]
    Sol = quadprog(zeros(N),Sigma,Aeq,'=',beq,zeros(N),ones(N),IpoptSolver(print_level=1))
    w_i = Sol.sol
    if Sol.status == :Optimal
      w_p[i,:] = w_i'
      Var_p[i] = (w_i'*Sigma*w_i)[1]
    end
  end

  return  Var_p,w_p

end
#------------------------------------------------------------------------------
