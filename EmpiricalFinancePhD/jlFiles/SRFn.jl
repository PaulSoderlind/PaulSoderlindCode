function SRFn(par)
# SRFn    Calculates the Sharpe ratio from the mean and 2nd moment

  mu   = par[1]          #E(Z)
  mu_2 = par[2]          #E(Z^2)
  s2   = mu_2 - mu.^2

  SR   = mu./sqrt(s2)

  return SR

end