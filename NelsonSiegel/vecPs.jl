function vecPs(x)

  if length(x) == 1          #if scalar or 1x1 array -> 1x1 array
    y = vec(collect(x))
  else
    y = deepcopy(vec(x))
  end

  return y

end