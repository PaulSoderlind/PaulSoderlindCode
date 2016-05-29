function vecPs(x)

  if length(x) == 1          #if scalar - > 1x1 array
    y = vec(collect(x))
  else
    y = vec(x)
  end

  return y

end