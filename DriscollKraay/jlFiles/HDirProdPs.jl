"""
    HDirProdPs(x,y)

Calculates horizontal direct product of two matrices with equal number of rows.
              z[i,:] is the Kronecker product of x[i,:] and y[i,:]

# Usage
z = HDirProdPs(x,y)

# Input
- `x::Array`:      T x Kx matrix
- `y::Array`:      T x Ky matrix

# Output
- `z::Array`:      T x (Kx*Ky)


# Example
```julia
julia>  x = [ 1 2;
              3 4]
julia>  y = [5 6 1;
             7 8 1]
julia> HDirProdPs(x,y)
2x6 Array{Int64,2}:
  5   6  1  10  12  2
 21  24  3  28  32  4
```

"""
function HDirProdPs(x,y)

  Base.require_one_based_indexing(x,y)

  Kx = size(x,2)       #columns in x
  Ky = size(y,2)       #columns in y

  z = repeat(y,1,Kx) .* kron(x,ones(Int,1,Ky))  #Int: more general, small perf penalty

  return z

end
#----------------------------------------------------------------------------
