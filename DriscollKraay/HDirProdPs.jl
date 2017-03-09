function HDirProdPs(x,y)
#HDirProdFn    Calculates horizontal direct product of two matrices with equal number of rows.
#              z[i,:] is the Kronecker product of x[i:i,:] and y[i:i,:]
#
#
#  Usage:    z = HDirProdPs(x,y)
#
#  Input:    x      T x Kx matrix
#            y      T x Ky matrix
#
#  Output:   z      T x (Kx*Ky)
#
#
#  Example:   x = [ 1 2;
#                   3 4]
#             y = [5 6 1;
#                  7 8 1]
#
#             then z = HDirProdPs(x,y) gives
#
#             z = [ 5  6  1  10  12  2;
#                  21 24  3  28  32  4]
#
# Paul.Soderlind@unisg.ch, Oct 2015
#----------------------------------------------------------------------------

  Kx = size(x,2)       #columns in x
  Ky = size(y,2)       #columns in y

  z = repmat(y,1,Kx) .* kron(x,ones(Int,1,Ky))  #Int: more general, small perf penalty

  return z

end
#----------------------------------------------------------------------------
