function excisePs(x)
#excisePs    Remove rows with NaNs
#
#
#
#
#  Usage:    (x,vvNaN,vvNaNRow,vvNoNaNRow) = excisePs(x)
#
#  Input:   x           Txm   matrix
#
#  Output:  x           Sxm   matrix, S<=T, rows with any NaN are cut out
#           vvNaN       Txm   matrix, true if element (i,j) in x is a NaN
#           vvNaNRow    Tx1   vector, indices of rows that have some NaNs
#           vvNoNaNRow  Tx1   vector, indices of rows that have no NaNs
#
#
#
#
#
#  Paul.Soderlind@unisg.ch   1 March 2005, to Julia Oct 2015
#------------------------------------------------------------------------------

  vvNaN      = isnan(x)             #1 if (i,j) is NaN
  vvRow      = any(vvNaN,2)
  vvNaNRow   = find(vvRow)          #rows with some NaNs
  vvNoNaNRow = find(!vvRow)         #rows without any NaNs
  z          = x[vvNoNaNRow,:]      #only keep rows with no NaNs. clumsy but works

  return z,vvNaN,vvNaNRow,vvNoNaNRow

end
#------------------------------------------------------------------------------
