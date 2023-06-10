"""
    IterationPrintPs(i,N,msgOld,freq=1,msgExtra="")

Prints iteration progress, overwrites old printing

# Input
- `i::Integer`:       iteration number
- `N::Integer`:       total number of iterations
- `msgOld::String`:   old string, exactly this space is overwritten
- `freq::Integer`:    (optional) indicating when to print
- `msgExtra::String`: (optional) extra string to print

# Output
- `msg::String`:msg, s tring printed to screen. Use as msgOld in the next iteration

# Uses
- Printf

Paul.Soderlind@unisg.ch

"""
function IterationPrintPs(i,N,msgOld,freq=1,msgExtra="")
  if mod(i,freq) == 0                      #print only if i/freq is an integer
    reverseStr = "\b" ^ length(msgOld)
    msg        = string(@sprintf("Processed %d/%d",i,N),msgExtra)
    msg        = rpad(msg,60," ")                         #fill at least 60 columns
    @printf("%s",reverseStr * msg)
  else
    msg = msgOld
  end
  return msg
end
#------------------------------------------------------------------------------
