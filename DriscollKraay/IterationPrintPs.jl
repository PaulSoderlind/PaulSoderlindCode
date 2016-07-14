function IterationPrintPs(i,N,msgOld,freq=1,msgExtra="")
#IterationPrintPs    Prints iteration progress, overwrites old printing.
#
#
#
#
#  Usage:  IterationPrintPs(i,N,msgOld,freq,msgExtra)
#
#  Input: i         scalar integer (1,2,3,...), iteration number
#         N         scalar, integer, total number of iterations
#         msgOld    string, old string, exactly this space is overwritten
#         freq      (optional) scalar, integer indicating when to print
#         msgExtra  (optional) string, extra string to print
#
#  Output:         msg, string printed to screen. Used as msgOld in the next iteration
#
#
#
#
#  Paul.Soderlind@unisg.ch   Sep 2012, to Julia Feb 2016
#------------------------------------------------------------------------------

  if mod(i,freq) == 0                      #print only if i/freq is an integer
    reverseStr = "\b" ^ length(msgOld)
    msg        = string(@sprintf("Processed %d/%d",i,N),msgExtra)
    msg        = rpad(msg,60," ")                         #fill at least 60 columns
    @printf("%s",reverseStr * msg)
  else
    msg = msgOld
  end

  return msg
  #fprintf("#s \n"," ");
end
#------------------------------------------------------------------------------

