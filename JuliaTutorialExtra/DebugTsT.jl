using Debug

#Instructions for how to run the debugger:
#s for stepping to next line
#c to next bp
#l to see 3 lines above and below the current one
#o to step out of currrent space (eg, a loop)
#q to quit
#
#At the debug propmt, you can print a variable by typing its name
#
#There are more and better examples at https://github.com/toivoh/Debug.jl

@debug function rosenbrock(x)
  x1 = x[1]
  x2 = x[2]
@bp
	LL = (1.0 - x1)^2 + 100.0 * (x2 - x1^2)^2
  return LL                      #before moving on, type LL to see the value
end
z = rosenbrock([0.0, 0.0])
#------------------------------------------------------------------------------


function rosenbrock2(x)
    x1 = x[1]
    x2 = x[2]
    LL = (1.0 - x1)^2 + 100.0 * (x2 - x1^2)^2
    return LL
end

#can also do this to debug code in the main file
@debug begin                                          #create begin ...  end
  x = collect(-3:3)                                   #add @debug and @bp
@bp
  y = cos(x)
  z = rosenbrock2([0.0, 0.0])    #before moving on, type y to see the value
end
#------------------------------------------------------------------------------
