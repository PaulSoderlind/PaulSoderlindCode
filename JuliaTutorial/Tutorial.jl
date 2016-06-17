#------------------ a hash sign makes the rest of line a comment -------------
#
#
# (0)  This file should work well with Julia 0.4. Open up Julia.exe.
#
# (1)  At the prompt, you may change directory by typing
#      cd("directory_name") (e.g. cd("e:/mystuff") and hit Return.
#
# (2)  Don't know in which directory you are? Type pwd() and hit Return.
#
# (3)  To see the files in the current directory, type readdir() and hit Return.
#
# (4)  Change this file in a text editor.
#
# (5)  To execute this program from the command line: type include("Tutorial.jl")
#      and hit Return. If you want to stop execution on line 35 (say), then
#      add error() on that line.
#
# (6)  Leave Julia: type exit() and hit Return.
#
# (7)  You find more information about Julia on, for instance,
#      https://en.wikibooks.org/wiki/Introducing_Julia
#
# (8)  This example file uses several (free) Julia packages. You find instructions
#      in the code below (close to "using..."). For documentation of these packages,
#      see
#      http://distributionsjl.readthedocs.org/en/latest/
#      https://github.com/JuliaLang/Roots.jl
#      https://github.com/JuliaOpt/Optim.jl
#      http://gadflyjl.org/
#
# (9)  Installing/compiling packages takes a bit of time. The next run of
#      the program is much faster. If you get warning messages while
#      installing packages, restart Julia and hope for the best.
#
#
#  Paul Söderlind (Paul.Soderlind at unisg.ch), October 2015, revised April 2016
#-------------------------------------------------------------------------


println("\n","-------------------------- matrices --------------------------","\n")
                                  #println(x) prints x, where x is a matrix or
                                  #a string (text) is like "dogs are better than cats"

Q = 1                             #create a scalar
q = [ 1 2 3;                      #create 2x3 matrix
      4 5 6 ]
println("\n","Q is a scalar. To print, use println()")    #"\n creates a line break"
println(Q)                        #print matrix
println("q is a matrix")
println(q)                        #case sensitive (q and Q are different)

println("\n","first line of q: ",q[1:1,:])
println("second column of q (printed compactly, with commas to indicate rows): ",q[:,2:2])

z = q'                            #transposing
println("\n","z is the transpose of q: ")
println(z)

p = collect(1:10:21)'
println("\n","p is a row vector with a sequence starting at 1 ending in 21: ",p)

z = [q;p]                         #overwriting old z with new
println("\n","stacking q and p vertically")
println(z)

z2 = [q q/10]                     #stacking q and q/10 horizontally
println("\n","stacking q and q/10 horizontally")
println(z2)

w = [ 10 11 12;                   #create another 2x3 matrix
      13 14 15 ]
println("\n","q (again, so you do not forget it)")
println(q)
println("w, another matrix")
println(w)

y1 = q + w                        #matrix addition
println("\n","q + w")
println(y1)

y2 = q.*w                         #element-by-element multiplication
println("\n","q.*w")
println(y2)

y3 = q'*w                          #matrix multiplication (q'w is the same as q'*w)
println("\n","q'*w")
println(y3)


println("\n","----------------------  Load data from ascii file ------------","\n")

#The following is a portion of MyData.csv (# signs added here):
#date,Mkt-RF,RF,SmallGrowth
#197901,4.18,0.77,10.96
#197902,-3.41,0.73,-2.09
#197903,5.75,0.81,11.71
#197904,0.05,0.8,3.27

xx   = readdlm("MyData.csv",',',header=true)
x    = xx[1]                       #xx[1] contains the data
println("Column headers: ",xx[2])  #xx[2] contains the headers
println("first four lines of x:")
println(x[1:4,:])

ym  = x[:,1]                    #yearmonth, like 200712
Rme = x[:,2]                    #picking out second column
Rf  = x[:,3]
R   = x[:,4] +                  #commands continue on the next line
      0.0
Re  = R - Rf


println("\n","--------------- OLS, random numbers and some stats -----------","\n")

ett = ones(size(Rme,1),1)       #a vector with ones, no. rows from variable
x   = [ett Rme]                 #x is Tx2 matrix
y   = Re                        #just to get standard OLS notation
b   = inv(x'x)*x'y              #OLS
b2  = x\y                       #also OLS, much quicker and numerically more stable
u   = y - x*b                   #OLS residuals
R2a = 1 - var(u)/var(y)         #avoid using name R2 (aleady in StatsBase package)
println("OLS coefficients, regressing Re on constant and Rme, different calculations")
println(round([b b2],3))        #round(,3) rounds to 3 digits. Easier to look at
println("R2: ",round(R2a,3))
println("no. of observations: ",size(Re,1))

x = randn(100,3)                  #matrix of random draws from N(0,1)
println("\n","mean and std of random draws from N(0,1): ")
mu    = mean(x,1)                 #mean of each column in matrix, gives row vector
sigma = std(x,1)
println(round([mu;sigma],3))

println("\n","cov(x): ")          #there are lots of statistics functions
println(round(cov(x),3))          #for more, see the package StatsBase.jl

using Distributions  #the first time, do Pkg.add("Distributions") to install the package
println("\n","5th percentile of N(0,1) and 95th of Chisquare(5)")      #lots of statistics functions
println(round([quantile(Normal(0,1),0.05) quantile(Chisq(5),0.95)],3))


println("\n","-------------------- Comparing things --------------------------","\n")

x = collect(linspace(-1.5,0.5,5))       #vector: -1.5,-1.0,-0.5,0,0.5
println("x values: ",x)
vv = -1 .< x .<= 0                      #true for x values (-1,0], boolean
println("x is in (-1,0]: ",vv)
x2 = x[vv]
println("x values in (-1,0]: ",x2)


println("\n","-------------------- Writing a loop --------------------------","\n")

i = 1
x = 0
while i <= 10
  x = x + i
  println("iteration and x value: ",[i x])
  i = i + 1
end


println("\n","-------------------- Writing nested loops --------------------","\n")
println("Pick out elements on and below diagonal from a square matrix x")

x = reshape(1:9,3,3)              #generate square matrix
(m,n) = size(x)                   #find no. rows and columns
println("\n","original matrix: ")
println(x)

nRows = round(Int,n*(n+1)/2)      #must be an integer, so use round()
v = fill(-999,(nRows,1))          #to put results in, initialized as -999
k = 1                             #fill(value,integer,integer) so must convert
for j = 1:n                       #loop over columns in x
  for i = j:n                     #loop over rows on and below diagonal
    v[k] = x[i,j]
    k = k + 1                     #update index in v
  end
end

println("vech of the matrix (stack elements on and below diagonal): ")
println(v)


println("\n","---------------------- Functions -----------------------------","\n")

function MathLecD(a,b)
  value = (a*b-1.1).^2 - 0.5
  return value
end

x = [1 1.5]
y = MathLecD(x,1)                   #calling on the function
println("result from the function MatLecD: ",round(y,4))


println("\n","----------------------------- Plotting -----------------------","\n")

t = collect(-3:6/99:6)

using Gadfly      #the first time, do Pkg.add("Gadfly") to install the package

set_default_plot_size(20cm, 13cm)       #size of plots

plot1 = plot(layer(x=t,y=MathLecD(t,1),Geom.line,Theme(default_color=colorant"red",line_width=2px)),
             layer(x=t,y=MathLecD(t,0.5),Geom.line,Theme(default_color=colorant"blue")),
             Guide.title("Results MathLecD"),
             Guide.xlabel("t"),
             Guide.ylabel("my output value"),
             Guide.manual_color_key(" ",
                          ["From MathLecD(t,1)","MathLecD(t,0.5)"],
                          ["red", "blue"]))
display(plot1)                               #show plot in browser
draw(SVG("Fig1.svg",15cm,10.5cm),plot1)      #save svg, use inkscape to covert (if needed)


println("\n","---------------- Solving (non-linear) equations --------------","\n")

using Roots        #the first time, do Pkg.add("Roots") to install the package
x1 = fzero(x->MathLecD(x,1),[-1;1])     #notice that x->MathLecD(x,1)
                                        #defines a new function (of x only)
                                        #[-1;1] searches roots in this interval
println("at which x is MathLecD(x,1) = 0? ",round(x1,3))

x1 = fzero(x->MathLecD(x,1),[1;3])  #now, try between 1 and 3
println("at which x is MathLecD(x,1) = 0? ",round(x1,3))
println("Yes, there are several roots. Just look at it (in plot).")


println("\n","------------------------- Optimization -----------------------","\n")

using Optim      #the first time, do Pkg.add("Optim") to install the package
Sol = optimize(x->MathLecD(x,1),-2.0,3.0)
println("argmin MathLecD(x,1), optim finds it. ")
println(Optim.minimizer(Sol))


#println("\n","------------------------- load xls file -------------------------------","\n")

#can be done using Taro (which needs the Java) or ExcelReaders (needs Python's xlrd),
#so it is for advanced users


println("\n","------------------------- More plots -------------------------","\n")

YearFrac = floor(ym/100) + (mod(ym,100)-1)/12    #year + (month-1)/12
plot3a = plot(x=YearFrac,y=Rme,Geom.line,Theme(default_color=colorant"blue"),
             Scale.x_continuous(minvalue=1978,maxvalue=2012),
             Guide.xticks(ticks=[1980,1990,2000,2010]),
             Guide.title("Time series plot: monthly US equity market excess return"),
             Guide.xlabel(""),
             Guide.ylabel("%"))
display(plot3a)

plot3b = plot(layer(x=[-40,60],y=[-40,60],Geom.line,Theme(default_color=colorant"black",line_width=1px)),
              layer(x=Rme,y=Re,Geom.point,Theme(default_color=colorant"blue")),
              Scale.x_continuous(minvalue=-40,maxvalue=60),
              Scale.y_continuous(minvalue=-40,maxvalue=60),
              Guide.title("Scatter plot: two monthly return series (and 45 degree line)"),
              Guide.xlabel("Market excess return, %"),
              Guide.ylabel("Excess returns on small growth stocks, %"))
display(plot3b)

plot3c = plot(x=Rme,Geom.histogram(bincount=25),
              Guide.title("Histogram: monthly US equity market excess return"),
              Guide.xlabel("Market excess return, %"),
              Guide.ylabel("Number of months"))
display(plot3c)

println("\n","---------------------- end of program ------------------------","\n")
