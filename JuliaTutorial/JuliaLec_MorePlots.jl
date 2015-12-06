#  This file demonstrates how to create plots in Julia by using the PyPlot package
#
#
#  PyPlot which relies on the matplotlib library, which is part of Python. If you have
#  Python installed, then it will be used as is. Otherwise, see PyPlot's homepage
#  https://github.com/stevengj/PyPlot.jl
#  for instruction on how to install what is needed. 
#
#  Notice: Restart Julia before running this file (at least if you have used 
#          another plotting package)
#
#  Paul Söderlind (Paul.Soderlind at unisg.ch), December 2015
#------------------------------------------------------------------------------

using PyPlot       #the first time, do Pkg.add("PyPlot") to install the package
#------------------------------------------------------------------------------

b = linspace(-3,3,20)                 #create some "data" to plot
c = linspace(1,7,25)' 

loss1 = 2*b.^2 + 0.5
loss2 = fill(NaN,(length(b),length(c)))  #to put results in, initialized as NaNs
for j = 1:length(c)                      #create loss2 column by column
  loss2[:,j] = 2*b.^2 + (c[j]-4)^2 - 0.0*b.*(c[j]-4) 
end
#------------------------------------------------------------------------------

t = collect(-3:6/99:6)           #A FIRST PLOT

figure(figsize=(12,8.5))
  plot(b,loss1,linestyle="--",color="b",linewidth=1.0)
  plot(b,log(loss1),linestyle="-",color="r",linewidth=3.0)
  title("a title")
  xlim(-3,3)               # set limits of the x-axis
  ylim(-1,20)              # set limits of the y-axis
  xlabel("b")
  ylabel("loss values")
  text(-2.5,0.9,"some text")
  legend(["Loss1";"log Loss1"])
  savefig("MathLec_MorePlots1.pdf")      #save pdf file of the plot
#------------------------------------------------------------------------------


figure(figsize=(12,8.5))        #BAR, STEP, SURFACE        
  bar(b,loss1,facecolor="red",edgecolor="white",align="center",width=0.3) 
  title("Bar chart")

figure(figsize=(12,8.5))
  step(b,loss1,linewidth=1.5)
  title("Stairs plot")

figure(figsize=(12,8.5))
  surf(c,b,loss2,rstride=1,cstride=1)    #rstride and cstride improves the look
  xlim(1,7)
  ylim(-3,3)
  zlim(0,30)
  title("Surface plot")
#------------------------------------------------------------------------------

N = 51
x = randn(N,1)                         #SCATTER, HISTOGRAM 
y = rand(N,1) 
areas = rand(51)*100                   #size of the scatter points

figure(figsize=(10,10))
  scatter(x,y,s=areas,alpha=0.5)
  title("Scatter plot")
  xlabel("x")
  ylabel("y")

figure(figsize=(10,10))
  plt[:hist](x,21)                 
  grid("on")
  title("Histogram")
  xlabel("x")
#------------------------------------------------------------------------------

                                     #TIME SERIES PLOT 
FirstDate = Date(2015,12,4)          #just faking some dates
dN = Array(Date,length(x))
for i = 1:length(dN)
  dN[i] = FirstDate + Dates.Day(i-1)       #add a day
end  

figure(figsize=(35,35/1.41))
  plot_date(dN,cumsum(x),"k-")
  title("Cumulative x",fontsize=18)
#------------------------------------------------------------------------------

using LaTeXStrings                 #add some LaTeX to the figure

t = collect(-3:6/99:6)          

PyPlot.matplotlib[:rc]("font", family="serif",size=10)  #font similar to LaTeX
figure(figsize=(12,8.5))
  plot(b,loss1)
  title("a title")
  xlabel("b")
  ylabel("log loss")
  text(-2.5,0.9,L"some text, $\ln(\mathrm{loss})$")
  text(-2.5,5,L"$\mu_2 = \int x^2 f(x) dx$")
#------------------------------------------------------------------------------
