Readme for JuliaTutorial
========================

This folder contains my Julia tutorials (aimed at students in finance and economics). 

*  Most files are jupyter notebooks. Click one of them to see it online. Sometimes it renders better if you activate the nbviewer (see the symbol in the upper right corner once you have opened the file).
*  To edit and run it online, use JuliaBox.com.
*  To edit and run it with your local Julia installation do the following: (a) start Julia; (b) cd(file location); (c) using IJulia; (d) notebook(dir=pwd()). You clearly need IJulia to be installed for this.

The files are:

1. Tutorial.ipynb: the main file. Instructions are found at the top of that file. It calls on the following files:
    * printmat.jl which is a function for printing arrays (matrices) in a pretty way
    * MyData.csv which contains data used in the main file

2. Tutorial_Finance.ipynb: some typical finance calculations
    * MyFunctions.jl contains some functions used in Tutorial_Finance.jl

4. Tutorial_TrickyStuff.ipynb: illustrates some tricky things with Julia (at least for an old Matlab user)

5. Tutorial_MorePlots.ipynb: contains more advanced plot examples (using another plotting package, PyPlot).
