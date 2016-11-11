Readme for JuliaTutorial
========================

This folder contains my Julia tutorial (aimed at students in finance and economics).

1. Tutorial.jl: the main file. Instructions are found at the top of that file. It calls on the following files:
    * printmat.jl which is a function for printing arrays (matrices) in a pretty way
    * MyData.csv which contains data used in the main file

2. Tutorial.ipynb is a jupyter notebook version of Tutorial.jl. You can click it to see it online.
   To use it with your local Julia installation do the following: (a) start Julia; 
   (b) cd(file location) ;(c) using IJulia; (d) notebook(dir=pwd())

3. Tutorial_Finance.jl: some typical finance calculations
    * MyFunctions.jl contains some functions used in Tutorial_Finance.jl

4. Tutorial_Finance.ipynb is a jupyter notebook version of Tutorial_Finance.jl

5. Tutorial_TrickyStuff.ipynb: illustrates some tricky things with Julia (at least for an old Matlab user)

6. Tutorial_MorePlots.jl: contains more advanced plot examples (using another plotting package, PyPlot)
