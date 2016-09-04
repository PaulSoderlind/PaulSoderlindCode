# Introduction
The .ipynb files are Julia notebooks. They use some data, which you find in the Data folder. The notebooks are meant to be used together with my lecture notes (pdf files). Actually, many of the examples in the lecture notes come from the calculations in the notebooks. 


# Instructions

You use IJulia to edit and run these Julia notebooks. You can either install Julia/IJulia on your local computer (see https://sites.google.com/site/paulsoderlindecon/home/software for instructions) or you can use JuliaBox.com. The rest of these instructions are for the JuliaBox.com alternative.

1. Log in at JuliaBox.com
2. In the menu, go to Files. Upload the notebook for Chapter 1.
3. In the menu, go to IJulia and open up the notebook by clicking on it. You can now start working. Using packages (Gadfly, Roots, Optim) may take some time on the first run. If those packages do not work, follow the instructions below.


## Installing more Packages

To install a package, do as follows:
In the menu, go to Console and type in 
>julia 

and hit enter. At the julia prompt type (to install Gadfly)
>Pkg.add("Gadfly") 

and hit enter. 
Then type 

>using Gadfly

and hit enter.
For the notebooks you also need the packages "Roots" and "Optim".
You can now leave julia (by typing exit() and hitting enter).

## Working with Files

The IJulia tab allows you to create folders and delete files. For instance, create a subfolder called Data and then go to Files, click on the Data folder, and then upload the csv files to this subfolder.
