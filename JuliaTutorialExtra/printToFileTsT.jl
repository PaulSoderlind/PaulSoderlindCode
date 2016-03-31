fh = open("newfile.txt", "w")  # "w" for writing

println(fh,"testing")
println(fh,"more testing")

x = rand(1:7,(4,3))
println(fh,x)

close(fh)

println("newfile.txt has been created in the current folder")