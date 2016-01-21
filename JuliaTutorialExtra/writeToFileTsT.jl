fh = open("newfile.txt", "w")  # "w" for writing

write(fh, "testing\n")
write(fh, "more testing\n")

x = rand(1:7,(4,3))
println(fh,x)

close(fh)