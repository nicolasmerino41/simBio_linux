numericaloutput = arg * 2
println("Hello from", arg)
write(open("output.txt","w"),string("Hello world from Julia #",arg))