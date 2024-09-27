numericaloutput = 2 * rand()
println("Hello from", numericaloutput)
write(open("output.txt","w"), string("Hello world from Julia #", numericaloutput))

