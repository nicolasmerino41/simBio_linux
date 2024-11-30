using CSV, DataFrames, Dates
epsilon
file = CSV.File("gbif_sizes.csv") |> DataFrame
println(epsilon)
println(file[1, 1])
current_datetime = now()
println("You started at: 16:40h and the finish date and time is: ", current_datetime)