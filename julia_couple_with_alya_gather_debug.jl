#!/usr/bin/env julia
# julia_couple_with_alya_gather_debug.jl
# DEBUG version with extra synchronization points

using MPI

println("Julia coupling code starting...")
flush(stdout)

MPI.Init()

world = MPI.COMM_WORLD
rank = MPI.Comm_rank(world)
size = MPI.Comm_size(world)

println("[Julia rank $rank] World size: $size")
flush(stdout)

# Try an initial barrier to see if Fortran has reached MPI_Init
println("[Julia rank $rank] Trying initial barrier after MPI_Init")
flush(stdout)
try
    MPI.Barrier(world)
    println("[Julia rank $rank] Initial barrier completed")
    flush(stdout)
catch e
    println("[Julia rank $rank] Initial barrier FAILED: $e")
    flush(stdout)
end

# Prepare application name (128 characters, space-padded)
app_name_str = "JEXPRESSO"
send_buf = fill(Int8(' '), 128)
nb = min(length(codeunits(app_name_str)), 128)
@inbounds for i in 1:nb
    send_buf[i] = Int8(codeunits(app_name_str)[i])
end

println("[Julia rank $rank] Prepared send buffer, waiting before gather")
flush(stdout)

# Sleep to give Fortran time to initialize
sleep(2)

println("[Julia rank $rank] Participating in MPI_Gather to root (rank 0)")
flush(stdout)

MPI.Gather!(send_buf, nothing, 0, world)

println("[Julia rank $rank] MPI_Gather completed")
flush(stdout)

# Try another barrier to sync before finalize
println("[Julia rank $rank] Trying barrier before finalize")
flush(stdout)
try
    MPI.Barrier(world)
    println("[Julia rank $rank] Pre-finalize barrier completed")
    flush(stdout)
catch e
    println("[Julia rank $rank] Pre-finalize barrier FAILED: $e")
    flush(stdout)
end

println("[Julia rank $rank] Calling MPI_Finalize")
flush(stdout)

MPI.Finalize()

println("[Julia rank $rank] Done")
flush(stdout)
