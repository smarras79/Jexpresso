#!/usr/bin/env julia
# julia_couple_with_alya_gather.jl
# Couples with Fortran code that uses MPI_Gather on MPI_COMM_WORLD
# In MPMD mode: Fortran gets ranks 0-N, Julia gets ranks N+1 onwards
# Only Fortran rank 0 receives the gathered data

using MPI

println("Julia coupling code starting...")
flush(stdout)

MPI.Init()

world = MPI.COMM_WORLD
rank = MPI.Comm_rank(world)
size = MPI.Comm_size(world)

println("[Julia rank $rank] World size: $size")
flush(stdout)

# Prepare application name (128 characters, space-padded) as Int8 (MPI_CHARACTER compatible)
app_name_str = "JEXPRESSO"
send_buf = fill(Int8(' '), 128)
nb = min(length(codeunits(app_name_str)), 128)
@inbounds for i in 1:nb
    send_buf[i] = Int8(codeunits(app_name_str)[i])
end

# Julia ranks are NEVER rank 0 in MPMD mode (Fortran gets rank 0)
# So Julia always participates as non-root sender
println("[Julia rank $rank] Participating in MPI_Gather to root (rank 0)")
flush(stdout)

MPI.Gather!(send_buf, nothing, 0, world)

println("[Julia rank $rank] MPI_Gather completed")
flush(stdout)

MPI.Finalize()
println("[Julia rank $rank] Done")
flush(stdout)
