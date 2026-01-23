#!/usr/bin/env julia
# julia_couple_with_alya_gather.jl
# Couples with Fortran code that uses MPI_Gather on MPI_COMM_WORLD
# Both codes share the same MPI_COMM_WORLD and participate in collective operations

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

# Participate in MPI_Gather
# All ranks send to rank 0 of MPI_COMM_WORLD
# Rank 0 receives from all ranks
if rank == 0
    # Root process: receive from all ranks (including itself)
    recv_buf = Matrix{Int8}(undef, 128, size)
    MPI.Gather!(send_buf, recv_buf, 0, world)

    # Print received names
    println("\n[Julia rank 0] Received application names from all ranks:")
    for i in 1:size
        name_chars = recv_buf[:, i]
        name_str = String(collect(Char.(Int.(name_chars)))) |> x->rstrip(x)
        println("  Rank $(i-1): $name_str")
        flush(stdout)
    end
else
    # Non-root processes: just send
    MPI.Gather!(send_buf, nothing, 0, world)
end

MPI.Barrier(world)
println("[Julia rank $rank] Done")
flush(stdout)

MPI.Finalize()
