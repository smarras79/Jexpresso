#!/usr/bin/env julia
# Jexpresso-mini-coupled.jl — Julia side with dynamic remote-leader discovery.
# FIXED: Tag mismatch corrected
#
#  mpirun --tag-output -np 2 -x APPID=0 ./alya/alya_mini_coupler : -np 2 -x APPID=1 julia --project=. Jexpresso-mini-coupled.jl
#
#

using MPI

println("Jexpresso-mini starting..."); flush(stdout)
MPI.Init()
println("init"); flush(stdout)

world = MPI.COMM_WORLD
wrank = MPI.Comm_rank(world)
wsize = MPI.Comm_size(world)
println("[Jexpresso rank $wrank] World size: $wsize"); flush(stdout)

# Read APPID (0 or 1) from environment
appid = try parse(Int, get(ENV, "APPID", "2")) catch; 2 end
println("[appid: $appid"); flush(stdout)
if appid < 0
    if wrank == 0
        println("[Jexpresso] ERROR: APPID not set. Launch with -x APPID=0 (Fortran) and -x APPID=1 (Julia).")
    end
    MPI.Abort(world, 1)
end

# Split WORLD into per-app local comms
println("[Split before $wrank"); flush(stdout)
local_comm = MPI.Comm_split(world, appid, wrank)
println("[Split after $wrank"); flush(stdout)
lrank = MPI.Comm_rank(local_comm)
lsize = MPI.Comm_size(local_comm)

# --- Dynamic remote leader discovery via Allgather on WORLD ---

#local_chars = Vector{UInt8}(repeat("popo", 32))  # 128 bytes
local_chars = Vector{UInt8}(rpad("popo", 128, ' '))
recv_buffer = nothing
MPI.Gather!(local_chars, recv_buffer, 0, world)

@info "CALL JEXPRESSO"

#push!(empty!(ARGS), "CompEuler", "wave1d")
#include("./src/Jexpresso.jl")

for i=1:5
    println("hola ", i)
end

MPI.Finalize()
