#!/usr/bin/env julia
# Jexpresso-mini-coupled.jl — Julia side with dynamic remote-leader discovery.
# FIXED: Tag mismatch corrected
#
# mpirun --tag-output -np 2 ./alya/Alya.x : \
#         -np 2 julia --project=. Jexpresso-mini-coupled.jl \
#         false --gather-coupling --coupling-test-only --code-name "Jexpresso"
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

# Set coupling mode to prevent auto-execution on module load
ENV["JEXPRESSO_COUPLING_MODE"] = "true"

# Set command line arguments for Jexpresso
push!(empty!(ARGS), "CompEuler", "wave1d")

#@info " AAAA ",
#@info " lrank: " ,  lrank , " Comm_size: ", MPI.Comm_size(local_comm)

# TO DO:
# PASS local_comm as global var to all functions that use MPI.COMM_WORLD
#

MPI.COMM_WORLD = local_comm
@info MPI.COMM_WORLD

# Load Jexpresso module (setup doesn't run yet because of JEXPRESSO_COUPLING_MODE)
include("./src/Jexpresso.jl")

# Set the custom local communicator
Jexpresso.set_mpi_comm(local_comm)

# Now run Jexpresso with the configured communicator
Jexpresso.jexpresso_main()



MPI.Finalize()
