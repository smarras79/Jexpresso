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
#local_chars = Vector{UInt8}(lpad("JEXPRESSO", 128, ' '))
#recv_buffer = (wrank == 0) ? Vector{UInt8}(undef, 128 * wsize) : nothing
#MPI.Gather!(local_chars, recv_buffer, 0, world)
#if wrank == 0
#    # parse recv_buffer into chunks of 128 bytes per rank
#end

local_chars = Vector{UInt8}(lpad("JEXPRESSO", 128, ' '))
recv_buffer = nothing
MPI.Gather!(local_chars, recv_buffer, 0, world)

#local_chars  = Vector{UInt8}(lpad("JEXPRESSO", 128, ' '))
#all_chars    = Vector{UInt8}(undef, 128 * wsize)
#MPI.Allgather!(local_chars, all_chars, world)
@info "CALL JEXPRESSO"

# Receive ndime from Alya via Bcast on COMM_WORLD (all ranks must participate)
# Alya: call MPI_Bcast(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
# Use Int32 to match Fortran's MPI_INTEGER (4 bytes)
ndime_buf = Vector{Int32}(undef, 1)
MPI.Bcast!(ndime_buf, 0, world)
ndime = ndime_buf[1]
println("[Jexpresso rank $wrank] Received ndime = $ndime from Alya"); flush(stdout)

# Receive per-dimension arrays: rem_min (REAL), rem_max (REAL), rem_nx (INTEGER)
# Alya broadcasts one element at a time inside a do loop over idime=1,ndime
# Use Float32 for Fortran MPI_REAL (4 bytes) and Int32 for MPI_INTEGER (4 bytes)
rem_min = Vector{Float32}(undef, ndime)
rem_max = Vector{Float32}(undef, ndime)
rem_nx  = Vector{Int32}(undef, ndime)
for idime in 1:ndime
    MPI.Bcast!(@view(rem_min[idime:idime]), 0, world)
    MPI.Bcast!(@view(rem_max[idime:idime]), 0, world)
    MPI.Bcast!(@view(rem_nx[idime:idime]),  0, world)
end
println("[Jexpresso rank $wrank] Received rem_min=$rem_min, rem_max=$rem_max, rem_nx=$rem_nx"); flush(stdout)

# Set coupling mode to prevent auto-execution on module load
ENV["JEXPRESSO_COUPLING_MODE"] = "true"

# Set command line arguments for Jexpresso
push!(empty!(ARGS), "CompEuler", "wave1d")

# Load Jexpresso module (setup doesn't run yet because of JEXPRESSO_COUPLING_MODE)
println("[Jexpresso rank $wrank] Loading Jexpresso module (JIT compilation may take minutes)..."); flush(stdout)
include("./src/Jexpresso.jl")
println("[Jexpresso rank $wrank] Jexpresso module loaded."); flush(stdout)

# Set the custom local communicator
Jexpresso.set_mpi_comm(local_comm)
println("[Jexpresso rank $wrank] MPI communicator set (local_comm, size=$lsize)."); flush(stdout)

# Now run Jexpresso with the configured communicator
println("[Jexpresso rank $wrank] Starting jexpresso_main()..."); flush(stdout)
Jexpresso.jexpresso_main()
println("[Jexpresso rank $wrank] jexpresso_main() finished."); flush(stdout)



MPI.Finalize()
