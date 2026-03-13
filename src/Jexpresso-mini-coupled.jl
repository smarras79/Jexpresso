#!/usr/bin/env julia
# Jexpresso-mini-coupled.jl — Julia side with dynamic remote-leader discovery.
#
# Launch coupled Jexpresso + Alya via:
#   mpirun --tag-output -np 2 ./alya/Alya_enhanced.x : -np 2 julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya
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

lrank        = MPI.Comm_rank(local_comm)
lsize        = MPI.Comm_size(local_comm)
nranks1      = lsize                      # Jexpresso
nranks2      = wsize - lsize              # Other code
local_chars  = Vector{UInt8}(rpad("JEXPRESSO", 128, ' '))
recv_buffer  = nothing
is_jexpresso = (appid == 2)

MPI.Gather!(local_chars, recv_buffer, 0, world)

# Set coupling mode to prevent auto-execution on module load
ENV["JEXPRESSO_COUPLING_MODE"] = "false"

# Set command line arguments for Jexpresso
push!(empty!(ARGS), "CompEuler", "thetaAlya")

# Load Jexpresso module (setup doesn't run yet because of JEXPRESSO_COUPLING_MODE)
println("[Jexpresso rank $wrank] Loading Jexpresso module (JIT compilation may take minutes)..."); flush(stdout)
include("./Jexpresso.jl")
println("[Jexpresso rank $wrank] Jexpresso module loaded."); flush(stdout)

# Set the custom local communicator (for Jexpresso-internal MPI: mesh partitioning, etc.)
Jexpresso.set_mpi_comm(local_comm)
println("[Jexpresso rank $wrank] MPI communicator set (local_comm, size=$lsize)."); flush(stdout)

# Store the world communicator for inter-code coupling with Alya
Jexpresso.set_mpi_comm_world(world)
println("[Jexpresso rank $wrank] World communicator stored for coupling."); flush(stdout)

#--------------------------------------------------------------------------------------------
# Receive ndime from Alya via Bcast on COMM_WORLD (all ranks must participate)
# Alya: call MPI_Bcast(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
# Use Int32 to match Fortran's MPI_INTEGER (4 bytes)
#--------------------------------------------------------------------------------------------
ndime_buf = Vector{Int32}(undef, 1)
MPI.Bcast!(ndime_buf, 0, world)
ndime = ndime_buf[1]
println("[Jexpresso rank $wrank] Received ndime = $ndime from Alya"); flush(stdout)

rem_min = zeros(Float64, 3)
rem_max = zeros(Float64, 3)
rem_nx  = zeros(Int32,   3)
MPI.Bcast!(rem_min, 0, world)
MPI.Bcast!(rem_max, 0, world)
MPI.Bcast!(rem_nx,  0, world)
println("[Jexpresso rank $wrank] Received ndime      = $ndime      from Alya"); flush(stdout)
println("[Jexpresso rank $wrank] Received rem_min    = $rem_min    from Alya"); flush(stdout)
println("[Jexpresso rank $wrank] Received rem_max    = $rem_max    from Alya"); flush(stdout)
println("[Jexpresso rank $wrank] Received rem_nx     = $rem_nx     from Alya"); flush(stdout)

alya2world_l = zeros(Int32, nranks2)
alya2world   = MPI.Allreduce(alya2world_l, MPI.SUM, world)

a_l = zeros(Int32, wsize, wsize)
a   = MPI.Allreduce(a_l, MPI.SUM, world)

println("[Jexpresso rank $wrank] Received alya2world = $alya2world from Alya"); flush(stdout)
#--------------------------------------------------------------------------------------------
# END Receive ndime from Alya
#--------------------------------------------------------------------------------------------

# Store coupling data so that je_couplingSetup() inside the module can access it
Jexpresso.set_coupling_data(Dict{Symbol,Any}(
    :ndime      => ndime,
    :rem_min    => rem_min,
    :rem_max    => rem_max,
    :rem_nx     => rem_nx,
    :alya2world => alya2world,
    :couple_matrix => a,
))
println("[Jexpresso rank $wrank] Coupling data stored."); flush(stdout)

# Now run Jexpresso with the configured communicator
Jexpresso.jexpresso_main()

MPI.Finalize()
