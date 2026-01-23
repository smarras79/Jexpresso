#!/usr/bin/env julia
# jexpresso_coupled.jl
# Lightweight entry point for intercommunicator coupling
# Initializes coupling before loading heavy Jexpresso modules

println("[Jexpresso] Starting coupling initialization...")
flush(stdout)

using MPI
using ArgParse

function parse_coupling_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "eqs"
        help = "Equations directory"
        default = "CompEuler"
        required = false

        "eqs_case"
        help = "Case name"
        default = "wave1d"
        required = false

        "CI_MODE"
        help = "CI_MODE: true or false"
        default = "false"
        required = false

        "--code-name"
        help = "Code name for coupling"
        arg_type = String
        default = "Jexpresso"
    end
    return parse_args(s)
end

# Parse args before MPI init
args = parse_coupling_args()
code_name = args["code-name"]

println("[Jexpresso] Initializing MPI...")
flush(stdout)

MPI.Init()

world = MPI.COMM_WORLD
wrank = MPI.Comm_rank(world)
wsize = MPI.Comm_size(world)

if wrank == 0
    println("[Jexpresso] World size: $wsize")
    flush(stdout)
end

# Read APPID from environment
appid = try
    parse(Int, get(ENV, "APPID", "-1"))
catch
    -1
end

if appid < 0
    if wrank == 0
        error("APPID environment variable not set. Launch with -x APPID=1")
    end
    MPI.Abort(world, 1)
end

if wrank == 0
    println("[Jexpresso] APPID=$appid")
    flush(stdout)
end

# Split by APPID to create local communicator
if wrank == 0
    println("[Jexpresso] Creating local communicator...")
    flush(stdout)
end

local_comm = MPI.Comm_split(world, appid, wrank)
lrank = MPI.Comm_rank(local_comm)
lsize = MPI.Comm_size(local_comm)

if wrank == 0
    println("[Jexpresso] Local communicator: rank=$lrank/$lsize")
    flush(stdout)
end

# Discover remote leader
my_appid = Int32(appid)
world_appids = Vector{Int32}(undef, wsize)
MPI.Allgather!(Ref(my_appid), world_appids, world)

remote_leader_world = -1
for i in 1:wsize
    if world_appids[i] != my_appid
        global remote_leader_world = i - 1
        break
    end
end

if remote_leader_world < 0
    if wrank == 0
        error("Could not find remote leader")
    end
    MPI.Abort(world, 2)
end

if wrank == 0
    println("[Jexpresso] Remote leader: rank $remote_leader_world in world")
    flush(stdout)
end

# Create intercommunicator
if wrank == 0
    println("[Jexpresso] Creating intercommunicator...")
    flush(stdout)
end

libmpi = isdefined(MPI, :libmpi) ? MPI.libmpi : MPI.MPI_LIBRARY
newcomm_ref = Ref{MPI.MPI_Comm}()
tag = Int32(12345)

ierr = ccall((:MPI_Intercomm_create, libmpi), Cint,
             (MPI.MPI_Comm, Cint, MPI.MPI_Comm, Cint, Cint, Ref{MPI.MPI_Comm}),
             local_comm.val, Cint(0), world.val, Cint(remote_leader_world), Cint(tag), newcomm_ref)

if ierr != 0
    error("MPI_Intercomm_create failed with error code $ierr")
end

inter_comm = MPI.Comm(newcomm_ref[])

if wrank == 0
    println("[Jexpresso] Intercommunicator created successfully")
    flush(stdout)
end

# Exchange names
if lrank == 0
    send_name = fill(Int8(' '), 128)
    name_bytes = codeunits(code_name)
    nb = min(length(name_bytes), 128)
    @inbounds for i in 1:nb
        send_name[i] = Int8(name_bytes[i])
    end

    recv_name = Vector{Int8}(undef, 128)
    MPI.Sendrecv!(send_name, 0, 101, recv_name, 0, 100, inter_comm)

    partner_name = String(collect(Char.(Int.(recv_name)))) |> x->rstrip(x)
    println("[Jexpresso] ✓ Coupled with: $partner_name")
    println("")
    flush(stdout)
end

# ============================================================
# COUPLING COMPLETE - Now load and run Jexpresso simulation
# ============================================================

if lrank == 0
    println("[Jexpresso] Coupling initialization complete")
    println("[Jexpresso] Loading Jexpresso simulation modules...")
    println("")
    flush(stdout)
end

# Set up environment for Jexpresso to use local communicator
# Store coupling context in global variables that Jexpresso can access
global JEXPRESSO_COUPLED_MODE = true
global JEXPRESSO_WORLD_COMM = world
global JEXPRESSO_WORLD_RANK = wrank
global JEXPRESSO_LOCAL_COMM = local_comm
global JEXPRESSO_LOCAL_RANK = lrank
global JEXPRESSO_LOCAL_SIZE = lsize
global JEXPRESSO_INTER_COMM = inter_comm

# Now include Jexpresso with coupling already set up
include("src/Jexpresso.jl")

if lrank == 0
    println("[Jexpresso] Simulation modules loaded")
    println("[Jexpresso] Starting simulation with local communicator (size=$lsize)")
    println("")
    flush(stdout)
end

# Run Jexpresso simulation using local communicator
# The simulation will use JEXPRESSO_LOCAL_COMM for all internal MPI operations
# ... simulation code here ...

if lrank == 0
    println("")
    println("[Jexpresso] Simulation complete")
    println("[Jexpresso] Cleaning up coupling...")
    flush(stdout)
end

# Cleanup
MPI.free(inter_comm)
MPI.free(local_comm)

if lrank == 0
    println("[Jexpresso] Done")
    flush(stdout)
end

MPI.Finalize()
