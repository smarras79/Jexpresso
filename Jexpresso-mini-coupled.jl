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

world = MPI.COMM_WORLD
wrank = MPI.Comm_rank(world)
wsize = MPI.Comm_size(world)
println("[Jexpresso rank $wrank] World size: $wsize"); flush(stdout)

# Read APPID (0 or 1) from environment
appid = try parse(Int, get(ENV, "APPID", "-1")) catch; -1 end
if appid < 0
    if wrank == 0
        println("[Jexpresso] ERROR: APPID not set. Launch with -x APPID=0 (Fortran) and -x APPID=1 (Julia).")
    end
    MPI.Abort(world, 1)
end

# Split WORLD into per-app local comms
local_comm = MPI.Comm_split(world, appid, wrank)
lrank = MPI.Comm_rank(local_comm)
lsize = MPI.Comm_size(local_comm)

# --- Dynamic remote leader discovery via Allgather on WORLD ---
my = Int32(appid)
world_appids = Vector{Int32}(undef, wsize)
MPI.Allgather!(Ref(my), world_appids, world)   # each rank contributes 1 Int32

remote_leader_world = -1
for i in 1:wsize
    if world_appids[i] != my
        global remote_leader_world = i - 1  # Julia arrays 1-based, ranks 0-based
        break
    end
end
if remote_leader_world < 0
    if wrank == 0
        println("[Jexpresso] ERROR: did not find remote leader")
    end
    MPI.Abort(world, 2)
end

if wrank == 0
    println("[Jexpresso] world_size=$wsize appid=$appid remote_leader_world=$remote_leader_world"); flush(stdout)
end

# Create intercommunicator via ccall
libmpi = isdefined(MPI, :libmpi) ? MPI.libmpi : MPI.MPI_LIBRARY
newcomm_ref = Ref{MPI.MPI_Comm}()
const TAG = 12345

ierr = ccall((:MPI_Intercomm_create, libmpi), Cint,
             (MPI.MPI_Comm, Cint, MPI.MPI_Comm, Cint, Cint, Ref{MPI.MPI_Comm}),
             local_comm.val, Cint(0), world.val, Cint(remote_leader_world), Cint(TAG), newcomm_ref)
if ierr != 0
    error("MPI_Intercomm_create failed with error code $ierr")
end
inter_comm = MPI.Comm(newcomm_ref[])

# Prepare 128 signed bytes (space-padded) as Int8 → MPI_SIGNED_CHAR
name = "JULIA"
send = fill(Int8(' '), 128)
nb = min(length(codeunits(name)), 128)
@inbounds for i in 1:nb
    send[i] = Int8(codeunits(name)[i])
end

if lrank == 0
    recv = Vector{Int8}(undef, 128)
    # FIXED: Julia sends tag 101, receives tag 100 (opposite of Fortran)
    # Fortran sends tag 100, receives tag 101
    MPI.Sendrecv!(send, 0, 101, recv, 0, 100, inter_comm)
    partner = String(collect(Char.(Int.(recv)))) |> x->rstrip(x)
    println("[Jexpresso root] Partner name = ", partner); flush(stdout)
end

for i in 1:5
    MPI.Barrier(inter_comm)
    if lrank == 0
        println("[Jexpresso] Step $i/5"); flush(stdout)
    end
end

MPI.free(inter_comm)
MPI.free(local_comm)
MPI.Finalize()
println("[Jexpresso rank $wrank] Done"); flush(stdout)
