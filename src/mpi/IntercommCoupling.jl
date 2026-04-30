# Intercommunicator Coupling Module for Jexpresso
# Provides true independence between coupled codes using MPI intercommunicators

module IntercommCoupling

using MPI

"""
    IntercommContext

Context for intercommunicator-based coupling.

Fields:
- `world_comm`: Original MPI_COMM_WORLD
- `world_rank`: Rank in MPI_COMM_WORLD
- `world_size`: Size of MPI_COMM_WORLD
- `appid`: Application ID from environment (0, 1, 2, ...)
- `local_comm`: Local communicator for this application
- `local_rank`: Rank within local communicator
- `local_size`: Size of local communicator
- `inter_comm`: Intercommunicator to other application(s)
- `remote_leader_world`: World rank of remote leader
"""
struct IntercommContext
    world_comm::MPI.Comm
    world_rank::Int
    world_size::Int
    appid::Int
    local_comm::MPI.Comm
    local_rank::Int
    local_size::Int
    inter_comm::MPI.Comm
    remote_leader_world::Int
end

"""
    initialize_intercomm_coupling(; code_name::String="Jexpresso")

Initialize intercommunicator-based coupling.

Requires APPID environment variable to be set (0 for first code, 1 for second, etc.)

Returns IntercommContext with local and inter communicators.

# Example
```julia
ctx = IntercommCoupling.initialize_intercomm_coupling(code_name="Jexpresso")

# Use ctx.local_comm for internal Jexpresso operations
# Use ctx.inter_comm to communicate with external code
```
"""
function initialize_intercomm_coupling(; code_name::String="Jexpresso")
    world = MPI.COMM_WORLD
    wrank = MPI.Comm_rank(world)
    wsize = MPI.Comm_size(world)

    if wrank == 0
        @info "Initializing intercommunicator coupling for $code_name"
    end

    # Read APPID from environment
    appid = try
        parse(Int, get(ENV, "APPID", "-1"))
    catch
        -1
    end

    if appid < 0
        if wrank == 0
            error("APPID environment variable not set. Launch with -x APPID=0 (first code) and -x APPID=1 (second code)")
        end
        MPI.Abort(world, 1)
    end

    if wrank == 0
        @info "APPID=$appid, world_rank=$wrank, world_size=$wsize"
    end

    # Split MPI_COMM_WORLD by APPID to create local communicator
    local_comm = MPI.Comm_split(world, appid, wrank)
    lrank = MPI.Comm_rank(local_comm)
    lsize = MPI.Comm_size(local_comm)

    if wrank == 0
        @info "Local communicator: rank=$lrank/$lsize"
    end

    # Discover remote leader via Allgather on WORLD
    my_appid = Int32(appid)
    world_appids = Vector{Int32}(undef, wsize)
    MPI.Allgather!(Ref(my_appid), world_appids, world)

    remote_leader_world = -1
    for i in 1:wsize
        if world_appids[i] != my_appid
            remote_leader_world = i - 1  # Convert to 0-based rank
            break
        end
    end

    if remote_leader_world < 0
        if wrank == 0
            error("Could not find remote leader (no other APPID found in MPI_COMM_WORLD)")
        end
        MPI.Abort(world, 2)
    end

    if wrank == 0
        @info "Remote leader in world: rank $remote_leader_world"
    end

    # Create intercommunicator
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
        @info "Intercommunicator created successfully"
    end

    # Exchange application names for verification
    if lrank == 0
        send_name = fill(Int8(' '), 128)
        name_bytes = codeunits(code_name)
        nb = min(length(name_bytes), 128)
        @inbounds for i in 1:nb
            send_name[i] = Int8(name_bytes[i])
        end

        recv_name = Vector{Int8}(undef, 128)

        # Fixed tags: Julia sends 101, receives 100
        MPI.Sendrecv!(send_name, 0, 101, recv_name, 0, 100, inter_comm)

        partner_name = String(collect(Char.(Int.(recv_name)))) |> x->rstrip(x)
        @info "Coupled with: $partner_name"
    end

    return IntercommContext(
        world, wrank, wsize,
        appid,
        local_comm, lrank, lsize,
        inter_comm,
        remote_leader_world
    )
end

"""
    finalize_intercomm_coupling(ctx::IntercommContext)

Clean up intercommunicator coupling context.
"""
function finalize_intercomm_coupling(ctx::IntercommContext)
    if ctx.world_rank == 0
        @info "Finalizing intercommunicator coupling"
    end

    MPI.free(ctx.inter_comm)
    MPI.free(ctx.local_comm)
end

end # module IntercommCoupling
