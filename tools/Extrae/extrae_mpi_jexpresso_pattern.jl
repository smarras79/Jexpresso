#==============================================================================
# extrae_mpi_jexpresso_pattern.jl  --  MPI Extrae.jl instrumentation example
#
# A small, self-contained distributed kernel that mirrors the communication
# pattern Jexpresso uses in its real solver (see
# src/kernel/mpi/mpi_communications.jl):
#
#   * the global 1-D domain is partitioned across MPI ranks,
#   * every iteration computes a local stencil update (the "RHS"),
#   * neighbouring ranks exchange halo/ghost values with non-blocking
#     MPI.Isend / MPI.Irecv! + MPI.Waitall  (exactly the Jexpresso idiom),
#   * a global residual is reduced with MPI.Allreduce.
#
# This is deliberately the same shape of workload analysed in the Extrae.jl
# paper (arXiv:2504.12087v1, Section 4), where MPI_Waitall/MPI_Allreduce show
# up as the communication hot spots in the Paraver timeline.
#
# Each phase is wrapped in an Extrae user-code region and tagged with an event
# so Paraver shows a labelled timeline per rank.
#
# Run on your laptop (a 4-core MacBook Air can host up to ~8 ranks; 4 is a
# good default):
#
#     julia --project=. -e 'using MPI; run(`$(mpiexec()) -n 4 \
#         $(Base.julia_cmd()) --project=. tools/Extrae/extrae_mpi_jexpresso_pattern.jl`)'
#
#   ...or simply:   ./tools/Extrae/run_extrae_example.sh 4
#
# To get a REAL Paraver trace you must run on Linux with Extrae installed and
# MPI auto-instrumentation enabled via LD_PRELOAD / MPIPreferences (see README).
#==============================================================================

using MPI

include(joinpath(@__DIR__, "ExtraeShim.jl"))
using .ExtraeShim

# --- event type codes (extrae_type_t / UInt32) ------------------------------
const EV_PHASE = 90000   # which phase of the iteration we are in
const EV_ITER  = 90001   # current iteration number
const EV_RES   = 90002   # (scaled) global residual

# phase value codes
const PHASE_END      = 0
const PHASE_RHS      = 1
const PHASE_HALO     = 2
const PHASE_ALLREDUCE = 3

function register_events(rank)
    register(EV_PHASE, "Solver phase",
             (PHASE_END, PHASE_RHS, PHASE_HALO, PHASE_ALLREDUCE),
             ("idle", "compute RHS", "halo exchange", "allreduce"))
    register(EV_ITER, "Iteration")
    register(EV_RES,  "Global residual (x1e6)")
    return nothing
end

"""
    halo_exchange!(u_local, left, right, comm, rank, nranks)

Exchange one ghost cell with each neighbour using the same non-blocking
idiom as Jexpresso's assembler (Isend/Irecv!/Waitall). `u_local` carries one
ghost cell on each side: index 1 (left ghost) and end (right ghost). Returns
nothing; ghosts are updated in place.
"""
function halo_exchange!(u_local, comm, rank, nranks)
    n = length(u_local)
    left  = rank == 0          ? MPI.PROC_NULL : rank - 1
    right = rank == nranks - 1 ? MPI.PROC_NULL : rank + 1

    # interior first/last (the values our neighbours need)
    send_left  = u_local[2]
    send_right = u_local[n - 1]
    recv_left  = Ref(zero(eltype(u_local)))
    recv_right = Ref(zero(eltype(u_local)))

    reqs = MPI.Request[]
    push!(reqs, MPI.Irecv!(recv_left,  comm; source = left,  tag = 1))
    push!(reqs, MPI.Irecv!(recv_right, comm; source = right, tag = 0))
    push!(reqs, MPI.Isend(Ref(send_right), comm; dest = right, tag = 1))
    push!(reqs, MPI.Isend(Ref(send_left),  comm; dest = left,  tag = 0))
    MPI.Waitall(reqs)

    if left != MPI.PROC_NULL
        u_local[1] = recv_left[]
    end
    if right != MPI.PROC_NULL
        u_local[n] = recv_right[]
    end
    return nothing
end

function main()
    MPI.Init()
    comm   = MPI.COMM_WORLD
    rank   = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)

    # In a Distributed.jl/MPI run Extrae maps ranks to Paraver TASK objects
    # automatically when the MPI library is instrumented. init() is still the
    # call that opens the trace.
    init()
    register_events(rank)

    # ---- problem setup: 1-D heat equation, explicit Jacobi smoothing -------
    nglobal = 4096
    base    = div(nglobal, nranks)
    nlocal  = base + (rank < nglobal % nranks ? 1 : 0)

    # local array with one ghost cell on each side
    u = zeros(Float64, nlocal + 2)
    # initial condition: a bump that depends on rank so there is real data flow
    for i in 2:(nlocal + 1)
        u[i] = sin(2pi * (rank * base + (i - 1)) / nglobal)
    end
    unew = copy(u)

    niter = 200
    if rank == 0
        println("Extrae MPI example: $nranks ranks, $nglobal global points, " *
                "$niter iterations (active trace = $(is_active())).")
    end

    for it in 1:niter
        emit(EV_ITER, it)

        # ---- halo exchange (communication) --------------------------------
        emit(EV_PHASE, PHASE_HALO)
        @user_function halo_exchange!(u, comm, rank, nranks)

        # ---- local RHS / stencil update (computation) ---------------------
        emit(EV_PHASE, PHASE_RHS)
        local_res = 0.0
        @user_function begin
            @inbounds for i in 2:(nlocal + 1)
                unew[i] = 0.5 * (u[i - 1] + u[i + 1])
                d = unew[i] - u[i]
                local_res += d * d
            end
            # swap interior; keep ghosts to be refreshed next exchange
            @inbounds for i in 2:(nlocal + 1)
                u[i] = unew[i]
            end
        end

        # ---- global residual reduction (collective) -----------------------
        emit(EV_PHASE, PHASE_ALLREDUCE)
        global_res = @user_function MPI.Allreduce(local_res, +, comm)
        emit(EV_RES, round(Int, sqrt(global_res) * 1e6))

        emit(EV_PHASE, PHASE_END)   # back to "idle" between iterations

        if rank == 0 && (it % 50 == 0 || it == niter)
            println("  iter $it  ||residual|| = $(sqrt(global_res))")
        end
    end

    finish()

    MPI.Barrier(comm)
    if rank == 0
        if is_active()
            println("Done. Extrae trace written; merge with mpi2prv (see README).")
        else
            println("Done (no-op shim mode: no Paraver trace produced on this " *
                    "platform). The instrumented MPI code ran correctly on all ranks.")
        end
    end
    MPI.Finalize()
end

main()
