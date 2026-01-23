#!/usr/bin/env julia
# jexpresso_gather_continuous.jl
# Julia simulation that continues running after gather-based coupling
# Simulates Jexpresso running in parallel with external code

using MPI
using ArgParse

function parse_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--code-name"
        help = "Application name to send during coupling"
        arg_type = String
        default = "Jexpresso"

        "--max-time"
        help = "Maximum simulation time in seconds"
        arg_type = Float64
        default = 10.0
    end
    return parse_args(s)
end

args = parse_args()
code_name = args["code-name"]
max_time = args["max-time"]

println("Jexpresso (continuous mode) starting...")
flush(stdout)

MPI.Init()

world = MPI.COMM_WORLD
rank = MPI.Comm_rank(world)
size = MPI.Comm_size(world)

println("[Jexpresso rank $rank] World size: $size")
flush(stdout)

# ============================================================
# COUPLING INITIALIZATION
# ============================================================

# Participate in MPI_Comm_split with MPI_UNDEFINED
if rank == 0
    println("[Jexpresso rank 0] Participating in MPI_Comm_split")
    flush(stdout)
end

MPI_UNDEFINED = Int32(-32766)
par_comm_final = MPI.Comm_split(world, MPI_UNDEFINED, rank)

if rank == 0
    println("[Jexpresso rank 0] MPI_Comm_split completed")
    flush(stdout)
end

# Prepare application name for MPI_Gather (128 bytes)
send_buf = fill(Int8(' '), 128)
nb = min(length(codeunits(code_name)), 128)
@inbounds for i in 1:nb
    send_buf[i] = Int8(codeunits(code_name)[i])
end

# Participate in MPI_Gather
if rank == 0
    println("[Jexpresso rank 0] Gathering application names")
    flush(stdout)
end

MPI.Gather!(send_buf, nothing, 0, world)

if rank == 0
    println("[Jexpresso rank 0] Coupling initialization complete. Starting simulation...")
    println("")
    flush(stdout)
end

# ============================================================
# SIMULATION LOOP - Runs in parallel with Fortran code
# ============================================================

start_time = time()
time_step = 0

while true
    current_time = time() - start_time
    if current_time >= max_time
        break
    end

    time_step += 1

    # Simulate work
    sleep(1)

    # Synchronize with all ranks (including Fortran ranks)
    MPI.Barrier(world)

    if rank == 0
        println("[Jexpresso rank 0] Time step $time_step, elapsed: $(round(current_time, digits=2))s")
        flush(stdout)
    end
end

# ============================================================
# FINALIZATION - Must synchronize with Fortran before finalizing
# ============================================================

if rank == 0
    println("")
    println("[Jexpresso rank 0] Simulation complete after $time_step steps")
    println("[Jexpresso rank 0] Waiting for all processes before finalize")
    flush(stdout)
end

# Final barrier to ensure all ranks finish together
MPI.Barrier(world)

MPI.Finalize()

if rank == 0
    println("[Jexpresso rank 0] Done")
    flush(stdout)
end
