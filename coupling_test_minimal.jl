#!/usr/bin/env julia
# coupling_test_minimal.jl
# Minimal test for gather-based coupling without loading full Jexpresso
# Use this for testing coupling with external codes that exit quickly

using MPI
using ArgParse

function parse_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--code-name"
        help = "Application name to send during coupling"
        arg_type = String
        default = "Jexpresso"
    end
    return parse_args(s)
end

args = parse_args()
code_name = args["code-name"]

println("Minimal coupling test starting...")
flush(stdout)

MPI.Init()

world = MPI.COMM_WORLD
rank = MPI.Comm_rank(world)
size = MPI.Comm_size(world)

println("[Jexpresso rank $rank] World size: $size")
flush(stdout)

# Participate in MPI_Comm_split with MPI_UNDEFINED
println("[Jexpresso rank $rank] Participating in MPI_Comm_split")
flush(stdout)

MPI_UNDEFINED = Int32(-32766)
par_comm_final = MPI.Comm_split(world, MPI_UNDEFINED, rank)

println("[Jexpresso rank $rank] MPI_Comm_split completed")
flush(stdout)

# Prepare application name for MPI_Gather (128 bytes)
send_buf = fill(Int8(' '), 128)
nb = min(length(codeunits(code_name)), 128)
@inbounds for i in 1:nb
    send_buf[i] = Int8(codeunits(code_name)[i])
end

# Participate in MPI_Gather
println("[Jexpresso rank $rank] Participating in MPI_Gather to root (rank 0)")
flush(stdout)

MPI.Gather!(send_buf, nothing, 0, world)

println("[Jexpresso rank $rank] MPI_Gather completed")
flush(stdout)

if rank == 0
    println("✓ Coupling test successful: Participated in split and gather, sending '$code_name'")
end

MPI.Finalize()
println("[Jexpresso rank $rank] Done")
flush(stdout)
