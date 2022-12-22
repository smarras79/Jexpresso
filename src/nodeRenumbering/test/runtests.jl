# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

using Base.Test
using NodeNumbering

test_files = String[]
push!(test_files, "test_calc_bw.jl")
push!(test_files, "test_create_adjacency_graph.jl")
push!(test_files, "test_node_degrees.jl")
push!(test_files, "test_RCM.jl")
push!(test_files, "test_create_RCM_adjacency.jl")
push!(test_files, "test_renumbering.jl")
push!(test_files, "test_adjacency_visualization.jl")

using TimerOutputs
const to = TimerOutput()
@testset "NodeNumbering.jl" begin
    for fn in test_files
        timeit(to, fn) do
            include(fn)
        end
    end
end
println()
println("Test statistics:")
println(to)
