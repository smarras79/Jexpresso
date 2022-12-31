# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

using NodeNumbering

@testset "Calculate bandwidth of adjacency graph" begin

    G = Dict(
        1 => [3, 8, 9],
        2 => [3, 8, 7],
        3 => [1, 2],
        4 => [8, 9],
        5 => [7, 8],
        6 => [2, 7],
        7 => [5, 2, 6],
        8 => [1, 2, 4, 5],
        9 => [1, 4])

    bw = NodeNumbering.calculate_bandwidth(G)

    @test bw == 17

end
