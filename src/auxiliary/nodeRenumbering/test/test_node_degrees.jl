# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

using NodeNumbering

adjacency = Dict(
    1 => [3, 8, 9],
    2 => [3, 8, 7, 6],
    3 => [1, 2],
    4 => [8, 9],
    5 => [7, 8],
    6 => [2, 7],
    7 => [5, 2, 6],
    8 => [1, 2, 4, 5],
    9 => [1, 4])

@testset "node degrees" begin
    E = Dict(
        1 => 3,
        2 => 4,
        3 => 2,
        4 => 2,
        5 => 2,
        6 => 2,
        7 => 3,
        8 => 4,
        9 => 2)
    result = NodeNumbering.node_degrees(adjacency)
    @test result == E
end
