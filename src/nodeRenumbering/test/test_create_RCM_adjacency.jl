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
            9 => [1, 4]);

finalorder =  Dict(
              7 => 2,
              9 => 9,
              4 => 8,
              2 => 3,
              3 => 5,
              5 => 4,
              8 => 6,
              6 => 1,
              1 => 7);

@testset "Create RCM adjacency" begin
    expected =  Dict(
                7 => [5, 6, 9],
                9 => [7, 8],
                4 => [2, 6],
                2 => [4, 3, 1],
                3 => [5, 6, 2, 1],
                8 => [6, 9],
                5 => [7, 3],
                6 => [7, 3, 8, 4],
                1 => [3, 2])
    result = NodeNumbering.create_RCM_adjacency(adjacency, finalorder)
    @test result == expected
end
