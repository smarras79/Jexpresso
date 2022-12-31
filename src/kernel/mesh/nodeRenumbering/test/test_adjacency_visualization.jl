# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

using NodeNumbering

RCM_adjacency = Dict{Int, Vector{Int}}(
                7 => [5, 6, 9],
                9 => [7, 8],
                4 => [2, 6],
                2 => [4, 3, 1],
                3 => [5, 6, 2, 1],
                8 => [6, 9],
                5 => [7, 3],
                6 => [7, 3, 8, 4],
                1 => [3, 2]);

@testset "Visualize an adjacency graph Dict." begin
    expected =
    [1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0;
     1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0;
     1.0  1.0  1.0  0.0  1.0  1.0  0.0  0.0  0.0;
     0.0  1.0  0.0  1.0  0.0  1.0  0.0  0.0  0.0;
     0.0  0.0  1.0  0.0  1.0  0.0  1.0  0.0  0.0;
     0.0  0.0  1.0  1.0  0.0  1.0  1.0  1.0  0.0;
     0.0  0.0  0.0  0.0  1.0  1.0  1.0  0.0  1.0;
     0.0  0.0  0.0  0.0  0.0  1.0  0.0  1.0  1.0;
     0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0]
    result = NodeNumbering.adjacency_visualization(RCM_adjacency)
    @test result == expected
end
