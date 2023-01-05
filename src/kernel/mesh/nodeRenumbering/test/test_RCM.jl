# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

using NodeNumbering

P = 9

Q = 12

adjacency = Dict(
                   1 => [3, 8, 9],
                   2 => [6, 7, 3, 8],
                   3 => [1, 2],
                   4 => [12, 13],
                   5 => [13, 16],
                   6 => [2, 7],
                   7 => [2, 6, 10],
                   8 => [1, 2, 10, 15],
                   9 => [1, 15],
                  10 => [8, 7],
                  11 => [12, 13, 14],
                  12 => [4, 11],
                  13 => [5, 16, 4, 11],
                  14 => [11, 16],
                  15 => [9, 8],
                  16 => [13, 5, 14]);

degrees = Dict(
                 1 => 3,
                 2 => 4,
                 3 => 2,
                 4 => 2,
                 5 => 2,
                 6 => 2,
                 7 => 3,
                 8 => 4,
                 9 => 2,
                10 => 2,
                11 => 3,
                12 => 2,
                13 => 4,
                14 => 2,
                15 => 2,
                16 => 3);

@testset "The Reverse Cuthill-McKee algorithm" begin
    expected = [16,5,14,13,11,4,12,6,7,2,10,3,8,1,15,9]
    result = NodeNumbering.RCM(adjacency, degrees, P, Q)
    @test result == expected
end
