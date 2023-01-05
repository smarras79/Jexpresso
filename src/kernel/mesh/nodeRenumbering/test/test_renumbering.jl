# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

using NodeNumbering

@testset "Renumbering the RCM ordered nodes." begin

    neworder = Int[6, 7, 2, 5, 3, 8, 1, 4, 9]

    expected = Dict(
               7 => 2,
               9 => 9,
               4 => 8,
               2 => 3,
               3 => 5,
               5 => 4,
               8 => 6,
               6 => 1,
               1 => 7)
    result = NodeNumbering.renumbering(neworder)
    @test result == expected
end
