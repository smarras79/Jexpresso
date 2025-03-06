# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

using NodeNumbering



@testset "create adjacency graph" begin
    expected = Dict(
        1 => [3, 8, 9],
        2 => [3, 8, 7, 6],
        3 => [1, 2],
        4 => [8, 9],
        5 => [7, 8],
        6 => [2, 7],
        7 => [5, 2, 6],
        8 => [1, 2, 4, 5],
        9 => [1, 4])

       elements=Dict(
            1 => [9, 1, 8, 4],
            2 => [1, 3, 2, 8],
            3 => [8, 2, 7, 5],
            4 => [2, 6, 7])

        element_types = Dict{Int, Symbol}(
            1 => :Quad4,
            2 => :Quad4,
            3 => :Quad4,
            4 => :Tri3)

    result = NodeNumbering.create_adjacency_graph(elements, element_types)
    result_bool = true
    for (k, v) in expected
        if (length(intersect(v, result[k]))==length(v))==false
            result_bool = false
            break
        end
    end
    #@test result == expected
    @test result_bool==true
end
