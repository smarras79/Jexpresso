# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

"""
    RCM(adjacency::Dict{Int, Vector{Int}}, degrees::Dict{Int, Int}, P::Int)

Calculate the Reverse Cuthill-McKee Algorithm for the adjacency graph.

# Examples

```jldoctest
julia> adjacency = Dict(
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

julia> degrees = Dict(
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

julia> P = 9;

julia> Q = 12;

julia> RCM(adjacency, degrees, P)
16-element Array{Int64,1}:
 16
 5
 14
 13
 11
 4
 12
 6
 7
 2
 10
 3
 8
 1
 15
 9

```

# References

* Wikipedia contributors. "Cuthill–McKee algorithm". Wikipedia, The Free Encyclopedia. Wikipedia, The Free Encyclopedia, 11 Jul. 2017. Web. 17 Jul. 2017. https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm
"""
function RCM(adjacency::Dict{Int, Vector{Int}}, degrees::Dict{Int, Int}, P::Int, Q::Int)
    n = length(adjacency)
    R = Int[P]
    for i=1:n
        if i > length(R)
            push!(R, Q)
            A = sort(adjacency[Q], by=j->degrees[j])
            for b in A
                if !(b in R)
                    push!(R, b)
                end
            end
        else
            t = sort(adjacency[R[i]], by=j->degrees[j])
            for T in t
                if !(T in R)
                    push!(R, T)
                end
            end
        end
    end
    new_order = reverse(R)
    return new_order
end
