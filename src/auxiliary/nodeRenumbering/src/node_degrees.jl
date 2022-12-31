# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

"""
    node_degrees(adjacency::Dict{Int, Vector{Int}})

Calculate the degrees for the adjacency graph nodes. Node degree = number of nodes a certain node is adjacent to.

# Examples

```jldoctest
julia> adjacency = Dict{Int, Vector{Int}}(
                   1 => [3, 8, 9],
                   2 => [3, 8, 7],
                   3 => [1, 2],
                   4 => [8, 9],
                   5 => [7, 8],
                   6 => [2, 7],
                   7 => [5, 2, 6],
                   8 => [1, 2, 4, 5],
                   9 => [1, 4]);

julia> node_degrees(adjacency::Dict{Int, Vector{Int}})
Dict{Int64,Int64} with 9 entries:
  7 => 3
  4 => 2
  9 => 2
  2 => 3
  3 => 2
  5 => 2
  8 => 4
  6 => 2
  1 => 3

```

# References

* Wikipedia contributors. "Degree (graph theory)". Wikipedia, The Free Encyclopedia. 24 Nov. 2016. Web. 18 Jul. 2017. https://en.wikipedia.org/wiki/Degree_(graph_theory)
"""
function node_degrees(adjacency::Dict{Int, Vector{Int}})
    degrees = Dict{Int, Int}()
    for (k, v) in adjacency
        degree = length(v)
        degrees[k] = degree
    end
    return degrees
end
