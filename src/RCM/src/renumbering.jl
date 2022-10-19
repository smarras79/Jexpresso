# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

"""
    renumbering(neworder::Array{Int})

Renumber the RCM ordered nodes. A part of the Reverse Cuthill-McKee Algorithm.

# Examples

```jldoctest
julia> neworder = Int[6, 7, 2, 5, 3, 8, 1, 4, 9];

julia> renumbering(neworder)
Dict{Int64,Int64} with 9 entries:
  7 => 2
  9 => 9
  4 => 8
  2 => 3
  3 => 5
  5 => 4
  8 => 6
  6 => 1
  1 => 7

```

# References

* Wikipedia contributors. "Cuthill–McKee algorithm". Wikipedia, The Free Encyclopedia. Wikipedia, The Free Encyclopedia, 11 Jul. 2017. Web. 20 Jul. 2017. https://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm
"""
function renumbering(neworder::Array{Int})
    finalnodes = Dict{Int, Int}()
    n=length(neworder)
    for i=1:n
      neworder[i]
      finalnodes[i] = neworder[i]
    end
Final = Dict((v => k for (k,v) in finalnodes))
return Final
end
