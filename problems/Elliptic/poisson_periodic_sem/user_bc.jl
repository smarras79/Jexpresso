# Doubly periodic domain: every boundary edge carries a periodicx/periodicy tag,
# whose nodes are identified with their periodic images when the mesh is read.
# There is therefore no Dirichlet data to impose; this hook is never called.
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    return nothing
end
