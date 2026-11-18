#-------------------------------------------------------------------------------------------------
# Insulated ends: nothing to do. The Dirichlet hook must exist but sets no value, so no node is
# constrained, and with :bdy_fluxes => false no boundary flux is added.
#-------------------------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, qe, ::TOTAL)
    return nothing
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, qe, ::PERT)
    return nothing
end
