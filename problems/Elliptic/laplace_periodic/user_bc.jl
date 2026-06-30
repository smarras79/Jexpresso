# Doubly-periodic problem: there are no Dirichlet boundaries, so this hook is a
# no-op. (Periodicity is carried by the mesh; the FFT solver is periodic by
# construction.) The function is kept so the case has the same structure as the
# other Elliptic problems and the setup pipeline can include it.
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    return qbdy
end
