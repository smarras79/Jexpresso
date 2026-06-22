"""
    Dirichlet data for the non-constant diffusivity manufactured-solution test.

    g = u_ex on the whole boundary (every tag), consistent with the source f and
    the exact field qe = u_ex set in initialize.jl.
"""
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)

    qbdy[1] = ncd_u(coords[1], coords[2])

end
