# Dirichlet data g = u_ex on every boundary edge (tags left_right, bottom_top).
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    qbdy[1] = user_fft_exact(coords[1], coords[2])
    return nothing
end
