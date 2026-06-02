"""
    CompEuler/case1 uses :lperiodic_1d => true, so apply_periodicity!
    is invoked instead of the Dirichlet path; the user_bc_dirichlet!
    methods below are kept as no-ops so the file remains structurally
    consistent with the other 1D problems.
"""
function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::TOTAL)
    nothing
end

function user_bc_dirichlet!(q, coords, t, tag::String, qbdy, qe, ::PERT)
    nothing
end

function user_bc_neumann(q::AbstractArray, gradq, coords, t, inputs)
    flux = zeros(size(q, 2), 1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, lpert)
    T = eltype(q)
    return T(q[1]), T(q[2]), T(q[3])
end
