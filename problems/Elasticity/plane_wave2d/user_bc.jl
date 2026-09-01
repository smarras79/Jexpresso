#---------------------------------------------------------------------------------
# There are no boundaries in this problem.
#
# The mesh tags both pairs of sides "periodicx" and "periodicy", and
# build_custom_bcs_dirichlet!(::NSD_2D, ...) skips periodic edges outright —
# this hook is never reached. Leaving qbdy untouched (it arrives prefilled with
# a sentinel) would be a no-op even if it were.
#
# That is the point of shipping this case: it exercises the 2D flux operator,
# the metric terms and the periodic node matching against an exact solution,
# with the boundary treatment taken entirely out of the picture.
#---------------------------------------------------------------------------------
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    return qbdy
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                            qbdy::AbstractArray, nx, ny, qe, ::PERT)
    return qbdy
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords,
                         t::AbstractFloat, tag::String, inputs)
    return zeros(size(q, 2), 1)
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords,
                         t::AbstractFloat, inputs)
    return zeros(size(q, 2), 1)
end

function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)
    T = eltype(q)
    return T(q[1]), T(q[2]), T(q[3]), T(q[4]), T(q[5])
end
