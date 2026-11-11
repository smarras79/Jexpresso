#---------------------------------------------------------------------------------
# What reaches the VTK files: every field the solver hands over (qvars in
# initialize.jl) — displacement, stress, the von Mises stress, the
# transformation strain (quadrature-point history projected onto the nodes),
# the thickness strain εzz of plane stress and the nodal reactions. (In the
# verification mode the error against the exact solution goes to
# solid_history.csv.)
#---------------------------------------------------------------------------------
function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    for i = 1:length(uout)
        uout[i] = u[i]
    end
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    user_uout!(ip, TOTAL(), uout, u, qe; kwargs...)
end
