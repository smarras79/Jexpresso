#---------------------------------------------------------------------------------
# What reaches the VTK files: every field the solver hands over (qvars in
# initialize.jl) — displacement, stress, the von Mises stress, the
# accumulated plastic strain ēᵖ (a quadrature-point history value projected
# onto the nodes) and the nodal reactions rx, ry.
#---------------------------------------------------------------------------------
function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    for i = 1:length(uout)
        uout[i] = u[i]
    end
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    user_uout!(ip, TOTAL(), uout, u, qe; kwargs...)
end
