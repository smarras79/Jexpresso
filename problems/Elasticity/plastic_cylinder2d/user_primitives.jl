#---------------------------------------------------------------------------------
# What reaches the VTK files: every field the solver hands over (qvars in
# initialize.jl) — displacement, stress, the von Mises stress and the
# accumulated plastic strain ēᵖ, the latter a quadrature-point history value
# projected onto the nodes. There is no closed form in the plastic range, so
# nothing is compared here (user_analytic.jl returns `nothing`).
#---------------------------------------------------------------------------------
function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    for i = 1:length(uout)
        uout[i] = u[i]
    end
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    user_uout!(ip, TOTAL(), uout, u, qe; kwargs...)
end
