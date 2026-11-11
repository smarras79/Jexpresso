#---------------------------------------------------------------------------------
# What reaches the VTK files. The quasi-static solver hands write_output its six
# nodal fields (ux, uy, sxx, syy, sxy, szz — qvars in initialize.jl) and the
# closed-form solution at the same nodes as `qe` (user_analytic.jl, Heaviside
# load). Written, in the order of qoutvars:
#
#   ux, uy                 displacement (on the x axis ux IS u_r)
#   sxx, syy, sxy, szz     stress: the in-plane ones stay the Lamé stresses at
#                          every t, only szz creeps
#   ux_exact, uy_exact     the closed form
#   ux_err, uy_err         SEM minus closed form. With JEXPRESSO_VE_LOAD =
#                          first_increment this is not an error: the closed form
#                          is for a load that is on at t = 0, and the difference
#                          is the first-increment ramp.
#---------------------------------------------------------------------------------
function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    for i = 1:6
        uout[i] = u[i]
    end
    uout[7]  = qe[1]
    uout[8]  = qe[2]
    uout[9]  = u[1] - qe[1]
    uout[10] = u[2] - qe[2]
end

function user_uout!(ip, ::PERT, uout, u, qe; kwargs...)
    user_uout!(ip, TOTAL(), uout, u, qe; kwargs...)
end
