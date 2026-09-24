#---------------------------------------------------------------------------------
# PODbenchmark — the flux of the linear advection equation,
#
#   ∂u/∂t + c ∂u/∂x = 0 ,   F = c u ,   c = 1.
#
# c MUST agree with the `c` of user_inputs.jl and initialize.jl: the closed-form
# POD of the README is the decomposition of a wave translating at that speed
# over exactly one revolution of the domain, and the deck derives the revolution
# time T = L/c from it.
#---------------------------------------------------------------------------------
function user_flux!(F, G, SD::NSD_1D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ip=1)

    c = 1.0
    F[1] = c*q[1]

end
