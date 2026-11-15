#---------------------------------------------------------------------------------
# Sources (units of user_flux.jl, g = 1, y up), as problems/MHD/fluxEmergenceLeake2022:
#
#   gravity            S_ρv = −ρ,  S_E = −ρv
#   GLM damping        S_ψ = −(α_p c_h/Δh) ψ  (Dedner mixed cleaning)
#   velocity damping   the paper's Eqs. 21-22 in the plane: Dv/Dt gets −N v where
#                      |x| > x_d (side, x_d = 3/8 of the full width) and y > y_d (top,
#                      y_d = 7/8 of the top height), with N = (|x| − x_d)/(max|x| − x_d)
#                      and (y − y_d)/(max y − y_d) rising linearly to 1 at the boundary
#                      (in 1/t₀); the two add up. In conserved form S_ρv = −N ρv for the
#                      three components and S_E = −N ρ|v|² (the kinetic energy removed;
#                      p is unchanged).
#---------------------------------------------------------------------------------
@inline function lk2d_damping_rate(x, y)
    xd = lk2d_damp_xd[]; yd = lk2d_damp_yd[]
    N = 0.0
    abs(x) > xd && (N += (abs(x) - xd)/(lk2d_xmax[] - xd))
    y > yd      && (N += (y - yd)/(lk2d_ymax[] - yd))
    return N
end

function user_source!(S,
                      q,
                      qe,
                      npoin::TInt,
                      ::CL, ::TOTAL;
                      neqs=9, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    for ieq = 1:neqs
        S[ieq] = 0.0
    end
    ρ = q[1]

    S[3] = -ρ                 # gravity on the vertical momentum
    S[4] = -q[3]              # work of gravity

    S[9] = -(lk2d_glm_alpha[]*c_h_mhd[]/lk2d_glm_dh[])*q[9]

    N = lk2d_damping_rate(x, y)
    if N > 0.0
        S[2] -= N*q[2]
        S[3] -= N*q[3]
        S[5] -= N*q[5]
        S[4] -= N*(q[2]*q[2] + q[3]*q[3] + q[5]*q[5])/ρ
    end
end

function user_source!(S, q, qe, npoin::TInt, ::CL, ::PERT; neqs=9, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)
    error(" problems/MHD/fluxEmergenceLeake2022_2D: PERT() solution variables are not supported.")
end
