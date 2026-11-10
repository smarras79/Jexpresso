#---------------------------------------------------------------------------------
# Sources (units of user_flux.jl, g = 1, z up):
#
#   gravity            S_ρw = −ρ,  S_E = −ρw
#   GLM damping        S_ψ = −(α_p c_h/Δh) ψ  (Dedner mixed cleaning)
#   velocity damping   Eqs. 21-22: every step v is divided by F = 1 + Δt N, where
#                      |κ| > κ_d for κ = x, y (κ_d = 3L_x/8, L_x the full width: the
#                      outer L_x/8 of each side) and z (κ_d = 7 max(z)/8), with
#                      N = (|κ| − κ_d)/(max|κ| − κ_d) rising linearly from 0 to 1
#                      at the boundary (in 1/t₀): to first order in Δt the source
#                      −N v below; the terms of x, y and z add up.
#                      In conserved form S_ρv = −N ρv and S_E = −N ρ|v|² (the
#                      kinetic energy removed; p is unchanged).
#---------------------------------------------------------------------------------
@inline function lk22_damping_rate(x, y, z)
    xd = lk22_damp_xd[]; zd = lk22_damp_zd[]; X = lk22_xmax[]; Z = lk22_zmax[]
    N = 0.0
    abs(x) > xd && (N += (abs(x) - xd)/(X - xd))
    abs(y) > xd && (N += (abs(y) - xd)/(X - xd))
    z > zd      && (N += (z - zd)/(Z - zd))
    return N
end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=9, x=0.0, y=0.0, z=0.0, xmax=0.0, xmin=0.0, zmax=0.0)

    for ieq = 1:neqs
        S[ieq] = 0.0
    end
    ρ = q[1]

    S[5] = -ρ                 # gravity on the vertical momentum
    S[4] = -q[5]              # work of gravity

    S[9] = -(lk22_glm_alpha[]*c_h_mhd[]/lk22_glm_dh[])*q[9]

    N = lk22_damping_rate(x, y, z)
    if N > 0.0
        S[2] -= N*q[2]
        S[3] -= N*q[3]
        S[5] -= N*q[5]
        S[4] -= N*(q[2]*q[2] + q[3]*q[3] + q[5]*q[5])/ρ
    end
end

function user_source!(S, q, qe, npoin::Int64, ::CL, ::PERT; neqs=9, x=0.0, y=0.0, z=0.0, xmax=0.0, xmin=0.0, zmax=0.0)
    error(" problems/MHD/fluxEmergenceLeake2022: PERT() solution variables are not supported.")
end
