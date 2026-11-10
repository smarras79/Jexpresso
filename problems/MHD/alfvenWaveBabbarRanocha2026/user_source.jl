#---------------------------------------------------------------------------------
# GLM damping source (Dedner et al., JCP 175:645-673, 2002 - "mixed" GLM),
# the one of problems/MHD/smoothVortex and orszagTangBormanis2024.
#
# Ideal MHD has no source here: the only right-hand side is the damping of the
# divergence-cleaning field,
#
#     S_ψ = -(c_h / c_r) ψ,       c_r = 0.18 (Dedner's ratio).
#
# Jexpresso's continuous-Galerkin weak form has no interface dissipation, so
# without it the ψ waves of the discretization's own ∇·B error keep bouncing
# on the doubly periodic box and alias into grid-scale oscillations. The exact
# solution has ψ = 0, so S_ψ vanishes on it and does not bias the error.
#---------------------------------------------------------------------------------

# Dedner damping ratio c_r (a Ref so it can be tuned from the REPL without
# re-including the case). Damping rate = c_h/c_r; larger c_r = weaker damping.
if !@isdefined(glm_cr_mhd)
    const glm_cr_mhd = Ref{Float64}(0.18)
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

    # ψ damping: S_ψ = -(c_h/c_r) ψ
    S[9] = -(c_h_mhd[]/glm_cr_mhd[])*q[9]

end
