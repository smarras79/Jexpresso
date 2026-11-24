#---------------------------------------------------------------------------------
# GLM damping source (Dedner et al., JCP 175:645-673, 2002 - "mixed" GLM).
#
# The ideal MHD equations of Bormanis et al. (Eqs. 1-4) have no source
# terms: the only right-hand side here is the divergence-cleaning damping
# that the GLM reformulation of the ∇·B = 0 constraint brings with it.
#
# The conservative GLM-MHD flux implements PURELY HYPERBOLIC divergence
# cleaning: (ψ, ∇·B) form a linear wave system with speed c_h and no
# dissipation. In a DG setting the interface Riemann fluxes upwind-damp
# those waves; Jexpresso's continuous-Galerkin weak form has no interface
# dissipation, so on a doubly periodic domain the ψ waves keep bouncing,
# alias into grid-scale (LGL odd-even) oscillations and slowly grow, while
# every other field stays clean.
#
# The standard cure is Dedner's mixed hyperbolic-parabolic cleaning: damp ψ
# with the source
#
#     S_ψ = -(c_h² / c_p²) ψ = -(c_h / c_r) ψ,       c_p² = c_h c_r,
#
# where c_r ~ 0.18 is Dedner's recommended ratio. Only ψ is damped: since
# the total energy carries ½ψ², removing ψ while leaving E untouched turns
# the cleaned ψ-energy into heat, which is the entropy-consistent behavior
# (Derigs et al., JCP 364:420-467, 2018).
#
# NOTE: the non-conservative Galilean GLM transport (v·∇ψ) is omitted; see
# user_flux.jl.
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
