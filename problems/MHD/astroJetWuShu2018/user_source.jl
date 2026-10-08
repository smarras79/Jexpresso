#---------------------------------------------------------------------------------
# GLM damping source (Dedner et al., JCP 175:645-673, 2002 — the "mixed"
# hyperbolic-parabolic GLM).
#
# The ideal MHD equations of this problem have no physical source terms: no
# gravity, no resistivity, no radiation. The only right-hand side is the
# divergence-cleaning damping that the GLM reformulation of ∇·B = 0 brings with
# it.
#
# The conservative GLM-MHD flux of user_flux.jl implements PURELY HYPERBOLIC
# cleaning: (ψ, ∇·B) form a linear wave system of speed c_h with no
# dissipation. In a DG setting the interface Riemann fluxes upwind-damp those
# waves; Jexpresso's continuous-Galerkin weak form has no interface
# dissipation, so the ψ waves alias into grid-scale (LGL odd-even) oscillations
# and grow. The standard cure is Dedner's mixed cleaning: damp ψ with
#
#     S_ψ = -(c_h² / c_p²) ψ = -(c_h / c_r) ψ,        c_p² = c_h c_r,
#
# with c_r ≈ 0.18 Dedner's recommended ratio. ONLY ψ is damped: the total energy
# carries ½ψ², so removing ψ while leaving E untouched turns the cleaned
# ψ-energy into heat, which is the entropy-consistent behaviour (Derigs et al.,
# JCP 364:420-467, 2018).
#
# Note the rate here: c_h/c_r ≈ 812/0.18 ≈ 4.5e3 per unit time, i.e. an
# e-folding in 2.2e-4 — about 1/10 of the run. That is stiff relative to
# nothing else in the problem but is still 440 explicit steps at Δt = 5e-7, so
# it integrates without any special treatment. JEXPRESSO_AJ_GLMCR raises c_r
# (weaker damping) if the ψ field ever needs to be left alone to diagnose it.
#
# NOTE: the non-conservative Galilean GLM transport (v·∇ψ) is omitted; see the
# header of user_flux.jl.
#---------------------------------------------------------------------------------

# Dedner damping ratio c_r (a Ref so it can be tuned from the REPL without
# re-including the case). Damping rate = c_h/c_r; larger c_r = weaker damping.
if !@isdefined(glm_cr_mhd)
    const glm_cr_mhd = Ref{Float64}(
        something(tryparse(Float64, get(ENV, "JEXPRESSO_AJ_GLMCR", "")), 0.18))
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
