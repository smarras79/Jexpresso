#---------------------------------------------------------------------------------
# Source terms of the flux-emergence problem (Son, Jang & Magara 2025, ApJS
# 277:46), i.e. the vector S of their Eq. (7),
#
#   S = (0, ρ g_x, ρ g_y, ρ g_z, 0, 0, 0, ρ (V·g), -(c_h²/c_p²) ψ)ᵀ
#
# mapped onto Jexpresso's slot ordering (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)
# with g = (0, -g₀, 0), g₀ = 1/γ in the nondimensional units of the paper's
# Table 1 (H₀ = C_s = ρ₀ = 1), plus the absorbing layer of the top boundary
# (paper Section 2.1: "free conditions with an absorbing layer at the top
# boundary ... to prevent artificial wave reflections (Machida & Matsumoto
# 2003)").
#
# 1. Gravity:          S[3] = -ρ g₀,   S[4] = -ρ v g₀.
#
# 2. GLM ψ damping (mixed hyperbolic-parabolic cleaning, Dedner et al. 2002):
#
#        S_ψ = -(c_h²/c_p²) ψ
#
#    parametrized as in the paper (its Eq. 12, after Mignone & Tzeferacos
#    2010) by the dimensionless α_p = Δh c_h / c_p², so that
#
#        c_h²/c_p² = α_p c_h / Δh,        Δh = smallest mesh spacing.
#
#    The paper adopts α_p = 0.2 for all WENO schemes (its Section 4.2). Δh is
#    the smallest LGL nodal spacing of the mesh, measured in initialize.jl.
#
# 3. Well-balanced correction: the LGL interpolant of the magnetostatic
#    state is not a discrete equilibrium where its tanh transitions are
#    under-resolved; the residual of the vertical momentum balance of the
#    initial state, tabulated against height by initialize.jl, is subtracted
#    as a static body force (S[3] += -R(z)). See initialize.jl, item 5.
#
# 4. Absorbing layer: Rayleigh damping of the departure from the initial
#    (magnetostatic) state above z_s,
#
#        S[i] -= σ(z) (q[i] - q_e[i]),   σ(z) = σ_max sin²[ (π/2) (z - z_s)/(Z_max - z_s) ],
#
#    applied to every field. Machida & Matsumoto (2003) do not tabulate their
#    layer; z_s = 30 H₀ (a 5 H₀ = one coronal-sound-crossing-time deep layer)
#    and σ_max = 2/τ₀ absorb an outgoing coronal wave (speed ≈ 5 C_s) over
#    ~2 e-folding times before it can reflect. Both are Refs, tunable from the
#    REPL; z_s is compared against the GLOBAL top of the mesh recorded by
#    initialize.jl (the `ymax` keyword is the rank-local extent under MPI).
#
# The CL/TOTAL signature is the one the 2D kernel calls (rhs.jl,
# _build_rhs!): x, y are the point coordinates.
#---------------------------------------------------------------------------------

# Mignone's damping parameter α_p (paper: 0.2) and the smallest mesh spacing
# Δh it refers to (filled by initialize.jl).
if !@isdefined(glm_alpha_p_mhd)
    const glm_alpha_p_mhd = Ref{Float64}(0.2)
end
if !@isdefined(glm_dh_mhd)
    const glm_dh_mhd = Ref{Float64}(1.0)
end

# Absorbing layer: base height z_s, top of the domain Z_max (set by
# initialize.jl from the mesh) and maximum damping rate σ_max [1/τ₀].
if !@isdefined(sponge_zs_mhd)
    const sponge_zs_mhd = Ref{Float64}(30.0)
end
if !@isdefined(sponge_ztop_mhd)
    const sponge_ztop_mhd = Ref{Float64}(35.0)
end
if !@isdefined(sponge_sigma_mhd)
    const sponge_sigma_mhd = Ref{Float64}(2.0)
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

    ρ  = q[1]
    ρv = q[3]

    # Gravity, g = (0, -g₀): momentum and energy (ρ V·g = -ρ v g₀)
    S[3] = -ρ*g_mhd
    S[4] = -ρv*g_mhd

    # Well-balanced correction: minus the discrete residual of the vertical
    # balance of the initial state (initialize.jl, header item 5)
    if fe_well_balanced[]
        S[3] += fe_wb_lookup(y)
    end

    # GLM ψ damping: S_ψ = -(c_h²/c_p²) ψ = -(α_p c_h/Δh) ψ
    S[9] = -(glm_alpha_p_mhd[]*c_h_mhd[]/glm_dh_mhd[])*q[9]

    # Absorbing layer at the top of the domain
    zs = sponge_zs_mhd[]
    if y > zs
        zt = sponge_ztop_mhd[]
        s  = sin(0.5*π*min((y - zs)/(zt - zs), 1.0))
        σ  = sponge_sigma_mhd[]*s*s
        for ieq = 1:neqs
            S[ieq] -= σ*(q[ieq] - qe[ieq])
        end
    end

end
