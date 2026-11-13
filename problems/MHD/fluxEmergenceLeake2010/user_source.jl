#---------------------------------------------------------------------------------
# Sources of the flux-emergence problem of Leake, Linton & Antiochos (ApJ 722:550,
# 2010), in the code units of user_flux.jl (p = ρT, g = LK10_G, y up):
#
# 1. Gravity (Eq. 3):            S_ρv = −ρ g,   S_E = −ρv g.
#
# 2. Energy relaxation (Eqs. 5-6): the specific internal energy ε = p/((γ−1)ρ)
#    relaxes to "its initial equilibrium values ε₀(ρ)", the internal energy that
#    the initial atmosphere has at the height where its density equals ρ, on
#
#        τ = 23 s (ρ/ρ_ph)^(−1.43)  =  τ₀ ρ^(−1.43),   τ₀ = 23 s/t0 = 0.992,
#
#    i.e. S_E = −ρ(ε − ε₀(ρ))/τ = −(p − ρ T₀(ρ))/((γ−1)τ), with T₀(ρ) the
#    background temperature tabulated against ln ρ by initialize.jl. τ uses the
#    local density, as the authors' Lare runs do (Arber, Haynes & Leake 2007, Eq.
#    16: τ = 0.1 (ρ/ρ_ph)^(−1.7) after Gudiksen & Nordlund 2005): it is ~1 s in
#    the convection zone, 23 s at the photosphere, 2.6e3 s at 0.5 Mm and longer
#    than the run above ~1 Mm. Below ρ = lk10_relax_rhomin (10⁻⁵ ρ_ph, a rate
#    of 7e-8 per t0) the term is skipped. Applied everywhere, it also acts on the
#    buoyant tube in the convection zone: the tube is in thermal equilibrium at
#    a lower density than its surroundings, so ε₀(ρ) is that of a higher, cooler
#    layer and the tube cools within ~1 s, keeping a density deficit of ~21%
#    instead of 33% at its axis (tube 2). The paper does not say where the term
#    acts; :lk10_relax_ymin (default −∞) restricts it to y > y_min.
#
# 3. GLM ψ damping (Dedner et al. 2002):   S_ψ = −(α_p c_h/Δh) ψ.
#
# 4. Velocity damping near the side and top boundaries (not in the paper, which
#    does not state its boundary conditions; README, "Deviations"): as in the 3D
#    Lare runs of the same authors (Leake, Linton & Török 2013, Eq. 10), Dv/Dt
#    gets −N v with N = N_max (|x| − x_d)/(X − x_d) for |x| > x_d plus
#    N_max (y − y_d)/(Y − y_d) for y > y_d; x_d = 250 (37.5 Mm), y_d = 500
#    (75 Mm), N_max = 1/t0. In conserved form S_ρv = −N ρv, S_E = −N ρ|v|².
#
# The explicit viscosity and the anomalous resistivity of the paper are not
# sources: they are handed to the DynSGS kernel (:mhd_dissipation,
# :mhd_anomalous_eta in user_inputs.jl).
#---------------------------------------------------------------------------------
if !@isdefined(lk10_relax)
    const lk10_relax         = Ref{Bool}(true)
    const lk10_tau0          = Ref{Float64}(23.0/LK10_T0)
    const lk10_tau_exp       = Ref{Float64}(1.43)
    const lk10_relax_rhomin  = Ref{Float64}(1.0e-5)
    const lk10_lnrho_lo      = Ref{Float64}(0.0)          # T₀(ρ) table: ln ρ of its first entry,
    const lk10_dlnrho        = Ref{Float64}(1.0)          # its uniform step in ln ρ,
    const lk10_Ttab          = Ref{Vector{Float64}}(Float64[])   # and the background T there
    const lk10_glm_alpha     = Ref{Float64}(1.0)
    const lk10_glm_dh        = Ref{Float64}(1.0)
    const lk10_damp_xd       = Ref{Float64}(250.0)
    const lk10_damp_yd       = Ref{Float64}(500.0)
    const lk10_damp_rate     = Ref{Float64}(1.0)
    const lk10_xmax          = Ref{Float64}(300.0)
    const lk10_ymax          = Ref{Float64}(600.0)
end
if !@isdefined(lk10_relax_ymin)
    const lk10_relax_ymin    = Ref{Float64}(-Inf)         # relaxation only above this height (:lk10_relax_ymin)
end

# Background temperature T₀ at the height where the initial density equals ρ
@inline function lk10_T0_of_rho(ρ)
    tab = lk10_Ttab[]
    n   = length(tab)
    s   = (log(ρ) - lk10_lnrho_lo[])/lk10_dlnrho[]
    s <= 0.0     && return tab[1]
    s >= n - 1.0 && return tab[n]
    k = floor(Int, s)
    w = s - k
    return (1.0 - w)*tab[k+1] + w*tab[k+2]
end

@inline function lk10_damping_rate(x, y)
    N  = 0.0
    xd = lk10_damp_xd[]; yd = lk10_damp_yd[]
    ax = abs(x)
    ax > xd && (N += (ax - xd)/(lk10_xmax[] - xd))
    y  > yd && (N += (y - yd)/(lk10_ymax[] - yd))
    return lk10_damp_rate[]*N
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

    # gravity
    S[3] = -ρ*LK10_G
    S[4] = -q[3]*LK10_G

    # energy relaxation to ε₀(ρ)
    if lk10_relax[] && ρ > lk10_relax_rhomin[] && y > lk10_relax_ymin[]
        p    = pressure_mhd(ρ, q[2], q[3], q[5], q[4], q[6], q[7], q[8], q[9])
        rate = ρ^lk10_tau_exp[]/lk10_tau0[]
        S[4] -= rate*(p - ρ*lk10_T0_of_rho(ρ))/(γ_mhd - 1.0)
    end

    # GLM damping
    S[9] = -(lk10_glm_alpha[]*c_h_mhd[]/lk10_glm_dh[])*q[9]

    # velocity damping near the side and top boundaries
    N = lk10_damping_rate(x, y)
    if N > 0.0
        S[2] -= N*q[2]
        S[3] -= N*q[3]
        S[5] -= N*q[5]
        S[4] -= N*(q[2]*q[2] + q[3]*q[3] + q[5]*q[5])/ρ
    end
end
