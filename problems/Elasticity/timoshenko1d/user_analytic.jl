#---------------------------------------------------------------------------------
# Exact standing mode of the simply supported Timoshenko beam, overlaid on the
# 1D output plots.
#
# The shared 1D plotter (src/io/plotting/jeplots.jl) calls
#
#     user_analytic_solution(x, t, outvar, inputs) -> npoin × nvar matrix
#
# if — and only if — the case being run ships this file. Entries that cannot be
# supplied should be NaN: Plots skips non-finite points.
#
# With k = nπ/L, and Ω the flexural root of the Timoshenko frequency equation
# (see timoshenko_mode() in user_flux.jl),
#
#     w = W sin(kx) cos(Ωt)              φ = Φ cos(kx) cos(Ωt)
#     v = ẇ = -W Ω sin(kx) sin(Ωt)       ω = φ̇ = -Φ Ω cos(kx) sin(Ωt)
#     γ = w_x - φ = (Wk - Φ) cos(kx) cos(Ωt)
#     χ = φ_x     = -Φ k sin(kx) cos(Ωt)
#
# This is an exact solution of the continuous system, not a series truncation,
# so the gap between the two curves is entirely the discretisation's.
#---------------------------------------------------------------------------------
function user_analytic_solution(x, t, outvar, inputs)

    p = timoshenko_properties()
    m = timoshenko_mode()

    nvar = length(outvar)
    q    = fill(NaN, length(x), nvar)

    cosΩt = cos(m.Ω*t)
    sinΩt = sin(m.Ω*t)

    for i in eachindex(x)
        kx = m.k*x[i]

        w = m.W*sin(kx)*cosΩt
        φ = m.Φ*cos(kx)*cosΩt
        v = -m.W*m.Ω*sin(kx)*sinΩt
        ω = -m.Φ*m.Ω*cos(kx)*sinΩt
        γ = (m.W*m.k - m.Φ)*cos(kx)*cosΩt
        χ = -m.Φ*m.k*sin(kx)*cosΩt

        for ivar = 1:nvar
            name = string(outvar[ivar])
            q[i, ivar] = name == "ρAv" ? p.ρA*v :
                         name == "ρIω" ? p.ρI*ω :
                         name == "γ"   ? γ      :
                         name == "χ"   ? χ      :
                         name == "v"   ? v      :
                         name == "ω"   ? ω      :
                         name == "w"   ? w      :
                         name == "φ"   ? φ      :
                         name == "Q"   ? p.κGA*γ :
                         name == "M"   ? p.EI*χ  : NaN
        end
    end

    return q
end
