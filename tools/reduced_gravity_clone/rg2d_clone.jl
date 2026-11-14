#!/usr/bin/env julia
#
# Standalone 2D CG-SEM clone of
# problems/ShallowWater/ReducedGravity (Jexpresso). See README.md.
#
# usage: julia -O3 rg2d_clone.jl key=value ...
#   nel=25 N=4 dt=600 tend_days=2160 filter=1 mu_filter=0.05 limiter=1 hmin=1e-3
#   mu=300 visc=cons|vel wind=smoothstep|linear|off hwind=0.1 hwind_lin=1.0
#   sigma=2.0833e-4 sigma_mode=smooth|hard report_days=5 umax_stop=1e3 tag=run dump=file
#
using LinearAlgebra, Printf

# ----------------------------------------------------------------------------
function lgl(N)
    N1 = N + 1
    x = cos.(pi .* (0:N) ./ N)
    P = zeros(N1, N1)
    xold = fill(2.0, N1)
    while maximum(abs.(x .- xold)) > 10eps(Float64)
        xold .= x
        P[:, 1] .= 1.0
        P[:, 2] .= x
        for k = 2:N
            P[:, k+1] .= ((2k - 1) .* x .* P[:, k] .- (k - 1) .* P[:, k-1]) ./ k
        end
        x .= xold .- (x .* P[:, N1] .- P[:, N]) ./ (N1 .* P[:, N1])
    end
    w = 2 ./ (N * N1 .* P[:, N1] .^ 2)
    idx = sortperm(x)
    return x[idx], w[idx]
end

function diffmat(x)
    n = length(x)
    wb = [1 / prod(x[j] - x[k] for k in 1:n if k != j) for j in 1:n]
    D = zeros(n, n)
    for i in 1:n, j in 1:n
        i != j && (D[i, j] = wb[j] / wb[i] / (x[i] - x[j]))
    end
    for i in 1:n
        D[i, i] = -sum(D[i, j] for j in 1:n if j != i)
    end
    return D
end

function vandeven_modal(kk, ngl, p)   # verbatim port of filter.jl
    pe = 0.3275911; a1 = 0.254829592; a2 = -0.284496736; a3 = 1.421413741; a4 = -1.453152027; a5 = 1.061405429
    n = ngl - 1; k = kk - 1; i = 2 * n / 3; ep = 1.0e-10
    if k <= i
        return 1.0
    elseif k > i && k < n
        x = Float64(k - i) / Float64(n - i)
        omega = abs(x) - 0.5
        xlog = log(1.0 - 4.0 * omega^2)
        c = 4.0 * omega^2
        diff = abs(x - 0.5)
        square_root = diff < ep ? 1.0 : sqrt(-xlog / c)
        z = 2.0 * sqrt(p) * omega * square_root
        zc = abs(z)
        t = 1.0 / (1.0 + pe * zc)
        c = 1.0 - (a1 * t + a2 * t^2 + a3 * t^3 + a4 * t^4 + a5 * t^5) * exp(-zc * zc)
        c = zc < ep ? 0.0 : c * z / zc
        return 0.5 * (1.0 - c)
    else
        return 0.0
    end
end

function erf_filter_matrix(N, xgl, mu)
    ngl = N + 1
    leg = zeros(ngl, ngl)
    for i in 1:ngl
        leg[i, 1] = 1.0
        ngl > 1 && (leg[i, 2] = xgl[i])
        for j in 3:ngl
            k = j - 2
            leg[i, j] = ((2k + 1) * xgl[i] * leg[i, j-1] - k * leg[i, j-2]) / (k + 1)
        end
    end
    leg2 = copy(leg)
    for i in 1:ngl
        leg2[i, 1] = 0.5 * (1 - xgl[i]); leg2[i, 2] = 0.5 * (1 + xgl[i])
        for j in 3:ngl
            leg2[i, j] = leg[i, j] - leg[i, j-2]
        end
    end
    leginv = inv(leg2)
    weight = [vandeven_modal(k, ngl, 12) for k in 1:ngl]
    f = mu .* (leg2 * Diagonal(weight) * leginv) .+ (1 - mu) .* Matrix(I, ngl, ngl)
    return f, weight
end

# ----------------------------------------------------------------------------
Base.@kwdef mutable struct Cfg
    nel::Int = 25; N::Int = 4; dt::Float64 = 600.0; tend_days::Float64 = 6 * 360.0
    filter::Bool = true; mu_filter::Float64 = 0.05
    limiter::Bool = true; hmin::Float64 = 1.0e-3
    mu::Float64 = 300.0; visc::String = "cons"; mu_h::Float64 = 0.0
    wind::String = "smoothstep"; hwind::Float64 = 0.1; hwind_lin::Float64 = 1.0; hrelax::Float64 = 0.1
    sigma::Float64 = 0.25 / 1200; sigma_mode::String = "smooth"
    report_days::Float64 = 5.0; umax_stop::Float64 = 1.0e3; tag::String = "run"; dump::String = ""
    dump_every_days::Float64 = 0.0; dump_after_days::Float64 = 0.0; report_after_days::Float64 = 0.0; report_fine_days::Float64 = 0.0
    check_every_step::Bool = false; ring_len::Int = 0; ring_every::Int = 6
    ckpt_days::String = ""; restart::String = ""
    eps_u::Float64 = 1.0e-3      # desingularization scale of u = Hu/H (deck: = _H_WET_SWE = 1 mm)
    ucap::Float64 = 0.0          # stage velocity cap |Hu| <= ucap*H (0 = off)
    press::String = "cons"       # "cons": flux g'(H²-He²)/2 (deck); "noncons": g' H ∇H
    diag_umax::Float64 = 0.0     # print the momentum budget at the max-|u| node when max|u| > diag_umax (0 = off)
    diag_max::Int = 6            # at most this many budgets
    ent_He::Float64 = 0.0        # entrainment thickness H_e [m] (McCreary & Kundu 1988); 0 = off
    ent_te::Float64 = 86400.0    # entrainment time scale t_e [s]
    ent_comp::Bool = false       # remove the entrained volume uniformly (volume-conserving)
    # physics (ReducedGravity/user_flux.jl, user_source.jl)
    g::Float64 = 0.015; Hwet::Float64 = 1.0e-3; He::Float64 = 150.0
    f0::Float64 = 8.37e-5; beta::Float64 = 1.8e-11; tau0::Float64 = 0.03; rho::Float64 = 1000.0; kappa::Float64 = 1.0e-8
    L::Float64 = 2.0e6
end

function parse_cfg(args)
    c = Cfg()
    for a in args
        k, v = split(a, "=")
        s = Symbol(k)
        T = typeof(getfield(c, s))
        val = T == String ? String(v) : T == Bool ? (v in ("1", "true")) : parse(T, v)
        setfield!(c, s, val)
    end
    return c
end

@inline function smoothstep(H, Hlo, Hhi)
    s = clamp((H - Hlo) / (Hhi - Hlo), 0.0, 1.0)
    return s * s * (3.0 - 2.0 * s)
end
@inline wet_fraction(H, c::Cfg) = smoothstep(H, c.Hwet, c.hwind)        # wind ramp (deck: 1 mm -> 0.1 m)
@inline relax_fraction(H, c::Cfg) = smoothstep(H, c.Hwet, c.hrelax)     # thin-film relaxation ramp
@inline function uvel(Hc, Hm, eps)
    H4 = max(Hc, eps)
    return sqrt(2.0) * Hc * Hm / sqrt(Hc^4 + H4^4)
end

struct Solver
    c::Cfg
    ngl::Int; nel::Int; npoin::Int
    xi::Vector{Float64}; w::Vector{Float64}; D::Matrix{Float64}; F::Matrix{Float64}; fweights::Vector{Float64}
    dxidx::Float64; detady::Float64; J::Float64; wJ::Matrix{Float64}
    conn::Array{Int,3}; x::Vector{Float64}; y::Vector{Float64}; M::Vector{Float64}
    left::Vector{Int}; right::Vector{Int}; bottom::Vector{Int}; top::Vector{Int}
    dxmin::Float64
    # scratch
    Fx::Matrix{Float64}; Gy::Matrix{Float64}; Src::Matrix{Float64}; RHS::Matrix{Float64}; B::Matrix{Float64}
    prim::Matrix{Float64}; U::Vector{Float64}; V::Vector{Float64}; Hc::Vector{Float64}
    ql::Array{Float64,3}; tmp::Array{Float64,3}; fx::Array{Float64,3}; fy::Array{Float64,3}
end

function Solver(c::Cfg)
    N = c.N; ngl = N + 1; nelx = nely = c.nel; L = c.L
    xi, w = lgl(N)
    D = diffmat(xi)
    F, fweights = erf_filter_matrix(N, xi, c.mu_filter)
    hx = hy = L / nelx
    dxidx = 2 / hx; detady = 2 / hy; J = hx * hy / 4
    nxn = nelx * N + 1; nyn = nely * N + 1
    npoin = nxn * nyn
    nel = nelx * nely
    conn = zeros(Int, nel, ngl, ngl)
    x = zeros(npoin); y = zeros(npoin)
    for ey in 0:nely-1, ex in 0:nelx-1
        e = ey * nelx + ex + 1
        for j in 1:ngl, i in 1:ngl
            ip = (ey * N + j - 1) * nxn + (ex * N + i - 1) + 1
            conn[e, i, j] = ip
            x[ip] = (ex + 0.5 * (xi[i] + 1)) * hx
            y[ip] = (ey + 0.5 * (xi[j] + 1)) * hy
        end
    end
    wJ = [w[i] * w[j] * J for i in 1:ngl, j in 1:ngl]
    M = zeros(npoin)
    for e in 1:nel, j in 1:ngl, i in 1:ngl
        M[conn[e, i, j]] += wJ[i, j]
    end
    tol = 1e-6 * hx
    left = findall(abs.(x) .< tol); right = findall(abs.(x .- L) .< tol)
    bottom = findall(abs.(y) .< tol); top = findall(abs.(y .- L) .< tol)
    dxmin = (xi[2] - xi[1]) / 2 * hx
    Solver(c, ngl, nel, npoin, xi, w, D, F, fweights, dxidx, detady, J, wJ, conn, x, y, M,
           left, right, bottom, top, dxmin,
           zeros(3, npoin), zeros(3, npoin), zeros(3, npoin), zeros(3, npoin), zeros(3, npoin),
           zeros(3, npoin), zeros(npoin), zeros(npoin), zeros(npoin),
           zeros(3, ngl, ngl), zeros(3, ngl, ngl), zeros(3, ngl, ngl), zeros(3, ngl, ngl))
end

@inline qe_of(m, c::Cfg) = m == 1 ? c.He : 0.0

function apply_filter!(q, S::Solver)
    c = S.c; ngl = S.ngl; conn = S.conn; Fm = S.F; ql = S.ql; tmp = S.tmp; B = S.B
    fill!(B, 0.0)
    @inbounds for e in 1:S.nel
        for j in 1:ngl, i in 1:ngl, m in 1:3
            ql[m, i, j] = q[m, conn[e, i, j]] - qe_of(m, c)
        end
        for j in 1:ngl, i in 1:ngl, m in 1:3
            acc = 0.0
            for k in 1:ngl
                acc += Fm[i, k] * ql[m, k, j]
            end
            tmp[m, i, j] = acc
        end
        for j in 1:ngl, i in 1:ngl, m in 1:3
            acc = 0.0
            for l in 1:ngl
                acc += tmp[m, i, l] * Fm[j, l]
            end
            B[m, conn[e, i, j]] += acc * S.wJ[i, j]
        end
    end
    @inbounds for ip in 1:S.npoin, m in 1:3
        q[m, ip] = B[m, ip] / S.M[ip] + qe_of(m, c)
    end
end

function apply_wall_bc!(q, S::Solver)
    @inbounds for ip in S.left; q[2, ip] = 0.0; end
    @inbounds for ip in S.right; q[2, ip] = 0.0; end
    @inbounds for ip in S.bottom; q[3, ip] = 0.0; end
    @inbounds for ip in S.top; q[3, ip] = 0.0; end
end

# RHS_MODE: 0 = full RHS (normal); for the term breakdown: 1 = pressure flux divergence only,
# 2 = advective flux divergence only, 3 = viscous term only (no filter/BC applied in modes 1-3)
const RHS_MODE = Ref(0)

function rhs!(dq, q, S::Solver, t)
    c = S.c
    mode = RHS_MODE[]
    if mode == 0
        c.filter && apply_filter!(q, S)
        apply_wall_bc!(q, S)
    end
    noncons = c.press == "noncons"
    pf = (mode == 0 || mode == 1) && !noncons ? 1.0 : 0.0   # conservative pressure flux
    af = (mode == 0 || mode == 2) ? 1.0 : 0.0                # advective fluxes
    sf = mode == 0 ? 1.0 : 0.0                               # pointwise sources
    lvisc = mode == 0 || mode == 3
    ngl = S.ngl; conn = S.conn; D = S.D; RHS = S.RHS; Fx = S.Fx; Gy = S.Gy; Src = S.Src
    g = c.g; He = c.He; eps = c.eps_u; L = c.L
    @inbounds for ip in 1:S.npoin
        H = q[1, ip]; Hu = q[2, ip]; Hv = q[3, ip]
        Hc = max(H, 0.0)
        S.Hc[ip] = Hc
        u = uvel(Hc, Hu, eps); v = uvel(Hc, Hv, eps)
        p = 0.5 * g * (Hc * Hc - He * He)
        Fx[1, ip] = af * Hc * u; Fx[2, ip] = af * Hc * u * u + pf * p; Fx[3, ip] = af * Hc * u * v
        Gy[1, ip] = af * Hc * v; Gy[2, ip] = af * Hc * v * u; Gy[3, ip] = af * Hc * v * v + pf * p
        S.U[ip] = u; S.V[ip] = v
        # sources (user_source.jl)
        y = S.y[ip]
        f = c.f0 + c.beta * (y - 0.5 * L)
        taux = -c.tau0 * cospi(2.0 * y / L)
        wet = c.wind == "smoothstep" ? wet_fraction(H, c) :
              c.wind == "linear" ? clamp(Hc / c.hwind_lin, 0.0, 1.0) : 0.0
        sigma = c.sigma_mode == "hard" ? (H < c.Hwet ? c.sigma : 0.0) : (1.0 - relax_fraction(H, c)) * c.sigma
        Src[1, ip] = (sf > 0 && c.ent_He > 0 && H < c.ent_He) ? (c.ent_He - H)^2 / (c.ent_te * c.ent_He) : 0.0
        Src[2, ip] = sf * (f * Hv + wet * taux / c.rho - (c.kappa + sigma) * Hu)
        Src[3, ip] = sf * (-f * Hu - (c.kappa + sigma) * Hv)
        S.prim[1, ip] = H - He
        if c.visc == "cons"
            S.prim[2, ip] = Hu; S.prim[3, ip] = Hv
        else
            S.prim[2, ip] = u; S.prim[3, ip] = v
        end
    end
    if c.ent_comp && c.ent_He > 0 && mode == 0
        Etot = 0.0; Atot = 0.0
        @inbounds for ip in 1:S.npoin
            Etot += S.M[ip] * Src[1, ip]; Atot += S.M[ip]
        end
        wbar = Etot / Atot
        @inbounds for ip in 1:S.npoin
            Src[1, ip] -= wbar
        end
    end
    fill!(RHS, 0.0)
    wJ = S.wJ; dxidx = S.dxidx; detady = S.detady
    @inbounds for e in 1:S.nel
        for j in 1:ngl, i in 1:ngl
            ip = conn[e, i, j]
            for m in 1:3
                dF = 0.0; dG = 0.0
                for k in 1:ngl
                    dF += D[i, k] * Fx[m, conn[e, k, j]]
                    dG += D[j, k] * Gy[m, conn[e, i, k]]
                end
                RHS[m, ip] -= wJ[i, j] * (dxidx * dF + detady * dG - Src[m, ip])
            end
            if noncons && (mode == 0 || mode == 1)
                # wet-aware pressure gradient g' H_i (∇H)_i: vanishes on a dry node whatever
                # its neighbours, unlike the strong-form derivative of the nodal flux g'H²/2
                dHx = 0.0; dHy = 0.0
                for k in 1:ngl
                    dHx += D[i, k] * S.Hc[conn[e, k, j]]
                    dHy += D[j, k] * S.Hc[conn[e, i, k]]
                end
                RHS[2, ip] -= wJ[i, j] * g * S.Hc[ip] * dxidx * dHx
                RHS[3, ip] -= wJ[i, j] * g * S.Hc[ip] * detady * dHy
            end
        end
    end
    if lvisc && (c.mu > 0.0 || c.mu_h > 0.0)
        fx = S.fx; fy = S.fy; prim = S.prim
        mlo = c.mu_h > 0.0 ? 1 : 2           # slot 1 (H - He) is diffused only when mu_h > 0
        @inbounds for e in 1:S.nel
            for l in 1:ngl, k in 1:ngl
                ipkl = conn[e, k, l]
                coefm = c.mu * (c.visc == "cons" ? 1.0 : max(q[1, ipkl], 0.0))
                for m in mlo:3
                    coef = m == 1 ? c.mu_h : coefm
                    gx = 0.0; gy = 0.0
                    for ii in 1:ngl
                        gx += D[k, ii] * prim[m, conn[e, ii, l]]
                        gy += D[l, ii] * prim[m, conn[e, k, ii]]
                    end
                    fx[m, k, l] = coef * gx * dxidx * dxidx * wJ[k, l]
                    fy[m, k, l] = coef * gy * detady * detady * wJ[k, l]
                end
            end
            for l in 1:ngl, i in 1:ngl, m in mlo:3
                acc = 0.0
                for k in 1:ngl
                    acc += D[k, i] * fx[m, k, l]
                end
                RHS[m, conn[e, i, l]] -= acc
            end
            for k in 1:ngl, i in 1:ngl, m in mlo:3
                acc = 0.0
                for l in 1:ngl
                    acc += D[l, i] * fy[m, k, l]
                end
                RHS[m, conn[e, k, i]] -= acc
            end
        end
    end
    @inbounds for ip in 1:S.npoin, m in 1:3
        dq[m, ip] = RHS[m, ip] / S.M[ip]
    end
    return nothing
end

function clip!(q, S::Solver)
    c = S.c
    n = 0
    if c.limiter
        @inbounds for ip in 1:S.npoin
            if q[1, ip] < c.hmin
                q[2, ip] = 0.0; q[3, ip] = 0.0; q[1, ip] = c.hmin
                n += 1
            end
        end
    end
    if c.ucap > 0.0
        # velocity cap: |Hu| <= ucap * max(H, 0), mass untouched
        @inbounds for ip in 1:S.npoin
            Hc = max(q[1, ip], 0.0)
            mm = hypot(q[2, ip], q[3, ip])
            lim = c.ucap * Hc
            if mm > lim
                r = mm > 0.0 ? lim / mm : 0.0
                q[2, ip] *= r; q[3, ip] *= r
            end
        end
    end
    return n
end

# momentum budget at node ip: the contribution of each term to d(Hu)/dt, d(Hv)/dt [m²/s²]
function breakdown(q, S::Solver, ip, t)
    c = S.c
    qq = copy(q)
    c.filter && apply_filter!(qq, S)
    # the filter is applied in place at every RK stage (5 per step): its increment
    # per application at this node, and the equivalent rate 5*Δq/Δt
    dqf = qq[:, ip] .- q[:, ip]
    apply_wall_bc!(qq, S)
    out = zeros(3, S.npoin)
    terms = Dict{String,Vector{Float64}}()
    for (m, name) in ((1, "pressure"), (2, "advection"), (3, "viscous"))
        RHS_MODE[] = m
        rhs!(out, qq, S, t)
        terms[name] = out[:, ip]
    end
    RHS_MODE[] = 0
    H = qq[1, ip]; Hu = qq[2, ip]; Hv = qq[3, ip]; Hc = max(H, 0.0)
    y = S.y[ip]; L = c.L
    f = c.f0 + c.beta * (y - 0.5 * L)
    taux = -c.tau0 * cospi(2.0 * y / L)
    wet = c.wind == "smoothstep" ? wet_fraction(H, c) : c.wind == "linear" ? clamp(Hc / c.hwind_lin, 0.0, 1.0) : 0.0
    sigma = c.sigma_mode == "hard" ? (H < c.Hwet ? c.sigma : 0.0) : (1.0 - relax_fraction(H, c)) * c.sigma
    terms["coriolis"] = [0.0, f * Hv, -f * Hu]
    terms["wind"] = [0.0, wet * taux / c.rho, 0.0]
    terms["relax+drag"] = [0.0, -(c.kappa + sigma) * Hu, -(c.kappa + sigma) * Hv]
    nxn = c.nel * c.N + 1
    @printf("#   BUDGET t=%.3f d node %d (x,y)=(%.1f,%.1f) km: H=%.4e Hu=%+.4e Hv=%+.4e u=%+.3e v=%+.3e (after filter/BC)\n",
            t / 86400, ip, S.x[ip] / 1e3, S.y[ip] / 1e3, H, Hu, Hv, Hc > 0 ? Hu / Hc : 0.0, Hc > 0 ? Hv / Hc : 0.0)
    ix = (ip - 1) % nxn; iy = (ip - 1) ÷ nxn
    print("#     H in the 5x5 patch around it (rows: y from +2 to -2 nodes; cols: x from -2 to +2):\n")
    for dy in 2:-1:-2
        print("#      ")
        for dx in -2:2
            jx = ix + dx; jy = iy + dy
            if 0 <= jx < nxn && 0 <= jy < nxn
                @printf(" %10.3e", qq[1, jy * nxn + jx + 1])
            else
                print("          -")
            end
        end
        print("\n")
    end
    for name in ("pressure", "advection", "viscous", "coriolis", "wind", "relax+drag")
        @printf("#     %-11s d(Hu)/dt=%+.4e  d(Hv)/dt=%+.4e   dH/dt=%+.4e\n", name, terms[name][2], terms[name][3], terms[name][1])
    end
    tot = sum(values(terms))
    @printf("#     %-11s d(Hu)/dt=%+.4e  d(Hv)/dt=%+.4e   dH/dt=%+.4e\n", "TOTAL", tot[2], tot[3], tot[1])
    rf = 5.0 / c.dt
    @printf("#     %-11s d(Hu)/dt~%+.4e  d(Hv)/dt~%+.4e   dH/dt~%+.4e   (5 x per-application increment / dt; increment: %+.3e %+.3e %+.3e)\n",
            "filter", rf * dqf[2], rf * dqf[3], rf * dqf[1], dqf[2], dqf[3], dqf[1])
    flush(stdout)
end

# SSPRK54 (Spiteri & Ruuth) with the stage limiter, as OrdinaryDiffEq does it
function step!(u0, S::Solver, t, dt, k, k3, k4, u1, u2, u3, u4, u5)
    rhs!(k, u0, S, t)
    @. u1 = u0 + 0.391752226571890 * dt * k; clip!(u1, S)
    rhs!(k, u1, S, t + 0.391752226571890 * dt)
    @. u2 = 0.444370493651235 * u0 + 0.555629506348765 * u1 + 0.368410593050371 * dt * k; clip!(u2, S)
    rhs!(k, u2, S, t + 0.586079689311540 * dt)
    @. u3 = 0.620101851488403 * u0 + 0.379898148511597 * u2 + 0.251891774271694 * dt * k; clip!(u3, S)
    rhs!(k3, u3, S, t + 0.474542363121400 * dt)
    @. u4 = 0.178079954393132 * u0 + 0.821920045606868 * u3 + 0.544974750228521 * dt * k3; clip!(u4, S)
    rhs!(k4, u4, S, t + 0.935010630967653 * dt)
    @. u5 = 0.517231671970585 * u2 + 0.096059710526147 * u3 + 0.063692468666290 * dt * k3 +
            0.386708617503269 * u4 + 0.226007483236906 * dt * k4
    return clip!(u5, S)
end

function dump_state(fname, S::Solver, q, t)
    open(fname, "w") do io
        @printf(io, "# t = %.6e s\n# x y H Hu Hv u v\n", t)
        for ip in 1:S.npoin
            Hc = max(q[1, ip], 0.0)
            @printf(io, "%.6e %.6e %.10e %.10e %.10e %.10e %.10e\n", S.x[ip], S.y[ip], q[1, ip], q[2, ip], q[3, ip],
                    uvel(Hc, q[2, ip], S.c.eps_u), uvel(Hc, q[3, ip], S.c.eps_u))
        end
    end
end

function main()
    c = parse_cfg(ARGS)
    S = Solver(c)
    @printf("# %s: nel=%dx%d N=%d npoin=%d dx_min=%.2f km dt=%gs filter=%d(mu_f=%g, weights=%s) limiter=%d(hmin=%g) ucap=%g eps_u=%g press=%s mu=%g(%s) mu_h=%g wind=%s(hwind=%g, hwind_lin=%g) sigma=%.3e(%s, hrelax=%g) entrainment(He=%g, te=%g s, comp=%d)\n",
            c.tag, c.nel, c.nel, c.N, S.npoin, S.dxmin / 1e3, c.dt, c.filter, c.mu_filter, string(round.(S.fweights, digits=4)),
            c.limiter, c.hmin, c.ucap, c.eps_u, c.press, c.mu, c.visc, c.mu_h, c.wind, c.hwind, c.hwind_lin, c.sigma, c.sigma_mode, c.hrelax, c.ent_He, c.ent_te, c.ent_comp)
    println("#   t[days]   min(H)      x,y of min(H) [km]   max|u| [m/s] (H there)   max|Hu|   #H<1m  #clips   mass/mass0-1")
    flush(stdout)
    q = zeros(3, S.npoin); q[1, :] .= c.He
    bufs = [zeros(3, S.npoin) for _ in 1:8]
    day = 86400.0
    t = 0.0
    if c.restart != ""
        # binary checkpoint written by this script (ckpt_days=...): t, then q (3 x npoin, Float64)
        open(c.restart, "r") do io
            t = read(io, Float64)
            read!(io, q)
        end
        @printf("# restarted from %s at t = %.4f days\n", c.restart, t / day)
    end
    qprev = copy(q)
    # ring buffer of recent states (every ring_every steps), written out at the blow-up
    ring = [zeros(3, S.npoin) for _ in 1:max(c.ring_len, 1)]
    ring_t = fill(NaN, max(c.ring_len, 1)); ring_k = 0
    ckpt = c.ckpt_days == "" ? Float64[] : parse.(Float64, split(c.ckpt_days, ","))
    nsteps = round(Int, (c.tend_days * day - t) / c.dt)
    report_every = max(1, round(Int, c.report_days * day / c.dt))
    report_fine = c.report_fine_days > 0 ? max(1, round(Int, c.report_fine_days * day / c.dt)) : 0
    dump_every = c.dump_every_days > 0 ? max(1, round(Int, c.dump_every_days * day / c.dt)) : 0
    mass0 = sum(S.M .* q[1, :])
    nclip_tot = 0; t0 = time(); n = 0; ndiag = 0
    blown = false
    for n in 1:nsteps
        qprev .= q
        nclip = step!(q, S, t, c.dt, bufs...)
        # u5 is the last buffer
        q, bufs[8] = bufs[8], q
        nclip_tot += nclip
        t += c.dt
        bad = !all(isfinite, q)
        for td in ckpt
            if c.dump != "" && t - c.dt < td * day <= t
                open(@sprintf("%s_ckpt_d%07.1f.bin", c.dump, td), "w") do io
                    write(io, t); write(io, q)
                end
            end
        end
        if c.ring_len > 0 && t >= c.report_after_days * day && n % c.ring_every == 0 && !bad
            ring_k = ring_k % c.ring_len + 1
            ring[ring_k] .= q; ring_t[ring_k] = t
        end
        if dump_every > 0 && c.dump != "" && t >= c.dump_after_days * day && n % dump_every == 0 && !bad
            dump_state(@sprintf("%s_t%08.2f.dat", c.dump, t / day), S, q, t)
        end
        fine = report_fine > 0 && t >= c.report_after_days * day && n % report_fine == 0
        # cheap per-step speed check late in the run, so that the stop (and the _prev dump) lands
        # on the first step where max|u| exceeds umax_stop rather than on the next report
        if !fine && !bad && c.check_every_step && t >= c.report_after_days * day
            umax_chk = 0.0
            @inbounds for ip in 1:S.npoin
                Hc = max(q[1, ip], 0.0)
                umax_chk = max(umax_chk, hypot(uvel(Hc, q[2, ip], c.eps_u), uvel(Hc, q[3, ip], c.eps_u)))
            end
            fine = umax_chk > c.umax_stop
        end
        if n % report_every == 0 || fine || bad || n <= 2
            H = view(q, 1, :)
            umax = 0.0; iu = 1; im = 1; hmin = Inf; nthin = 0
            if !bad
                for ip in 1:S.npoin
                    Hc = max(H[ip], 0.0)
                    spd = hypot(uvel(Hc, q[2, ip], c.eps_u), uvel(Hc, q[3, ip], c.eps_u))
                    if spd > umax; umax = spd; iu = ip; end
                    if H[ip] < hmin; hmin = H[ip]; im = ip; end
                    H[ip] < 1.0 && (nthin += 1)
                end
            end
            mass = sum(S.M .* H)
            @printf("%9.2f %11.4e  (%7.1f,%7.1f)  %9.3e (%9.3e)  %9.3e  %5d  %7d  %+.3e\n",
                    t / day, hmin, S.x[im] / 1e3, S.y[im] / 1e3, umax, H[iu], maximum(abs, view(q, 2, :)), nthin, nclip_tot, mass / mass0 - 1)
            flush(stdout)
            if !bad && c.diag_umax > 0 && umax > c.diag_umax && ndiag < c.diag_max
                ndiag += 1
                breakdown(q, S, iu, t)
            end
            if bad
                @printf("# BLOW-UP (non-finite state) at t = %.2f days, step %d\n", t / day, n)
                blown = true
            elseif umax > c.umax_stop
                @printf("# BLOW-UP (|u| > %g m/s) at t = %.2f days, step %d: |u|=%.3e at H=%.3e, (x,y)=(%.1f,%.1f) km\n",
                        c.umax_stop, t / day, n, umax, H[iu], S.x[iu] / 1e3, S.y[iu] / 1e3)
                blown = true
            end
            if blown
                c.dump != "" && dump_state(c.dump * "_prev.dat", S, qprev, t - c.dt)
                c.dump != "" && dump_state(c.dump * "_last.dat", S, q, t)
                if c.dump != "" && c.ring_len > 0
                    for k in 1:c.ring_len
                        isnan(ring_t[k]) && continue
                        dump_state(@sprintf("%s_ring_t%010.4f.dat", c.dump, ring_t[k] / day), S, ring[k], ring_t[k])
                    end
                end
                break
            end
        end
    end
    if !blown
        @printf("# completed %d steps to t = %.1f days without blow-up; min(H) = %.4e\n", nsteps, t / day, minimum(view(q, 1, :)))
        c.dump != "" && dump_state(c.dump * "_last.dat", S, q, t)
    end
    @printf("# wall time %.1f s\n", time() - t0)
end

main()
