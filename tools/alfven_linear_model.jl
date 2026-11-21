#!/usr/bin/env julia
#---------------------------------------------------------------------------------
# THE ALFVÉN WAVE ON JEXPRESSO'S CG OPERATOR, LINEARIZED: a model that runs in
# seconds to minutes and reproduces the code's error on
# problems/MHD/alfvenWaveBabbarRanocha2026 in the linear regime to 3-4 digits
# (N = 4: in-plane error 4.02x and 12.87x the out-of-plane one on 16² and 32²;
# the code at amplitude 0.01: 4.01x, 12.87x).
#
# Around the constant background q0 (ρ = 1, u = 0, p = 0.1, B = e∥, ψ = 0) the
# collocated CG scheme dq/dt = −M⁻¹ Σ_e ∫ φ (∂x F(q_h) + ∂y G(q_h)) + S(q)
# linearizes EXACTLY to
#
#     dq'/dt = −(A_x ⊗ D_x + A_y ⊗ D_y) q' + (A_s ⊗ I) q',
#
# with A_x, A_y, A_s the Jacobians of the DECK'S OWN flux and source at q0
# (user_flux.jl and user_source.jl of the case, included unchanged; complex-step
# derivatives, exact to round-off for this rational flux) and D_x, D_y the
# assembled CG derivatives (tensor-product LGL, lumped mass, periodic [0, √2]²).
# Integrated with Vern9 at reltol 1e-12, so the time error is out of the way.
# The exact linear solution is the deck's wave with A → 0:
#     u⊥ = b⊥ = sin(2π(x∥ + t)),  w = b_z = cos(2π(x∥ + t)).
#
# Reported (RMS over the box at N + 4 Gauss points per direction, as the case):
#     b·e⊥   the in-plane wave component
#     b·e∥   the in-plane component along B0, exactly 0: what the discrete
#            divergence of the in-plane fields drives
#     Bz     the out-of-plane field: its equations are scalar advection along
#            e∥, so this is the plain CG-SEM error, the reference
#
# Switches (first argument): "full", or any of these joined by "_"
# (powell_noglm, noglm_nogasp_nomagp, ...). Each removes one coupling of the
# in-plane wave to the rest of the linear system, at the Jacobian level:
#     noglm      c_h = 0: no GLM cleaning (ψ stays 0)
#     powell     + the Godunov-Powell term −B0 (∇·b) on the momentum, which
#                cancels the spurious divergence force of the conservative form
#     nogasp     no gas-pressure force: ∂p/∂q (the deck's pressure_mhd) taken
#                out of the x-flux of ρu and the y-flux of ρv
#     nomagp     no magnetic-pressure force: ∂(½|B|²)/∂q taken out of the same
#     noinddiv   no field compression in the induction: + B0 (∇·u) cancels the
#                −B0 (∇·u) of the flux, the term that turns a discrete ∇·u into b∥
# With all five, u_t = (B0·D) b and b_t = (B0·D) u: the in-plane pair obeys the
# out-of-plane equations, so ⊥/Bz must come out ~1 and ∥/Bz ~0 (a check of the
# switches themselves).
#
#   julia --project=. tools/alfven_linear_model.jl full 4,5,6 8,16
#   julia --project=. tools/alfven_linear_model.jl noglm 4 8,16,32
#
# Options:
#   --t=T          final time (default 2, the deck's: one crossing of the box along e∥)
#   --field[=DIR]  where the error is, per (N, nel) (DIR default alfven_linear_fields):
#                  printed, the RMS error of b·e⊥, b·e∥, Bz, ψ, p', u·e⊥ split into the
#                  wave's own Fourier mode ±(1,1), the element-scale rings
#                  ±(1,1) + nel·(j,l) (on this uniform periodic mesh the error has
#                  no other modes), and the rest; the amplitude and phase error of
#                  the computed wave mode; and the time history at T/4 ... T.
#                  Written, DIR/<variant>_N<N>_nel<nel>_field.dat (the error fields on
#                  a uniform grid), _spectrum.dat (the Fourier coefficients) and .png.
#   julia --project=. tools/alfven_linear_model.jl full 4 16,32 --field
#
# What it found (Nov 2026). The excess of the in-plane error over the
# out-of-plane one for even N (⊥/Bz 3.2 and 10.3 at N = 4 on 16² and 32², at
# t = 2) is not an error of the in-plane wave. It is a cancellation that the
# out-of-plane field gets at t = 2 and the in-plane field does not:
#   - CG-SEM on scalar advection converges at N for even N and at about N + 1
#     for odd N (1D, exact in time, 16 → 32 elements: 3.97-4.09 for N = 4 and
#     5.86-6.04 for N = 6 at t = 0.25, 0.5, 0.75, 1.5). The error is
#     element-scale (the rings); its O(h^N) part cancels whenever the wave has
#     crossed the periodic box a whole number of times: at t = 1 and 2 the order
#     is 5.5-5.8 for N = 4 and 7.7-7.9 for N = 6.
#   - t = 2 is one crossing of [0, √2]² along e∥. Bz obeys scalar advection and
#     sits on the cancellation; b·e⊥ is coupled through −c_h∇ψ, −∇p and
#     −∇(½|B|²) to ψ, p and B0·b, and does not cancel. With those three couplings
#     off (noglm_nogasp_nomagp) b·e⊥ equals Bz at every time; GLM alone carries
#     most of the difference. Off the crossing times there is no excess: at
#     t = 1.5, ⊥/Bz = 0.42, 0.60 (N = 4 on 8², 16²), 1.23, 1.10 (N = 5) and 1.23,
#     0.94 (N = 6), and both converge at about N for N = 4.
#   - The wave itself is computed equally well in both polarizations: its
#     Fourier mode holds under 0.01% of the error, and its amplitude and phase
#     errors are nearly equal for b·e⊥ and Bz (N = 4 on 32²: phase −8.7e-13 and
#     −8.8e-13 rad at t = 2, growing linearly in t).
#---------------------------------------------------------------------------------
using LinearAlgebra, SparseArrays, Printf
using OrdinaryDiffEq
function lgl(N)                      # LGL nodes and weights on [-1, 1]
    x = -cos.(pi*(0:N)/N)            # Chebyshev-Lobatto start
    P = zeros(N+1, N+1); xold = 2x
    while maximum(abs.(x - xold)) > 1e-15
        xold = copy(x)
        P[:,1] .= 1; P[:,2] .= x
        for k = 2:N
            P[:,k+1] = ((2k-1)*x.*P[:,k] - (k-1)*P[:,k-1])/k
        end
        x = xold - (x.*P[:,N+1] - P[:,N]) ./ ((N+1)*P[:,N+1])
    end
    w = 2 ./ (N*(N+1)*P[:,N+1].^2)
    return x, w
end
function gauss(n)                    # Gauss-Legendre nodes/weights (Golub-Welsch)
    β = [k/sqrt(4k^2 - 1) for k = 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), β))
    return E.values, 2 .* E.vectors[1,:].^2
end
function lagrange(ξ, x)              # L[i,j] = ℓ_j(x_i), D[i,j] = ℓ_j'(ξ_i)
    n = length(ξ)
    L = [prod((x[i] - ξ[m])/(ξ[j] - ξ[m]) for m in 1:n if m != j) for i in eachindex(x), j in 1:n]
    D = zeros(n, n)
    for i = 1:n, j = 1:n
        if i != j
            D[i,j] = prod((ξ[i]-ξ[m])/(ξ[j]-ξ[m]) for m in 1:n if m != j && m != i) / (ξ[j] - ξ[i])
        else
            D[i,j] = sum(1/(ξ[i]-ξ[m]) for m in 1:n if m != i)
        end
    end
    return L, D
end


# stubs for the deck's signatures
struct NSD_2D end; struct St_mesh end; abstract type CL end; struct CL0 <: CL end; struct TOTAL end
const TInt = Int64
const DECK = joinpath(@__DIR__, "..", "problems", "MHD", "alfvenWaveBabbarRanocha2026")
include(joinpath(DECK, "user_flux.jl"))
include(joinpath(DECK, "user_source.jl"))

const α  = π/4
const L  = sqrt(2.0)
const CH_DECK = 0.1 + sqrt(5.0/3.0*0.1 + 1.01)      # the deck's c_h (initialize.jl): 1.1847

function jacobians(ch)
    c_h_mhd[] = ch
    q0 = [1.0, 0.0, 0.0, 0.1/(γ_mhd - 1.0) + 0.5, 0.0, cos(α), sin(α), 0.0, 0.0]
    fx(q) = (F = similar(q); G = similar(q); user_flux!(F, G, NSD_2D(), q, nothing, St_mesh(), CL0(), TOTAL()); F)
    fy(q) = (F = similar(q); G = similar(q); user_flux!(F, G, NSD_2D(), q, nothing, St_mesh(), CL0(), TOTAL()); G)
    fs(q) = (S = similar(q); user_source!(S, q, nothing, 1, CL0(), TOTAL()); S)
    # complex step: J[:,k] = Im f(q0 + i h e_k)/h, exact to round-off for an analytic f
    function cjac(f)
        J = zeros(9, 9); hc = 1.0e-30
        for k = 1:9
            qc = complex.(q0); qc[k] += im*hc
            J[:, k] .= imag.(f(qc)) ./ hc
        end
        return J
    end
    # gradients of the gas pressure and of the magnetic pressure ½|B|², the two
    # pressures in the momentum fluxes (F[2], G[3] of user_flux.jl)
    gas(q) = [pressure_mhd(q[1], q[2], q[3], q[5], q[4], q[6], q[7], q[8], q[9])]
    mag(q) = [0.5*(q[6]*q[6] + q[7]*q[7] + q[8]*q[8])]
    function cgrad(f)
        g = zeros(9); hc = 1.0e-30
        for k = 1:9
            qc = complex.(q0); qc[k] += im*hc
            g[k] = imag(f(qc)[1])/hc
        end
        return g
    end
    return cjac(fx), cjac(fy), cjac(fs), q0, cgrad(gas), cgrad(mag)
end

# assembled scalar CG derivatives on the periodic [0, L]² mesh
function cg_ops(N, nel)
    ξ, w = lgl(N); _, D = lagrange(ξ, ξ)
    h = L/nel; J1 = h/2; n1 = nel*N; np = n1^2
    g1(e, i) = mod((e-1)*N + (i-1), n1)
    gid(ex, ey, i, j) = g1(ey, j)*n1 + g1(ex, i) + 1
    Ix = Int[]; Jx = Int[]; Vx = Float64[]; Iy = Int[]; Jy = Int[]; Vy = Float64[]
    M = zeros(np); xn = zeros(np); yn = zeros(np)
    for ey = 1:nel, ex = 1:nel, j = 1:N+1, i = 1:N+1
        g = gid(ex, ey, i, j); wij = w[i]*w[j]*J1^2
        M[g] += wij
        xn[g] = (ex-1)*h + (ξ[i]+1)*J1; yn[g] = (ey-1)*h + (ξ[j]+1)*J1
        for m = 1:N+1
            push!(Ix, g); push!(Jx, gid(ex, ey, m, j)); push!(Vx, wij*D[i,m]/J1)
            push!(Iy, g); push!(Jy, gid(ex, ey, i, m)); push!(Vy, wij*D[j,m]/J1)
        end
    end
    Mi = spdiagm(1 ./ M)
    return Mi*sparse(Ix, Jx, Vx, np, np), Mi*sparse(Iy, Jy, Vy, np, np), xn, yn, gid, ξ, J1, h
end

# the exact linear wave (amplitude A) in the deck's conservative variables
function wave(x, y, t; A = 1.0)
    s, c = sin(α), cos(α)
    f = A*sin(2π*(x*c + y*s + t)); g = A*cos(2π*(x*c + y*s + t))
    return (0.0, -s*f, c*f, 0.0, g, -s*f, c*f, g, 0.0)    # (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)'
end

const NAMES = ("ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ")

# the linear system for a variant, integrated to the times ts (the last one is
# the end of the integration, the others are read off Vern9's interpolant)
function solve_linear(N, nel; variant = "full", ts = [2.0])
    ch = occursin("noglm", variant) ? 0.0 : CH_DECK
    Ax, Ay, As, q0, gp, gm = jacobians(ch)
    # the pressure forces: p and ½|B|² enter the momentum only through F[2], G[3]
    occursin("nogasp", variant) && (Ax[2, :] .-= gp; Ay[3, :] .-= gp)
    occursin("nomagp", variant) && (Ax[2, :] .-= gm; Ay[3, :] .-= gm)
    Dx, Dy, xn, yn, gid, ξ, J1, h = cg_ops(N, nel)
    np = length(xn)
    Lop = -(kron(sparse(Ax), Dx) + kron(sparse(Ay), Dy)) + kron(sparse(As), sparse(I, np, np))
    if occursin("noinddiv", variant)        # + B0 (∇·u) on Bx, By (u = ρu at ρ0 = 1)
        Qx = spzeros(9, 9); Qy = spzeros(9, 9)
        Qx[6,2] = q0[6]; Qx[7,2] = q0[7]; Qy[6,3] = q0[6]; Qy[7,3] = q0[7]
        Lop += kron(Qx, Dx) + kron(Qy, Dy)
    end
    if occursin("powell", variant)          # −B0 (∇·b) on ρu, ρv
        Px = spzeros(9, 9); Py = spzeros(9, 9)
        Px[2,6] = q0[6]; Px[3,6] = q0[7]; Py[2,7] = q0[6]; Py[3,7] = q0[7]
        Lop -= kron(Px, Dx) + kron(Py, Dy)
    end
    Q0 = zeros(np, 9)
    for g = 1:np
        Q0[g, :] .= wave(xn[g], yn[g], 0.0)
    end
    # saveat does not change the steps: the final state is the same as without it
    prob = ODEProblem((dq, q, p, s) -> mul!(dq, Lop, q), vec(Q0), (0.0, ts[end]))
    sol  = solve(prob, Vern9(); reltol = 1e-12, abstol = 1e-14, saveat = ts)
    return (Q = [reshape(u, np, 9) for u in sol.u], ts = sol.t, N = N, nel = nel,
            gid = gid, ξ = ξ, J1 = J1, h = h, gp = gp)
end

# the error fields at points (xa, yb) of [-1, 1]² in element (ex, ey) of state Q
# at time t, from the element polynomial: E[a, b, f] for f in FIELDS
const FIELDS = ("b·e⊥", "b·e∥", "Bz", "ψ", "p'", "u·e⊥")
function element_errors!(E, uq, ue, s, Q, t, ex, ey, Lp)
    N = s.N; sα, cα = sin(α), cos(α)
    for k = 1:9
        for j = 1:N+1, i = 1:N+1
            ue[i,j] = Q[s.gid(ex, ey, i, j), k]
        end
        uq[:, :, k] .= Lp*ue*Lp'
    end
    xp = (ex-1)*s.h .+ (Lp*s.ξ .+ 1) .* s.J1      # Lp reproduces the nodes' coordinates
    yp = (ey-1)*s.h .+ (Lp*s.ξ .+ 1) .* s.J1
    for b in axes(uq, 2), a in axes(uq, 1)
        qe = wave(xp[a], yp[b], t)
        e  = ntuple(k -> uq[a,b,k] - qe[k], 9)
        E[a,b,1] = -sα*e[6] + cα*e[7]                   # b·e⊥, the wave
        E[a,b,2] =  cα*e[6] + sα*e[7]                   # b·e∥, exactly 0
        E[a,b,3] = e[8]                                  # Bz, the scalar reference
        E[a,b,4] = e[9]                                  # ψ, exactly 0
        E[a,b,5] = sum(s.gp[k]*e[k] for k = 1:9)        # p' (linearized), exactly 0
        E[a,b,6] = -sα*e[2] + cα*e[3]                   # u·e⊥ (ρ0 = 1), the wave
    end
    return xp, yp
end

# RMS errors per component at N+4 Gauss points (as the case measures them)
function rms_errors(s, Q, t)
    N = s.N; nel = s.nel
    xg, wg = gauss(N + 4); Lg, _ = lagrange(s.ξ, xg)
    # entries 1-9: the components; 10, 11: b·e⊥ (the wave) and b·e∥ (exactly 0)
    s2 = zeros(11); uq = zeros(length(xg), length(xg), 9); ue = zeros(N+1, N+1)
    sα, cα = sin(α), cos(α)
    for ey = 1:nel, ex = 1:nel
        for k = 1:9
            for j = 1:N+1, i = 1:N+1
                ue[i,j] = Q[s.gid(ex, ey, i, j), k]
            end
            uq[:, :, k] .= Lg*ue*Lg'
        end
        for b in eachindex(xg), a in eachindex(xg)
            x = (ex-1)*s.h + (xg[a]+1)*s.J1; y = (ey-1)*s.h + (xg[b]+1)*s.J1
            qe = wave(x, y, t); wab = wg[a]*wg[b]*s.J1^2
            e = [uq[a,b,k] - qe[k] for k = 1:9]
            s2[1:9] .+= wab .* e.^2
            s2[10] += wab*(-sα*e[6] + cα*e[7])^2
            s2[11] += wab*( cα*e[6] + sα*e[7])^2
        end
    end
    return sqrt.(s2 ./ L^2)
end

#---------------------------------------------------------------------------------
# --field: where the error is, in space and in wavenumber.
# The mesh is uniform and periodic and the system linear with constant
# coefficients, so the discrete solution is a Bloch wave: started from the
# wave's mode k0 = ±(1,1) (in units of 2π/L), it only ever contains the modes
# k0 + nel·(j, l). The error therefore splits exactly into
#   the wave mode  ±(1,1)            an error of the wave itself: amplitude and
#                                     phase (speed) of the computed wave, smooth
#                                     in space; it grows with t if it is a
#                                     wrong speed or a damping
#   ring r         ±(1,1)+nel·(j,l),  structure on the element scale (r = 1:
#                  max(|j|,|l|) = r   wavelength ~ one element), the part of the
#                                     error that lives inside the elements
# with Σ |ê|² = RMS² (Parseval), ê(m,n) = (1/L²) ∫ e exp(−2πi(mx + ny)/L).
# The element polynomials oscillate ~(N+1)/2 times across an element, so the
# rings are computed up to R = N + 2; "rest" is RMS² minus all of them.
#---------------------------------------------------------------------------------
function wave_modes(nel, R)
    ks = Tuple{Int,Int}[]
    for sgn in (1, -1), j = -R:R, l = -R:R
        k = (sgn + nel*j, sgn + nel*l)
        k in ks || push!(ks, k)
    end
    return ks
end
# ring of mode k: 0 for the wave mode, r for ±(1,1) + nel·(j,l) with max(|j|,|l|) = r
function ring(k, nel)
    r = typemax(Int)
    for sgn in (1, -1)
        dm, dn = k[1] - sgn, k[2] - sgn
        (mod(dm, nel) == 0 && mod(dn, nel) == 0) && (r = min(r, max(abs(dm), abs(dn)) ÷ nel))
    end
    return r
end

# Fourier coefficients of the error fields at the modes ks, by Gauss quadrature
# of the element polynomials, separable: for each element Φx (W E) Φy over all
# the m and n of ks. A mode m turns by π m/nel across an element in reference
# coordinates; N + that + 12 Gauss points integrate it to round-off.
function error_spectrum(s, Q, t, ks)
    N = s.N; nel = s.nel; nf = length(FIELDS)
    mv = sort(unique(first.(ks))); nv = sort(unique(last.(ks)))
    xg, wg = gauss(N + ceil(Int, π*maximum(abs, [mv; nv])/nel) + 12)
    Lg, _ = lagrange(s.ξ, xg); nq = length(xg)
    E = zeros(nq, nq, nf); uq = zeros(nq, nq, 9); ue = zeros(N+1, N+1)
    G = zeros(ComplexF64, length(mv), length(nv), nf); tot = zeros(nf)
    for ey = 1:nel, ex = 1:nel
        xp, yp = element_errors!(E, uq, ue, s, Q, t, ex, ey, Lg)
        Φx = [cis(-2π*m*xp[a]/L)*wg[a]*s.J1 for m in mv, a = 1:nq]
        Φy = [cis(-2π*n*yp[b]/L)*wg[b]*s.J1 for b = 1:nq, n in nv]
        for f = 1:nf
            Ef = view(E, :, :, f)
            G[:, :, f] .+= Φx * Ef * Φy
            tot[f] += sum(wg[a]*wg[b]*Ef[a,b]^2 for a = 1:nq, b = 1:nq)*s.J1^2
        end
    end
    C = [G[searchsortedfirst(mv, k[1]), searchsortedfirst(nv, k[2]), f] for k in ks, f = 1:nf]
    return C ./ L^2, tot ./ L^2                     # Σ_all |ê|² = tot (Parseval)
end

# the computed wave's own mode against the exact one: amplitude ratio − 1 and
# phase error (rad; > 0: the computed wave is ahead) of b·e⊥ and Bz
function wave_mode_error(C, ks, t)
    i11 = findfirst(==((1, 1)), ks)
    xh  = cis(2π*t)                                  # exp(2πi(x∥ + t)) is mode (1,1)
    ref = (xh/(2im), xh/2)                           # exact b⊥ = sin, Bz = cos, at (1,1)
    r   = ((C[i11, 1] + ref[1])/ref[1], (C[i11, 3] + ref[2])/ref[2])
    return (abs(r[1]) - 1, angle(r[1])), (abs(r[2]) - 1, angle(r[2]))
end

# RMS² of field f in the wave mode, in each ring 1..R, and the rest (> R)
ring_parts(C, tot, rings, f, R) =
    (ps = [sum(abs2(C[i, f]) for i in eachindex(rings) if rings[i] == r; init = 0.0) for r = 0:R];
     [ps; max(tot[f] - sum(ps), 0.0)])

function spectrum_report(io, s, variant)
    R = s.N + 2
    ks = wave_modes(s.nel, R); rings = [ring(k, s.nel) for k in ks]
    C, tot = error_spectrum(s, s.Q[end], s.ts[end], ks)
    @printf(io, "#\n#  %s  N=%d nel=%d  t=%g   error by Fourier mode: RMS of each part, %% of RMS²\n",
            variant, s.N, s.nel, s.ts[end])
    @printf(io, "#  %-6s %11s  %11s %5s  %11s %5s  %11s %5s  %11s %5s  %9s %5s\n", "field", "RMS",
            "wave mode", "%", "ring 1", "%", "ring 2", "%", "rings 3-$R", "%", "rest", "%")
    for f in eachindex(FIELDS)
        ps = ring_parts(C, tot, rings, f, R)
        parts = [ps[1], ps[2], ps[3], sum(ps[4:R+1]), ps[end]]
        @printf(io, "#  %-6s %11.4e", FIELDS[f], sqrt(tot[f]))
        for (i, x) in enumerate(parts)
            sh = tot[f] > 0 ? 100x/tot[f] : 0.0
            i < 5 ? @printf(io, "  %11.4e %5.1f", sqrt(x), sh) : @printf(io, "  %9.2e %5.1f", sqrt(x), sh)
        end
        println(io)
    end
    (aw, pw), (az, pz) = wave_mode_error(C, ks, s.ts[end])
    @printf(io, "#  wave mode, computed/exact:  b·e⊥ amplitude %+.3e phase %+.3e rad;  Bz amplitude %+.3e phase %+.3e rad\n",
            aw, pw, az, pz)
    # in time: an error of the wave's speed or damping grows ~ t in the wave
    # mode; element-scale error (the rings) does not accumulate
    if length(s.ts) > 1
        kw = [(1, 1), (-1, -1)]
        @printf(io, "#  %5s   %-24s   %-24s   %s\n", "t", "b·e⊥: wave mode  rings", "Bz: wave mode  rings",
                "b·e⊥ wave-mode phase (rad)")
        for (it, t) in enumerate(s.ts)
            Ct, tt = error_spectrum(s, s.Q[it], t, kw)
            w⊥ = sum(abs2, Ct[:, 1]); wz = sum(abs2, Ct[:, 3])
            (_, φ), _ = wave_mode_error(Ct, kw, t)
            @printf(io, "#  %5.2f   %11.4e %11.4e   %11.4e %11.4e   %+.3e\n", t, sqrt(w⊥),
                    sqrt(max(tt[1] - w⊥, 0.0)), sqrt(wz), sqrt(max(tt[3] - wz, 0.0)), φ)
        end
    end
    return ks, rings, C, tot, R
end

# the error fields on a uniform grid (cell centres, nsub per element and
# direction), for plotting: x y and one column per FIELDS entry
function write_fields(file, s, Q, t)
    nsub = max(2, min(2*s.N, 256 ÷ s.nel))
    ξs = [-1 + (2a - 1)/nsub for a = 1:nsub]; Ls, _ = lagrange(s.ξ, ξs)
    E = zeros(nsub, nsub, length(FIELDS)); uq = zeros(nsub, nsub, 9); ue = zeros(s.N+1, s.N+1)
    n1 = s.nel*nsub; X = zeros(n1); Y = zeros(n1); F = zeros(n1, n1, length(FIELDS))
    for ey = 1:s.nel, ex = 1:s.nel
        xp, yp = element_errors!(E, uq, ue, s, Q, t, ex, ey, Ls)
        ia = (ex-1)*nsub .+ (1:nsub); ib = (ey-1)*nsub .+ (1:nsub)
        X[ia] .= xp; Y[ib] .= yp; F[ia, ib, :] .= E
    end
    open(file, "w") do io
        println(io, "# x y ", join(FIELDS, " "), "   (errors at t = $t; b·e∥, ψ, p' are 0 in the exact solution)")
        for b = 1:n1, a = 1:n1
            @printf(io, "%.6f %.6f", X[a], Y[b])
            for f in eachindex(FIELDS)
                @printf(io, " % .6e", F[a, b, f])
            end
            println(io)
        end
    end
    return X, Y, F
end

# the spectrum, one line per computed mode: m n ring, then Re Im of ê per field
function write_spectrum(file, ks, rings, C)
    open(file, "w") do io
        println(io, "# m n ring  then Re(ê) Im(ê) for ", join(FIELDS, " "),
                "   (modes ±(1,1) + nel·(j,l) in units of 2π/L; ring 0 = the wave mode)")
        for i in sortperm(rings)
            @printf(io, "%4d %4d %2d", ks[i][1], ks[i][2], rings[i])
            for f in eachindex(FIELDS)
                @printf(io, " % .6e % .6e", real(C[i, f]), imag(C[i, f]))
            end
            println(io)
        end
    end
end

# one PNG per case: the error fields (diverging: blue < 0 < red, gray at 0,
# symmetric limits per panel) and where their RMS sits in wavenumber
function plot_case(file, X, Y, F, rings, C, tot, R, variant, N, nel, t)
    div = Plots.cgrad(["#2a78d6", "#f0efec", "#e34948"])
    # each field over its own max |e|, given in the title
    panel(f) = (m = maximum(abs, view(F, :, :, f)); m = m > 0 ? m : 1.0;
                Plots.heatmap(X, Y, permutedims(F[:, :, f]) ./ m; c = div, clims = (-1, 1),
                              aspect_ratio = :equal, xlims = (0, L), ylims = (0, L),
                              title = @sprintf("%s / %.2e    (RMS %.2e)", FIELDS[f], m, sqrt(tot[f])),
                              titlefontsize = 10, tickfontsize = 7, framestyle = :box))
    # where the error is: RMS of the wave mode, of each ring and of the rest
    sp = Plots.plot(; yscale = :log10, xticks = (0:R+1, ["wave"; string.(1:R); "rest"]),
                    xlabel = "ring", title = "error by Fourier mode (RMS)", titlefontsize = 10,
                    tickfontsize = 7, guidefontsize = 8, legend = :bottomright, legendfontsize = 7,
                    framestyle = :box, grid = :y, gridalpha = 0.15)
    for (f, c) in ((1, "#2a78d6"), (3, "#eb6834"), (2, "#1baf7a"), (4, "#eda100"))
        tot[f] > 0 || continue
        v = sqrt.(max.(ring_parts(C, tot, rings, f, R), 1e-300))
        Plots.plot!(sp, 0:R+1, v; label = FIELDS[f], color = c, lw = 2, marker = :circle,
                    ms = 5, markerstrokecolor = "#fcfcfb", markerstrokewidth = 2)
    end
    p = Plots.plot(panel(1), panel(2), panel(4), panel(3), panel(5), sp; layout = (2, 3),
                   size = (1350, 860), plot_title = @sprintf("%s   N = %d, %d² elements, t = %g:  errors",
                                                             variant, N, nel, t),
                   plot_titlefontsize = 12, left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)
    Plots.savefig(p, file)
end

#---------------------------------------------------------------------------------
flags = filter(a -> startswith(a, "--"), ARGS)
pargs = filter(a -> !startswith(a, "--"), ARGS)
fielddir = nothing; tend = 2.0
for f in flags
    if f == "--field"
        global fielddir = "alfven_linear_fields"
    elseif startswith(f, "--field=")
        global fielddir = String(split(f, "="; limit = 2)[2])
    elseif startswith(f, "--t=")
        global tend = parse(Float64, split(f, "="; limit = 2)[2])
    else
        error("unknown option \"$f\": the options are --field[=DIR] and --t=T")
    end
end
variant = length(pargs) > 0 ? pargs[1] : "full"
# every word of the variant must be a switch: an unknown one is an error, not a
# silent "full" (an older copy of this tool, without nogasp/nomagp/noinddiv,
# ran those as full and gave full's numbers)
const SWITCHES = ("noglm", "powell", "nogasp", "nomagp", "noinddiv")
active = variant == "full" ? String[] : String.(split(variant, "_"))
for f in active
    f in SWITCHES || error("unknown switch \"$f\" in \"$variant\": use full, or any of " *
                           join(SWITCHES, ", ") * " joined by _")
end
println("# switches: ", isempty(active) ? "none (the deck, linearized)" : join(active, ", "))
Ns   = length(pargs) > 1 ? parse.(Int, split(pargs[2], ",")) : [4, 5, 6]
nels = length(pargs) > 2 ? parse.(Int, split(pargs[3], ",")) : [8, 16]
if fielddir !== nothing
    mkpath(fielddir)
    ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")     # GR without a window
    @eval using Plots
end
# with --field, also the states at T/4, T/2, 3T/4 for the time history
ts = fielddir === nothing ? [tend] : tend .* [0.25, 0.5, 0.75, 1.0]
println("# final time t = $tend", tend == 2.0 ? " (the deck's; one crossing of the box along e∥)" : "")
@printf("%-36s %2s %4s %11s %11s %11s %7s %7s %10s %10s\n", "variant", "N", "nel",
        "b·e⊥ (wave)", "b·e∥ (=0)", "Bz (scalar)", "⊥/Bz", "∥/Bz", "ψ", "ρE")
cases = []
for N in Ns, nel in nels
    sl = solve_linear(N, nel; variant = variant, ts = ts)
    e  = rms_errors(sl, sl.Q[end], sl.ts[end])
    @printf("%-36s %2d %4d %11.4e %11.4e %11.4e %7.2f %7.2f %10.3e %10.3e\n",
            variant, N, nel, e[10], e[11], e[8], e[10]/e[8], e[11]/e[8], e[9], e[4])
    flush(stdout)
    fielddir === nothing || push!(cases, sl)
end
for sl in cases
    ks, rings, C, tot, R = spectrum_report(stdout, sl, variant)
    base = joinpath(fielddir, @sprintf("%s_N%d_nel%d", variant, sl.N, sl.nel))
    X, Y, F = write_fields(base * "_field.dat", sl, sl.Q[end], sl.ts[end])
    write_spectrum(base * "_spectrum.dat", ks, rings, C)
    plot_case(base * ".png", X, Y, F, rings, C, tot, R, variant, sl.N, sl.nel, sl.ts[end])
    println("#  wrote ", base, "_field.dat, _spectrum.dat, .png")
    flush(stdout)
end
