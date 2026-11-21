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
# What it found (Nov 2026): the even-N excess of the in-plane error is linear
# and the Powell term does not touch it. It reaches the WAVE (b·e⊥) through
# three gradient forces, any one of which suffices: −∇p, −∇(½|B|²) and the GLM
# −c_h∇ψ. With all three off (noglm_nogasp_nomagp) b·e⊥ equals the scalar
# reference to every printed digit at N = 4, 5, 6 on 8² and 16², and the field
# compression (noinddiv) changes b∥ but not the wave; GLM alone
# (nogasp_nomagp) carries the whole excess (⊥/Bz on 16²: 3.08, 0.80, 2.57 for
# N = 4, 5, 6, against 3.23, 0.56, 2.11 in full). p', b∥ and ψ are zero in the
# exact solution: what those forces push back into the wave is driven by the
# O(h^N) discrete divergence of the in-plane fields. Why that comes back at
# O(h^N) for even N only is not settled.
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

function run(N, nel; variant = "full", t = 2.0)
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
    prob = ODEProblem((dq, q, p, s) -> mul!(dq, Lop, q), vec(Q0), (0.0, t))
    sol  = solve(prob, Vern9(); reltol = 1e-12, abstol = 1e-14, save_everystep = false)
    Q = reshape(sol.u[end], np, 9)
    # RMS errors per component at N+4 Gauss points
    xg, wg = gauss(N + 4); Lg, _ = lagrange(ξ, xg)
    # entries 1-9: the components; 10, 11: b·e⊥ (the wave) and b·e∥ (exactly 0)
    s2 = zeros(11); uq = zeros(length(xg), length(xg), 9); ue = zeros(N+1, N+1)
    sα, cα = sin(α), cos(α)
    for ey = 1:nel, ex = 1:nel
        for k = 1:9
            for j = 1:N+1, i = 1:N+1
                ue[i,j] = Q[gid(ex, ey, i, j), k]
            end
            uq[:, :, k] .= Lg*ue*Lg'
        end
        for b in eachindex(xg), a in eachindex(xg)
            x = (ex-1)*h + (xg[a]+1)*J1; y = (ey-1)*h + (xg[b]+1)*J1
            qe = wave(x, y, t); wab = wg[a]*wg[b]*J1^2
            e = [uq[a,b,k] - qe[k] for k = 1:9]
            s2[1:9] .+= wab .* e.^2
            s2[10] += wab*(-sα*e[6] + cα*e[7])^2
            s2[11] += wab*( cα*e[6] + sα*e[7])^2
        end
    end
    return sqrt.(s2 ./ L^2)
end

variant = length(ARGS) > 0 ? ARGS[1] : "full"
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
Ns   = length(ARGS) > 1 ? parse.(Int, split(ARGS[2], ",")) : [4, 5, 6]
nels = length(ARGS) > 2 ? parse.(Int, split(ARGS[3], ",")) : [8, 16]
@printf("%-36s %2s %4s %11s %11s %11s %7s %7s %10s %10s\n", "variant", "N", "nel",
        "b·e⊥ (wave)", "b·e∥ (=0)", "Bz (scalar)", "⊥/Bz", "∥/Bz", "ψ", "ρE")
for N in Ns, nel in nels
    e = run(N, nel; variant = variant)
    @printf("%-36s %2d %4d %11.4e %11.4e %11.4e %7.2f %7.2f %10.3e %10.3e\n",
            variant, N, nel, e[10], e[11], e[8], e[10]/e[8], e[11]/e[8], e[9], e[4])
    flush(stdout)
end
