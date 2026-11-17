#=
sw_ref_1d.jl : plane-parallel shortwave reference solver (single column), standard library only.

  μ dI/dτ = −I + (ω/2) ∫ p̄(μ,μ′) I dμ′ + (ω S₀/4π) p̄(μ,μ₀) e^(−τ/μ₀)          (diffuse radiance, azimuth-averaged)
  F_dir   = μ₀ S₀ e^(−τ/μ₀)                                                    (direct beam, Beer–Lambert)

Conventions: τ grows downward from TOA; μ > 0 is DOWNWARD travel; S₀ is the beam flux normal to the beam;
p normalised so (1/4π)∫p dΩ = 1 (the 3D code's Φ = p/4π). Henyey–Greenstein: χ_l = g^l.
Delta-M (Wiscombe 1977): f = χ_{2N}; σ′ = (1−f)σ, κ_abs′ = κ_abs, χ_l′ = (χ_l − f)/(1−f).
Surface: Lambertian, I↑ = (A/π)(F_dir + F↓_diffuse).

Two independent solvers:
  solve_column  : inhomogeneous column, diamond difference on a fine τ grid, one sparse direct solve.
  doubling_slab : homogeneous slab by doubling–adding (+ Lambertian surface by adding).
=#
module SWRef1D
using LinearAlgebra, SparseArrays, Printf

"Gauss–Legendre nodes/weights on [-1,1] (Golub–Welsch)."
function gauss_legendre(n::Int)
    β = [k / sqrt(4k^2 - 1) for k = 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), β))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end

"Double-Gauss: N points on (0,1) for μ>0; returns μ (2N: first N downward >0, next N upward <0) and weights (sum 2)."
function double_gauss(N::Int)
    x, w = gauss_legendre(N)
    μp = (x .+ 1) ./ 2; wp = w ./ 2
    return vcat(μp, -μp), vcat(wp, wp)
end

"Legendre polynomials P_0..P_L at x."
function legendre_all(L::Int, x::Float64)
    P = zeros(L + 1); P[1] = 1.0
    L >= 1 && (P[2] = x)
    for l = 1:L-1
        P[l+2] = ((2l + 1) * x * P[l+1] - l * P[l]) / (l + 1)
    end
    return P
end

"Azimuth-averaged phase matrix p̄(μ_i, μ_j) = Σ (2l+1) χ_l P_l(μ_i) P_l(μ_j), l = 0..L."
function phase_matrix(μi::Vector{Float64}, μj::Vector{Float64}, χ::Vector{Float64})
    L = length(χ) - 1
    Pi = hcat([legendre_all(L, m) for m in μi]...)   # (L+1) × ni
    Pj = hcat([legendre_all(L, m) for m in μj]...)
    c = [(2l + 1) * χ[l+1] for l = 0:L]
    return Pi' * Diagonal(c) * Pj
end

"Optical properties after delta-M (or none). Returns κ_abs, σ′, χ′ (length 2N) for one level."
function scaled_props(κa, σ, g, N; deltam::Bool=true, nleg::Int=2N)
    χ = [g^l for l = 0:nleg-1]
    f = deltam ? g^(nleg) : 0.0
    χs = deltam ? (χ .- f) ./ (1 - f) : χ
    return κa, (1 - f) * σ, χs
end

"""
    solve_column(z, κa, σ, g; S0, μ0, A, N, deltam, dτmax) -> NamedTuple

Column given at nodes z (any order) with κ_abs, σ [m⁻¹] and g at the nodes, linear between nodes.
Returns, at the ORIGINAL nodes (sorted top→bottom): z, G (W/m²·4π-normalised mean radiance ×4π), Fdn, Fup, Fdir, Q = κ_abs·G,
plus column budget: F_toa_in, F_toa_up, F_sfc_net, absorbed (∫Q dz) and absorbed_flux (F_net,top − F_net,sfc).
"""
function solve_column(z, κa, σ, g; S0=1.0, μ0=0.5, A=0.0, N=16, deltam::Bool=true, dτmax=0.05)
    o = sortperm(z; rev=true); z = Float64.(z[o]); κa = Float64.(κa[o]); σ = Float64.(σ[o]); g = Float64.(g[o])
    # fine grid: subdivide each interval so the scaled Δτ′ ≤ dτmax (properties linear in z)
    zf = Float64[z[1]]; ka_f = Float64[κa[1]]; sg_f = Float64[σ[1]]; g_f = Float64[g[1]]; orig = [1]
    for k = 1:length(z)-1
        f1 = deltam ? g[k]^(2N) : 0.0; f2 = deltam ? g[k+1]^(2N) : 0.0
        e1 = κa[k] + (1 - f1) * σ[k]; e2 = κa[k+1] + (1 - f2) * σ[k+1]
        dτ = 0.5 * (e1 + e2) * (z[k] - z[k+1])
        ns = max(1, ceil(Int, dτ / dτmax))
        for s = 1:ns
            t = s / ns
            push!(zf, z[k] + t * (z[k+1] - z[k])); push!(ka_f, κa[k] + t * (κa[k+1] - κa[k]))
            push!(sg_f, σ[k] + t * (σ[k+1] - σ[k])); push!(g_f, g[k] + t * (g[k+1] - g[k]))
        end
        push!(orig, length(zf))
    end
    K = length(zf)                                   # levels
    μ, w = double_gauss(N); nd = 2N
    # per-level scaled properties and ω′, p̄
    ωl = zeros(K); ext = zeros(K); Pl = Vector{Matrix{Float64}}(undef, K); pbeam = Vector{Vector{Float64}}(undef, K)
    cache = Dict{Float64,Tuple{Matrix{Float64},Vector{Float64},Float64}}()   # phase matrices keyed by g (rounded 1e-4)
    for k = 1:K
        gk = round(g_f[k], digits=4)
        Pm, pb, f = get!(cache, gk) do
            _, _, χs = scaled_props(0.0, 1.0, gk, N; deltam=deltam)
            (phase_matrix(μ, μ, χs), vec(phase_matrix(μ, [μ0], χs)), deltam ? gk^(2N) : 0.0)
        end
        ext[k] = ka_f[k] + (1 - f) * sg_f[k]; ωl[k] = ext[k] > 0 ? (1 - f) * sg_f[k] / ext[k] : 0.0
        Pl[k] = Pm; pbeam[k] = pb
    end
    τ = zeros(K); for k = 2:K; τ[k] = τ[k-1] + 0.5 * (ext[k-1] + ext[k]) * (zf[k-1] - zf[k]); end
    # direct beam with the SAME z-centred scheme as the diffuse field:
    #   μ0 (b_{k+1} - b_k)/dz = -½(e_k b_k + e_{k+1} b_{k+1})
    # so the beam's loss in each cell equals exactly the trapezoid of what it deposits (scattered + absorbed).
    beam = zeros(K); beam[1] = S0
    for k = 1:K-1
        dz = zf[k] - zf[k+1]
        beam[k+1] = beam[k] * (μ0 / dz - 0.5 * ext[k]) / (μ0 / dz + 0.5 * ext[k+1])
    end
    # unknown index: level k, direction i  -> (k-1)*nd + i
    idx(k, i) = (k - 1) * nd + i
    I_ = Int[]; J_ = Int[]; V_ = Float64[]; b = zeros(K * nd); row = 0
    addc!(r, c, v) = (push!(I_, r); push!(J_, c); push!(V_, v))
    # source J_k,i = (ω/2) Σ_j w_j p̄_ij I_k,j + (ω S0/4π) p̄_i0 e^{-τ/μ0}
    sg = ωl .* ext                                    # scaled scattering coefficient σ′ at the levels
    for k = 1:K-1, i = 1:nd                           # diamond difference IN z on cell [k,k+1] (depth increment dz > 0)
        row += 1; dz = zf[k] - zf[k+1]
        # μ_i (I_{k+1}-I_k)/dz + ½(e_k I_k + e_{k+1} I_{k+1}) - ½(S_k + S_{k+1}) = 0,
        # S_k = (σ′_k/2) Σ_j w_j p̄_k,ij I_k,j + σ′_k S0/(4π) p̄_k,i0 e^{-τ_k/μ0}
        # Summed over directions, the cell's diffuse absorption is exactly the trapezoid ½dz(κa_k G_k + κa_{k+1} G_{k+1}).
        addc!(row, idx(k+1, i), μ[i] / dz + 0.5 * ext[k+1]); addc!(row, idx(k, i), -μ[i] / dz + 0.5 * ext[k])
        for kk in (k, k + 1), j = 1:nd
            addc!(row, idx(kk, j), -0.5 * 0.5 * sg[kk] * w[j] * Pl[kk][i, j])
        end
        b[row] = 0.5 * sum(sg[kk] * pbeam[kk][i] / (4π) * beam[kk] for kk in (k, k + 1))
    end
    for i = 1:N                                       # TOA: no diffuse inflow downward
        row += 1; addc!(row, idx(1, i), 1.0); b[row] = 0.0
    end
    Fdir_s = μ0 * beam[K]
    for i = N+1:nd                                    # surface: Lambertian  I↑ = (A/π)(F_dir + 2π Σ μ_j w_j I↓_j)
        row += 1; addc!(row, idx(K, i), 1.0)
        for j = 1:N; addc!(row, idx(K, j), -(A / π) * 2π * μ[j] * w[j]); end
        b[row] = (A / π) * Fdir_s
    end
    Amat = sparse(I_, J_, V_, row, K * nd)
    x = Amat \ b
    Ik = reshape(x, nd, K)
    G = [2π * sum(w .* Ik[:, k]) + beam[k] for k = 1:K]
    Fdn = [μ0 * beam[k] + 2π * sum(μ[1:N] .* w[1:N] .* Ik[1:N, k]) for k = 1:K]
    Fup = [2π * sum(abs.(μ[N+1:nd]) .* w[N+1:nd] .* Ik[N+1:nd, k]) for k = 1:K]
    Q = ka_f .* G
    absorbed = sum(0.5 * (Q[k] + Q[k+1]) * (zf[k] - zf[k+1]) for k = 1:K-1)
    Fnet_top = Fdn[1] - Fup[1]; Fnet_sfc = Fdn[K] - Fup[K]
    return (z = zf[orig], G = G[orig], Fdn = Fdn[orig], Fup = Fup[orig], Fdir = μ0 .* beam[orig], Q = Q[orig],
            F_toa_in = μ0 * S0, F_toa_up = Fup[1], F_sfc_net = Fnet_sfc, absorbed = absorbed,
            absorbed_flux = Fnet_top - Fnet_sfc, nlevels = K, μ = μ, Idiff = Ik[:, orig])
end

"""
    doubling_slab(τ, ω, g; S0, μ0, A, N, deltam) -> NamedTuple

Homogeneous slab by doubling–adding (azimuth-averaged), same quadrature as solve_column.
Returns reflected flux at TOA, total downward flux at the surface, net surface flux, and column absorption.
"""
function doubling_slab(τtot, ω, g; S0=1.0, μ0=0.5, A=0.0, N=16, deltam::Bool=true)
    χ = [g^l for l = 0:2N-1]; f = deltam ? g^(2N) : 0.0
    χs = deltam ? (χ .- f) ./ (1 - f) : χ
    τs = (1 - ω * f) * τtot; ωs = (1 - f) * ω / (1 - ω * f)          # delta-M scaled slab
    x, wq = gauss_legendre(N); μp = (x .+ 1) ./ 2; wp = wq ./ 2
    Ppp = phase_matrix(μp, μp, χs); Ppm = phase_matrix(-μp, μp, χs)   # down→down, down→up
    b_dn = vec(phase_matrix(μp, [μ0], χs)); b_up = vec(phase_matrix(-μp, [μ0], χs))
    M⁻¹ = Diagonal(1 ./ μp); W = Diagonal(wp)
    nd = ceil(Int, log2(τs / 2.0^-20)); δτ = τs / 2.0^nd
    R = δτ * M⁻¹ * (ωs / 2) * Ppm * W
    T = I - δτ * M⁻¹ + δτ * M⁻¹ * (ωs / 2) * Ppp * W
    s_up = δτ * M⁻¹ * (ωs * S0 / 4π) * b_up; s_dn = δτ * M⁻¹ * (ωs * S0 / 4π) * b_dn
    e = exp(-δτ / μ0)
    for _ = 1:nd                                         # combine two identical layers
        X = inv(I - R * R)
        D = X * (s_dn + R * (e * s_up))                  # downward diffuse at the interface
        U = e * s_up + R * D                             # upward diffuse at the interface
        s_up, s_dn = s_up + T * U, e * s_dn + T * D
        R, T = R + T * X * R * T, T * X * T
        e = e^2
    end
    Fdir_s = μ0 * S0 * exp(-τs / μ0)
    # Lambertian surface: U_s = (A/π)(F_dir + 2π Σ μ w D_s);  D_s = s_dn + R U_s
    ρ = (A / π) * 2π * (μp .* wp)'                       # row vector: maps D_s to the isotropic upward radiance
    one_ = ones(N)
    # D_s = s_dn + R*one_*(A/π)(Fdir_s + 2π μw·D_s)
    Mx = I - R * one_ * ρ
    D_s = Mx \ (s_dn + R * one_ * (A / π) * Fdir_s)
    U_s = one_ * ((A / π) * Fdir_s + (ρ*D_s)[1])
    up_top = s_up + T * U_s
    F_toa_up = 2π * sum(μp .* wp .* up_top)
    F_sfc_dn = Fdir_s + 2π * sum(μp .* wp .* D_s)
    F_sfc_up = 2π * sum(μp .* wp .* U_s)
    absorbed = μ0 * S0 - F_toa_up - (F_sfc_dn - F_sfc_up)
    return (F_toa_up = F_toa_up, F_sfc_dn = F_sfc_dn, F_sfc_net = F_sfc_dn - F_sfc_up, absorbed = absorbed)
end

end # module
