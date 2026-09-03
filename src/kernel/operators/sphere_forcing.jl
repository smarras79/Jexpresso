#---------------------------------------------------------------------------------
# sphere_forcing.jl
#
# FORCING AND LARGE-SCALE DISSIPATION for the shallow water equations on the
# spherical shell — the forced-dissipative setup of
#
#   Scott, R. K. & Polvani, L. M. (2007), "Forced-dissipative shallow-water
#   turbulence on the sphere and the atmospheric circulation of the giant
#   planets", J. Atmos. Sci. 64, 3158-3176.
#
# The paper writes the equations in (vorticity ζ, divergence δ, height h),
# Eq. (9), and forces the VORTICITY equation with a random, isotropic,
# small-scale process F, confined to a band of spherical-harmonic total
# wavenumbers around n_f, Eq. (11). Large-scale dissipation is a linear term
# D_ξ on ζ, δ (Rayleigh friction) or on h (radiative relaxation), Eq. (15) with
# δ_l = 0. Jexpresso integrates the CONSERVATIVE CARTESIAN form of the same
# equations, q = [φ, φu, φv, φw] (see problems/ShallowWater/SWsphere), so both
# terms have to be translated into that form. That translation is this file.
#
#
# THE FORCING, IN PHYSICAL SPACE
# ------------------------------
# A vorticity forcing F is the curl of a velocity forcing. Writing
#
#   ψ_F  with  ∇ₛ²ψ_F = F ,        u_F = n̂ × ∇ₛψ_F ,
#
# u_F is the non-divergent velocity whose vorticity is F, and forcing the
# momentum with φu_F
#
#   ∂(φu)/∂t = … + φ u_F
#
# forces the vorticity of the flow with F and its divergence with nothing —
# exactly Eq. (9a)-(9b) with the forcing in the vorticity equation only.
#
# ψ_F is built directly from the paper's spectral band: for every total
# wavenumber n with |n - n_f| ≤ Δn/2 and every zonal wavenumber m = 0..n,
#
#   ψ_F(λ,φ) = Σ_{n,m} [ c_{nm}  Y_{nm}^c(λ,φ) + s_{nm} Y_{nm}^s(λ,φ) ]
#
# with Y^c = P̄_{nm}(sin φ) cos mλ /√(4π) and Y^s = P̄_{nm}(sin φ) sin mλ /√(4π)
# the REAL, ORTHONORMAL spherical harmonics (∫ Y² dΩ = 1 over the unit sphere;
# P̄ are the fully normalised associated Legendre functions, computed by the
# standard stable recurrence below). The paper's "random phase e^{iθ}" for the
# complex coefficient f̂_{mn} is (c, s) = A (cos θ, sin θ) here, so every
# (n, m) pair carries a FIXED amount of energy and a random orientation.
#
# The velocity u_F = n̂ × ∇ₛψ_F is NOT formed analytically: ψ_F is evaluated
# at the nodes and differentiated with the SAME spectral-element surface
# gradient the RHS uses (the metric terms a¹, a² of sphere_metrics.jl), then
# made continuous with the mass-weighted direct stiffness average. This keeps
# the forcing exactly tangential, exactly non-divergent in the discrete sense
# the solver understands, and free of the cos φ singularities at the poles
# that the analytic (u, v) components of a spherical harmonic carry — the
# cubed sphere puts a node on each pole.
#
#
# THE AMPLITUDE: A PRESCRIBED ENERGY INPUT RATE
# ---------------------------------------------
# The paper fixes the rate ε₀ at which the forcing injects kinetic energy,
# Eq. (12). For a white-in-time forcing held fixed over one step Δt, one step
# of ∂u/∂t = u_F changes ½|u|² by Δt u·u_F + ½Δt²|u_F|²; the first term
# averages to zero (u_F is independent of u), so the mean injection rate is
#
#   ε = ½ Δt ⟨|u_F|²⟩            (area mean, per unit mass)
#
# — which is why the amplitude in Eq. (12) scales as 1/√Δt. For one (n, m)
# pair of orthonormal harmonics with amplitude A,
#
#   ⟨|u_F|²⟩ = A² n(n+1) / (4π a²)
#
# (∫|∇Y|² = n(n+1)∫Y² on the unit sphere), so the pair injects
#
#   ε_{nm} = Δt A² n(n+1) / (8π a²) .
#
# Following Eq. (12), every total wavenumber n in the band receives the same
# share ε₀/n_band, spread evenly over its 2n+1 (complex) zonal modes — m = 0
# gets one share of 2n+1, each m ≥ 1 pair gets two — and A is solved for.
#
# MARKOVIAN FORCING, Eq. (13). With a decorrelation time τ > 0 the complex
# coefficient z = c + is evolves as the AR(1) process
#
#   z(t+Δt) = r z(t) + √(1-r²) A e^{iθ} ,     r = exp(-Δt/τ)
#
# (the paper's c_r = Δt/(1-r), i.e. r = 1 - Δt/c_r, to first order in Δt/τ).
# Its stationary variance is A², and the injection rate becomes
#
#   ε = Δt ⟨|u_F|²⟩ (1+r) / (2(1-r))
#
# — the correlated forcing keeps pushing the velocity it created — which the
# amplitude below accounts for. τ = 0 is the white-in-time (δ-correlated) case.
#
# THAT LINEAR-RESPONSE AMPLITUDE IS NOT ENOUGH ON A ROTATING SHELL. It assumes
# the velocity the forcing creates persists and stays aligned with the forcing
# for a time τ. It does not: Coriolis turns the velocity vector round in half a
# rotation, and a forced scale larger than L_D sheds most of what it was given
# to gravity waves during geostrophic adjustment. Measured on Jupiter with
# τ = 10 rotations, the fixed-amplitude Markovian forcing injected 0.19 ε₀.
# The paper meets exactly this with Eq. (14),
#
#   f̂_mn = ε₀ f̃_mn / ⟨f̃_mn ζ̃*_mn⟩ ,
#
# a renormalisation of the whole forcing field, every step, by the measured
# correlation between the forcing and the flow, so that the injection IS ε₀.
# That is what :forcing_normalize => true (the default) does, in the discrete
# energy budget of one step: with u_F = g ũ_F the kinetic energy per unit mass
# changes over the step by Δt [ g A + ½Δt g² B ], where
#
#   A = ⟨φ u·ũ_F⟩/⟨φ⟩ ,     B = ⟨φ |ũ_F|²⟩/⟨φ⟩ ,
#
# and g is the positive root of g A + ½Δt g² B = ε₀. The Itô term ½Δt g² B is
# what makes this regular where Eq. (14) is singular: from rest (A = 0) it gives
# the white-noise amplitude instead of dividing by zero, and it keeps g finite
# and positive when the flow happens to oppose the forcing. g is capped at ten
# times the white-noise gain √(2ε₀/(ΔtB)) as a safeguard (reported as `gain`).
# With :forcing_normalize => false the fixed amplitude above is used and the
# injection is whatever the dynamics make of it. Either way the injection rate
# actually produced is measured every update and reported as ε_inj/ε₀.
#
#
# LARGE-SCALE DISSIPATION, Eq. (15) with δ_l = 0
# ----------------------------------------------
#   Rayleigh friction     D_ζ = -ν_r ζ, D_δ = -ν_r δ  ⇒  ∂(φu)/∂t += -ν_r φu
#   Radiative relaxation  D_h = -ν_h h'               ⇒  ∂φ/∂t    += -ν_h (φ - φ_ref)
#
# with φ_ref the reference (rest) geopotential the case stores in q.qe. Both
# are pointwise and both may be on at once. The paper's third option,
# HYPODIFFUSION (-Δ)⁻¹, is a global inverse Laplacian and is not available in
# this element-local RHS.
#
# The SMALL-SCALE dissipation of the paper (∇⁸ hyperdiffusion) is likewise not
# reproduced: the shell already carries the artificial diffusion ν∇ₛ² and the
# modal filter of sphere_rhs.jl, and those are what the case decks use.
#
#
# WHERE IT ENTERS
# ---------------
#   * build_sphere_params (sphere_rhs.jl) builds the forcing from the deck;
#   * sphere_rhs! adds  φu_F - ν_r φu  and  -ν_h(φ - φ_ref)  to the assembled
#     RHS, node by node, after the M⁻¹ scaling (a nodal source under a diagonal
#     mass matrix is exactly that);
#   * the step limiter of sphere_time_loop.jl draws the next step's forcing.
#
# The random draw uses a private SplitMix64 generator, seeded from the deck, so
# a run is reproducible across Julia versions and platforms, and — because
# every rank draws the same sequence — identical on every MPI rank with no
# communication. The basis Y is npoin × nmodes and built once.
#
#
# DECK SWITCHES
# -------------
#   :lsphere_forcing        => false   master switch
#   :forcing_epsilon        => ε₀      energy input rate [m²/s³]
#   :forcing_nf             => n_f     centre of the forced band (total wavenumber)
#   :forcing_dn             => Δn      band width: |n - n_f| ≤ Δn/2 (default 4)
#   :forcing_tau            => τ       decorrelation time [s]; 0 = white in time
#   :forcing_seed           => 1234    seed of the random draw
#   :forcing_normalize      => true    renormalise every step so the injection is ε₀ (Eq. 14)
#   :rayleigh_friction      => ν_r     [1/s], 0 = off
#   :radiative_relaxation   => ν_h     [1/s], 0 = off
#   :sphere_Omega           => Ω       [1/s], only for the Rossby-number diagnostic
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

export St_sphere_forcing
export build_sphere_forcing
export sphere_forcing_update!
export sphere_forcing_apply!
export sphere_forcing_diagnostics
export sphere_zonal_mean
export sphere_harmonic_basis


#---------------------------------------------------------------------------------
# SplitMix64: a tiny, dependency-free, platform-independent generator. Good
# enough for a random phase, and — unlike the default RNG — guaranteed to give
# the same stream on every Julia version, which is what makes a forced run
# reproducible against a stored reference.
#---------------------------------------------------------------------------------
@inline function _splitmix64(state::UInt64)
    state += 0x9e3779b97f4a7c15
    z = state
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    z = z ⊻ (z >> 31)
    return state, z
end

# uniform in [0, 1)
@inline function _rand_uniform(state::UInt64)
    state, z = _splitmix64(state)
    return state, Float64(z >> 11) * (1.0/9007199254740992.0)
end


#---------------------------------------------------------------------------------
# Fully normalised associated Legendre functions P̄_{nm}(μ), μ = sin φ, for all
# n ≤ nmax and m ≤ n, by the standard forward recurrence (e.g. Holmes &
# Featherstone 2002, "geodesy" normalisation):
#
#   P̄_00 = 1
#   P̄_mm = √((2m+1)/(2m)) cos φ P̄_{m-1,m-1}           (P̄_11 = √3 cos φ)
#   P̄_nm = a_nm μ P̄_{n-1,m} - b_nm P̄_{n-2,m}          n > m, P̄_{m-1,m} ≡ 0
#
#   a_nm = √( (2n-1)(2n+1) / ((n-m)(n+m)) )
#   b_nm = √( (2n+1)(n+m-1)(n-m-1) / ((n-m)(n+m)(2n-3)) )
#
# normalised so that the area MEAN of (P̄_nm cos mλ)² over the sphere is 1 for
# every (n, m) — hence the /√(4π) in the harmonics, which makes them
# orthonormal, ∫ Y² dΩ = 1. The Condon-Shortley phase is omitted: only the sign
# of a basis function, irrelevant under a random phase.
#
# Returns P̄ as an (nmax+1) × (nmax+1) matrix indexed [n+1, m+1].
#---------------------------------------------------------------------------------
function _legendre_normalized!(P::AbstractMatrix{TF}, μ::TF, nmax::Int) where {TF}

    cφ = sqrt(max(zero(TF), one(TF) - μ*μ))

    fill!(P, zero(TF))
    P[1, 1] = one(TF)

    # sectoral: P̄_mm. The m = 1 start is √3, NOT the general √((2m+1)/2m) = √(3/2):
    # the fully normalised P̄_nm carry a √2 for every m ≥ 1 (cos mλ averages to
    # ½, not 1), and that factor enters the recurrence exactly here.
    nmax >= 1 && (P[2, 2] = sqrt(TF(3)) * cφ)
    for m = 2:nmax
        P[m+1, m+1] = sqrt(TF(2m+1)/TF(2m)) * cφ * P[m, m]
    end

    # the rest of each order m, going up in degree n
    for m = 0:nmax
        for n = m+1:nmax
            anm = sqrt(TF((2n-1)*(2n+1)) / TF((n-m)*(n+m)))
            bnm = n-2 >= m ?
                  sqrt(TF((2n+1)*(n+m-1)*(n-m-1)) / TF((n-m)*(n+m)*(2n-3))) : zero(TF)
            P[n+1, m+1] = anm*μ*P[n, m+1] - bnm*P[max(n-1, 1), m+1]
        end
    end
    return P
end


#---------------------------------------------------------------------------------
# The real orthonormal spherical harmonics of every (n, m) with n in `ns`,
# evaluated at the nodes of `mesh`.
#
# Returns (Y, pairs): Y is npoin × nmodes, and pairs is a vector of named
# tuples (n, m, kc, ks) giving, for each (n, m), the column kc of the cos
# harmonic and the column ks of the sin harmonic (ks = 0 for m = 0).
#---------------------------------------------------------------------------------
function sphere_harmonic_basis(mesh::St_mesh, ns::AbstractVector{Int}; TF = TFloat)

    npoin = Int(mesh.npoin)
    nmax  = maximum(ns)
    minimum(ns) >= 1 || error(" # ERROR sphere_forcing.jl: the forced band must have n ≥ 1.")

    # the (n, m) pairs and their columns
    pairs = NamedTuple{(:n, :m, :kc, :ks), NTuple{4, Int}}[]
    k = 0
    for n in ns, m = 0:n
        kc = k + 1
        ks = m == 0 ? 0 : k + 2
        push!(pairs, (n = n, m = m, kc = kc, ks = ks))
        k += m == 0 ? 1 : 2
    end
    nmodes = k

    Y   = zeros(TF, npoin, nmodes)
    P   = zeros(TF, nmax+1, nmax+1)
    crd = mesh.coords
    inv4π = one(TF)/sqrt(TF(4π))

    @inbounds for ip = 1:npoin
        x, y, z = crd[1, ip], crd[2, ip], crd[3, ip]
        r = sqrt(x*x + y*y + z*z)
        μ = clamp(z/r, -one(TF), one(TF))      # sin(latitude)
        λ = atan(y, x)

        _legendre_normalized!(P, μ, nmax)

        for p in pairs
            Pnm = P[p.n+1, p.m+1]*inv4π
            if p.m == 0
                Y[ip, p.kc] = Pnm
            else
                sm, cm = sincos(p.m*λ)
                Y[ip, p.kc] = Pnm*cm
                Y[ip, p.ks] = Pnm*sm
            end
        end
    end

    return Y, pairs
end


#---------------------------------------------------------------------------------
# The forcing object. Mutable for the per-step state (coefficients, RNG,
# running diagnostics); the arrays are allocated once.
#---------------------------------------------------------------------------------
mutable struct St_sphere_forcing{TF}

    #--- the band
    nf     ::Int
    dn     ::Int
    nband  ::Int                # number of total wavenumbers in the band
    n      ::Vector{Int}        # per (n,m) pair
    m      ::Vector{Int}
    kc     ::Vector{Int}        # column of the cos harmonic
    ks     ::Vector{Int}        # column of the sin harmonic (0 for m = 0)
    w      ::Vector{TF}         # energy share of the pair, Σw = 1
    A      ::Vector{TF}         # amplitude of the pair for the current Δt
    zc     ::Vector{TF}         # AR(1) state, real part   (= c_nm)
    zs     ::Vector{TF}         # AR(1) state, imag part   (= s_nm)

    #--- the basis and the fields
    Y      ::Matrix{TF}         # npoin × nmodes
    c      ::Vector{TF}         # nmodes  coefficients of ψ_F
    ψ      ::Vector{TF}         # npoin   ψ_F
    F      ::Vector{TF}         # npoin   vorticity forcing ∇ₛ²ψ_F = -Σ n(n+1)/a² c Y
    uF     ::Matrix{TF}         # npoin × neqs : [0, u_F] the Cartesian forcing velocity

    #--- parameters
    ε      ::TF                 # energy input rate [m²/s³]
    τ      ::TF                 # decorrelation time [s], 0 = white
    a      ::TF                 # sphere radius [m]
    Ω      ::TF                 # rotation rate [1/s], diagnostics only (0 = unknown)
    νr     ::TF                 # Rayleigh friction [1/s]
    νh     ::TF                 # radiative relaxation [1/s]
    seed   ::UInt64
    state  ::UInt64             # SplitMix64 state
    lnorm  ::Bool               # per-step renormalisation to ε₀ (Eq. 14)
    gain   ::TF                 # the g of the last update (1 when lnorm is off)
    Δt     ::TF                 # the step A and r were computed for
    r      ::TF                 # AR(1) coefficient for that step
    started::Bool               # has the AR(1) state been drawn

    #--- diagnostics
    Pinj   ::TF                 # injection rate at the last update [m²/s³]
    Psum   ::TF                 # Σ Pinj, for the running mean
    nupd   ::Int                # number of updates
end


function build_sphere_forcing(mesh::St_mesh, metrics::St_sphere_metrics, inputs;
                              neqs = 4, TF = TFloat, verbose = true)

    get(inputs, :lsphere_forcing, false) == true || return nothing

    neqs >= 4 || error(" # ERROR sphere_forcing.jl: the shell forcing needs the 4-equation state [φ, φu, φv, φw].")

    ε  = TF(get(inputs, :forcing_epsilon, 0.0))
    ε >= zero(TF) || error(" # ERROR sphere_forcing.jl: :forcing_epsilon must be ≥ 0.")
    nf = Int(get(inputs, :forcing_nf, 0))
    dn = Int(get(inputs, :forcing_dn, 4))
    τ  = TF(get(inputs, :forcing_tau, 0.0))
    τ >= zero(TF) || error(" # ERROR sphere_forcing.jl: :forcing_tau must be ≥ 0 (0 = white in time).")
    νr = TF(get(inputs, :rayleigh_friction, 0.0))
    νh = TF(get(inputs, :radiative_relaxation, 0.0))
    νr >= zero(TF) && νh >= zero(TF) ||
        error(" # ERROR sphere_forcing.jl: :rayleigh_friction and :radiative_relaxation must be ≥ 0.")
    Ω  = TF(get(inputs, :sphere_Omega, 0.0))
    seed = UInt64(get(inputs, :forcing_seed, 1234))
    lnorm = get(inputs, :forcing_normalize, true) == true

    nf >= 1 || error(" # ERROR sphere_forcing.jl: :forcing_nf (the centre of the forced band) must be ≥ 1.")
    dn >= 0 || error(" # ERROR sphere_forcing.jl: :forcing_dn must be ≥ 0.")

    nlo = max(1, nf - dn ÷ 2)
    nhi = nf + dn ÷ 2
    ns  = collect(nlo:nhi)

    a = TF(mesh.radius)
    a > zero(TF) || error(" # ERROR sphere_forcing.jl: mesh.radius is not set; the forcing needs the shell radius.")

    Y, pairs = sphere_harmonic_basis(mesh, ns; TF = TF)

    npair = length(pairs)
    n  = [p.n  for p in pairs]
    m  = [p.m  for p in pairs]
    kc = [p.kc for p in pairs]
    ks = [p.ks for p in pairs]

    # Eq. (12): each n gets ε₀/n_band, spread over its 2n+1 complex modes; the
    # real pair (cos, sin) of an m ≥ 1 stands for m and -m, i.e. two of them.
    nband = length(ns)
    w = [TF(m[p] == 0 ? 1 : 2)/(TF(2n[p]+1)*TF(nband)) for p = 1:npair]
    @assert abs(sum(w) - 1) < 1e-12

    npoin = Int(mesh.npoin)

    fc = St_sphere_forcing{TF}(nf, dn, nband, n, m, kc, ks, w,
                               zeros(TF, npair), zeros(TF, npair), zeros(TF, npair),
                               Y, zeros(TF, size(Y, 2)),
                               zeros(TF, npoin), zeros(TF, npoin), zeros(TF, npoin, neqs),
                               ε, τ, a, Ω, νr, νh, seed, seed, lnorm, one(TF),
                               zero(TF), zero(TF), false,
                               zero(TF), zero(TF), 0)

    if verbose
        println(" # SHELL FORCING (Scott & Polvani 2007) ..........................")
        @printf(" #   energy input ε₀ = %.4e m²/s³ into n ∈ [%d, %d] (n_f = %d, Δn = %d): %d (n,m) pairs, %d real modes\n",
                ε, nlo, nhi, nf, dn, npair, size(Y, 2))
        if τ > 0
            @printf(" #   Markovian in time, decorrelation time τ = %.4e s\n", τ)
        else
            println(" #   white in time (δ-correlated)")
        end
        νr > 0 && @printf(" #   Rayleigh friction   ν_r = %.4e 1/s on the momentum\n", νr)
        νh > 0 && @printf(" #   radiative relaxation ν_h = %.4e 1/s on φ - φ_ref\n", νh)
        νr > 0 || νh > 0 || println(" #   NO large-scale dissipation: the energy will grow without bound")
        println(" #   amplitude: ", lnorm ? "renormalised every step to inject exactly ε₀ (paper Eq. 14)" :
                                          "fixed (linear-response formula); the injection is reported, not enforced")
        @printf(" #   seed = %d\n", seed)
        println(" # SHELL FORCING ................................................ DONE")
    end

    return fc
end


#
# Amplitudes and AR(1) coefficient for a given step. Called on the first update
# and whenever Δt changes (it never does on the fixed-step shell loop, but the
# cost is nil).
#
function _forcing_set_step!(fc::St_sphere_forcing{TF}, Δt::TF) where {TF}
    Δt > zero(TF) || error(" # ERROR sphere_forcing.jl: Δt must be positive.")
    r = fc.τ > zero(TF) ? exp(-Δt/fc.τ) : zero(TF)
    @inbounds for p = 1:length(fc.n)
        nn = TF(fc.n[p]*(fc.n[p]+1))
        # ε_p = Δt A² n(n+1)/(8πa²) · (1+r)/(1-r)   ⇒   A
        fc.A[p] = sqrt(8π*fc.a*fc.a*fc.ε*fc.w[p]*(one(TF) - r) / (Δt*(one(TF) + r)*nn))
    end
    fc.Δt = Δt
    fc.r  = r
    return fc
end


#---------------------------------------------------------------------------------
# sphere_forcing_update!(fc, q, Δt, mesh, metrics, sp)
#
# Draw the forcing for the NEXT step: new random phases, the AR(1) update of the
# coefficients, ψ_F and its vorticity at the nodes, and u_F = n̂ × ∇ₛψ_F by the
# spectral-element surface gradient + direct stiffness average. `q` is the state
# at the START of the step the forcing will act on; it is only read, for the
# injection-rate diagnostic.
#---------------------------------------------------------------------------------
function sphere_forcing_update!(fc::St_sphere_forcing{TF}, q, Δt,
                                mesh::St_mesh, metrics::St_sphere_metrics,
                                sp) where {TF}

    Δt = TF(Δt)
    (fc.started && fc.Δt == Δt) || _forcing_set_step!(fc, Δt)

    #--- coefficients
    state = fc.state
    r     = fc.r
    s1r2  = sqrt(max(zero(TF), one(TF) - r*r))
    @inbounds for p = 1:length(fc.n)
        A = fc.A[p]
        if fc.m[p] == 0
            state, u = _rand_uniform(state)
            ξc = u < 0.5 ? -A : A          # random sign
            ξs = zero(TF)
        else
            state, u = _rand_uniform(state)
            θ  = TF(2π)*u                  # random phase
            ξc = A*cos(θ)
            ξs = A*sin(θ)
        end
        if fc.started
            fc.zc[p] = r*fc.zc[p] + s1r2*ξc
            fc.zs[p] = r*fc.zs[p] + s1r2*ξs
        else
            fc.zc[p] = ξc                  # start from the stationary distribution
            fc.zs[p] = ξs
        end
        fc.c[fc.kc[p]] = fc.zc[p]
        fc.ks[p] > 0 && (fc.c[fc.ks[p]] = fc.zs[p])
    end
    fc.state   = state
    fc.started = true

    #--- ψ_F and F = ∇ₛ²ψ_F at the nodes
    npoin = Int(mesh.npoin)
    _forcing_fields!(fc.ψ, fc.F, fc.Y, fc.c, fc.n, fc.m, fc.kc, fc.ks, fc.a, npoin)

    #--- u_F = n̂ × ∇ₛψ_F, element by element, then DSS + M⁻¹
    _forcing_velocity_kernel!(fc.uF, fc.ψ, mesh.connijk,
                              Int(mesh.nelem), Int(mesh.ngl),
                              metrics.dξdx, metrics.dξdy, metrics.dξdz,
                              metrics.dηdx, metrics.dηdy, metrics.dηdz,
                              metrics.nx, metrics.ny, metrics.nz,
                              metrics.Je, metrics.dψ, metrics.ω, Int(sp.neqs))
    _sphere_dss_scale!(fc.uF, metrics.Minv, npoin, Int(sp.neqs), sp.dss)

    #--- the injection rate this forcing would produce over the coming step,
    #    Δt(A + ½ΔtB), and — Eq. (14) — the gain g that makes it exactly ε₀
    A, B = _forcing_correlation(q, fc.uF, metrics.M, mesh.gip2owner,
                                MPI.Comm_rank(get_mpi_comm()), npoin)
    g = one(TF)
    if fc.lnorm && B > zero(TF)
        g0 = sqrt(2fc.ε/(Δt*B))                        # white-noise gain: ½Δt g0² B = ε₀
        g  = (-A + sqrt(A*A + 2Δt*B*fc.ε))/(Δt*B)       # positive root of gA + ½Δtg²B = ε₀
        g  = min(g, 10g0)
        if g != one(TF)
            @inbounds for ip = 1:npoin
                fc.uF[ip,2] *= g; fc.uF[ip,3] *= g; fc.uF[ip,4] *= g
                fc.ψ[ip]    *= g; fc.F[ip]    *= g
            end
        end
    end
    fc.gain  = g
    fc.Pinj  = g*A + TF(0.5)*Δt*g*g*B
    fc.Psum += fc.Pinj
    fc.nupd += 1

    return fc
end


function _forcing_fields!(ψ::Vector{TF}, F::Vector{TF}, Y::Matrix{TF}, c::Vector{TF},
                          n, m, kc, ks, a::TF, npoin::Int) where {TF}
    fill!(ψ, zero(TF))
    fill!(F, zero(TF))
    @inbounds for p = 1:length(n)
        λn = -TF(n[p]*(n[p]+1))/(a*a)      # eigenvalue of ∇ₛ² on Y_n
        k  = kc[p]
        ck = c[k]
        for ip = 1:npoin
            ψ[ip] += ck*Y[ip, k]
            F[ip] += λn*ck*Y[ip, k]
        end
        if ks[p] > 0
            k  = ks[p]
            ck = c[k]
            for ip = 1:npoin
                ψ[ip] += ck*Y[ip, k]
                F[ip] += λn*ck*Y[ip, k]
            end
        end
    end
    return ψ
end


#
# n̂ × ∇ₛψ at every node of every element, mass-weighted and accumulated into
# uF[:, 2:4]; column 1 (and any beyond 4) stays zero. The caller completes the
# assembly with _sphere_dss_scale!.
#
function _forcing_velocity_kernel!(uF::AbstractMatrix{TF}, ψ::Vector{TF},
                                   connijk, nelem::Int, ngl::Int,
                                   dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                                   nx, ny, nz, Je, dψ, ω, neqs::Int) where {TF}

    fill!(uF, zero(TF))

    @inbounds for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl

            dξ = dη = zero(TF)
            for k = 1:ngl
                dξ += dψ[k,i]*ψ[connijk[iel,k,j]]
                dη += dψ[k,j]*ψ[connijk[iel,i,k]]
            end

            # ∇ₛψ = a¹ ψ_ξ + a² ψ_η, tangential by construction
            gx = dξ*dξdx[iel,i,j] + dη*dηdx[iel,i,j]
            gy = dξ*dξdy[iel,i,j] + dη*dηdy[iel,i,j]
            gz = dξ*dξdz[iel,i,j] + dη*dηdz[iel,i,j]

            nxx, nyy, nzz = nx[iel,i,j], ny[iel,i,j], nz[iel,i,j]

            # n̂ × ∇ₛψ
            ux = nyy*gz - nzz*gy
            uy = nzz*gx - nxx*gz
            uz = nxx*gy - nyy*gx

            wJ = ω[i]*ω[j]*Je[iel,i,j]
            ip = connijk[iel,i,j]
            uF[ip, 2] += wJ*ux
            uF[ip, 3] += wJ*uy
            uF[ip, 4] += wJ*uz
        end
    end
    return uF
end


#
# The two area means that fix the kinetic energy injection of one step,
#
#   A = ⟨φ u·u_F⟩/⟨φ⟩        the correlation of the flow with the forcing
#   B = ⟨φ |u_F|²⟩/⟨φ⟩       the forcing's own energy
#
# so that the step injects Δt(A + ½ΔtB) per unit mass; the second, Itô-type
# term is the step's own contribution and is the WHOLE injection for
# white-in-time forcing. Owned nodes only, reduced across ranks.
#
function _forcing_correlation(q::AbstractMatrix{TF}, uF::AbstractMatrix{TF},
                              M::AbstractVector{TF}, gip2owner, rank::Int,
                              npoin::Int) where {TF}
    s1 = s2 = sm = zero(TF)
    @inbounds for ip = 1:npoin
        gip2owner[ip] == rank || continue
        φ  = q[ip,1]
        ux, uy, uz = q[ip,2]/φ, q[ip,3]/φ, q[ip,4]/φ
        fx, fy, fz = uF[ip,2], uF[ip,3], uF[ip,4]
        s1 += M[ip]*φ*(ux*fx + uy*fy + uz*fz)
        s2 += M[ip]*φ*(fx*fx + fy*fy + fz*fz)
        sm += M[ip]*φ
    end
    comm = get_mpi_comm()
    if MPI.Comm_size(comm) > 1
        s1 = MPI.Allreduce(s1, MPI.SUM, comm)
        s2 = MPI.Allreduce(s2, MPI.SUM, comm)
        sm = MPI.Allreduce(sm, MPI.SUM, comm)
    end
    return s1/sm, s2/sm
end


#---------------------------------------------------------------------------------
# sphere_forcing_apply!(RHS, q, qe, fc, npoin)
#
# Add the forcing and the large-scale dissipation to the ASSEMBLED RHS:
#
#   RHS[φu] += φ u_F - ν_r φu
#   RHS[φ]  += -ν_h (φ - φ_ref)
#
# Node by node, after the M⁻¹ scaling, which is exactly what a pointwise source
# under the diagonal mass matrix amounts to.
#---------------------------------------------------------------------------------
function sphere_forcing_apply!(RHS::AbstractMatrix{TF}, q, qe,
                               fc::St_sphere_forcing{TF}, npoin::Int) where {TF}
    uF = fc.uF
    νr = fc.νr
    νh = fc.νh
    @inbounds for ip = 1:npoin
        φ = q[ip,1]
        RHS[ip,2] += φ*uF[ip,2] - νr*q[ip,2]
        RHS[ip,3] += φ*uF[ip,3] - νr*q[ip,3]
        RHS[ip,4] += φ*uF[ip,4] - νr*q[ip,4]
        RHS[ip,1] -= νh*(φ - qe[ip,1])
    end
    return RHS
end

# nothing to add when the deck has no forcing
sphere_forcing_apply!(RHS, q, qe, ::Nothing, npoin::Int) = RHS


#---------------------------------------------------------------------------------
# Flow diagnostics for a forced run: kinetic energy per unit mass
# ⟨½φ|u|²⟩/⟨φ⟩, the velocity scale U = √(2 KE), the Rossby number U/(2aΩ), and
# the mean injection rate since the start relative to ε₀.
#---------------------------------------------------------------------------------
function sphere_forcing_diagnostics(q::AbstractMatrix{TF}, fc::St_sphere_forcing{TF},
                                    mesh::St_mesh, metrics::St_sphere_metrics) where {TF}

    comm  = get_mpi_comm()
    rank  = MPI.Comm_rank(comm)
    npoin = Int(mesh.npoin)

    ke, sm = _forcing_ke(q, metrics.M, mesh.gip2owner, rank, npoin)
    if MPI.Comm_size(comm) > 1
        ke = MPI.Allreduce(ke, MPI.SUM, comm)
        sm = MPI.Allreduce(sm, MPI.SUM, comm)
    end

    KE   = ke/sm
    U    = sqrt(2KE)
    Ro   = fc.Ω > 0 ? U/(2fc.a*fc.Ω) : TF(NaN)
    Pbar = fc.nupd > 0 ? fc.Psum/fc.nupd : zero(TF)
    rel  = fc.ε > 0 ? Pbar/fc.ε : TF(NaN)

    return (KE = KE, U = U, Ro = Ro, Pinj = fc.Pinj, Pmean = Pbar, Pmean_rel = rel, gain = fc.gain)
end

function _forcing_ke(q::AbstractMatrix{TF}, M::AbstractVector{TF}, gip2owner,
                     rank::Int, npoin::Int) where {TF}
    ke = sm = zero(TF)
    @inbounds for ip = 1:npoin
        gip2owner[ip] == rank || continue
        φ  = q[ip,1]
        ke += M[ip]*TF(0.5)*(q[ip,2]^2 + q[ip,3]^2 + q[ip,4]^2)/φ
        sm += M[ip]*φ
    end
    return ke, sm
end


#---------------------------------------------------------------------------------
# sphere_zonal_mean(q, mesh, metrics; nbins)
#
# The ZONAL MEAN of the zonal velocity, the meridional velocity and the
# geopotential in `nbins` latitude bands — the paper's main diagnostic (its
# Figs. 2, 4, 6, 7, 13 are all zonal-mean zonal velocity against latitude).
# Mass-weighted (M) averages over the owned nodes of each band, reduced across
# ranks. Returns (lat_centres [deg], ū, v̄, φ̄, count).
#
# u_zonal = u·e_λ with e_λ = (-sin λ, cos λ, 0), formed from the Cartesian
# momentum so it does not depend on any case's output ordering.
#---------------------------------------------------------------------------------
function sphere_zonal_mean(q::AbstractMatrix{TF}, mesh::St_mesh,
                           metrics::St_sphere_metrics; nbins::Int = 90) where {TF}

    comm  = get_mpi_comm()
    rank  = MPI.Comm_rank(comm)
    npoin = Int(mesh.npoin)

    su = zeros(TF, nbins); sv = zeros(TF, nbins); sφ = zeros(TF, nbins)
    sw = zeros(TF, nbins); cnt = zeros(Int, nbins)

    _zonal_mean_kernel!(su, sv, sφ, sw, cnt, q, mesh.coords, metrics.M,
                        mesh.gip2owner, rank, npoin, nbins)

    if MPI.Comm_size(comm) > 1
        su  = MPI.Allreduce(su,  MPI.SUM, comm)
        sv  = MPI.Allreduce(sv,  MPI.SUM, comm)
        sφ  = MPI.Allreduce(sφ,  MPI.SUM, comm)
        sw  = MPI.Allreduce(sw,  MPI.SUM, comm)
        cnt = MPI.Allreduce(cnt, MPI.SUM, comm)
    end

    lat = [(-90 + 180*(b - 0.5)/nbins) for b = 1:nbins]
    ū  = [sw[b] > 0 ? su[b]/sw[b] : TF(NaN) for b = 1:nbins]
    v̄  = [sw[b] > 0 ? sv[b]/sw[b] : TF(NaN) for b = 1:nbins]
    φ̄  = [sw[b] > 0 ? sφ[b]/sw[b] : TF(NaN) for b = 1:nbins]

    return lat, ū, v̄, φ̄, cnt
end

function _zonal_mean_kernel!(su, sv, sφ, sw, cnt, q::AbstractMatrix{TF},
                             crd::AbstractMatrix{TF}, M::AbstractVector{TF},
                             gip2owner, rank::Int, npoin::Int, nbins::Int) where {TF}
    @inbounds for ip = 1:npoin
        gip2owner[ip] == rank || continue
        x, y, z = crd[1, ip], crd[2, ip], crd[3, ip]
        r  = sqrt(x*x + y*y + z*z)
        φl = asin(clamp(z/r, -one(TF), one(TF)))
        λ  = atan(y, x)
        b  = clamp(1 + floor(Int, (φl + π/2)/π*nbins), 1, nbins)

        φ  = q[ip,1]
        ux, uy, uz = q[ip,2]/φ, q[ip,3]/φ, q[ip,4]/φ
        sλ, cλ = sincos(λ)
        sφl, cφl = sincos(φl)
        uz_ =  -ux*sλ     + uy*cλ
        vm  =  -ux*sφl*cλ - uy*sφl*sλ + uz*cφl

        su[b]  += M[ip]*uz_
        sv[b]  += M[ip]*vm
        sφ[b]  += M[ip]*φ
        sw[b]  += M[ip]
        cnt[b] += 1
    end
    return nothing
end
