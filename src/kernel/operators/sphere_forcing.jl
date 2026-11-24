#---------------------------------------------------------------------------------
# STOCHASTIC SMALL-SCALE FORCING AND LARGE-SCALE DISSIPATION ON THE SHELL
#
#   Scott, R. K. & Polvani, L. M. (2007), "Forced-dissipative shallow-water
#   turbulence on the sphere and the atmospheric circulation of the giant
#   planets", J. Atmos. Sci. 64, 3158-3176.
#
# This is the operator problems/ShallowWater/SWsphere_ScottPolvani asks for with
# :lsphere_forcing => true. Everything else in that case — the equations, the
# cubed sphere, the spectral elements, the Lagrange projection, the filter —
# is SWsphere's; the only thing that makes it a Scott & Polvani run is what is
# in this file, because the initial state is REST and every bit of motion in the
# answer is put there from here.
#
# WHY IT IS NOT IN user_source!. The forcing is ONE random field drawn per time
# step and shared by every node — a global object. user_source! is pointwise: it
# sees a single node's state and nothing else, so it cannot see a field. Hence
# this operator, which is added to the assembled tendency (after the direct
# stiffness summation and the M⁻¹ scaling) rather than integrated against the
# test functions. Adding a nodal tendency to ∂q/∂t is the same thing as adding
# its L² projection to the weak form, because the mass matrix is diagonal here.
#
#---------------------------------------------------------------------------------
# WHAT IS FORCED
#
# A random, isotropic VORTICITY field confined to a narrow band of total
# spherical-harmonic wavenumbers (the paper's Eq. 11, with Δn = 4),
#
#   ζ_F(λ,φ) = Σ_{|n - n_f| ≤ Δn/2} Σ_{m=-n}^{n} c_{nm} Y_n^m(λ,φ) ,
#
# and converted to the velocity the shallow-water momentum equation can carry.
# A vorticity field is realised by the non-divergent velocity u = n̂ × ∇ₛψ with
# ∇ₛ²ψ = ζ, and on the sphere the harmonics are eigenfunctions of the surface
# Laplacian, ∇ₛ²Y_n^m = -n(n+1)/a² Y_n^m, so the streamfunction is had by
# dividing coefficient by coefficient:
#
#   ψ_F = Σ  [ -a²/(n(n+1)) ] c_{nm} Y_n^m ,        u_F = n̂ × ∇ₛψ_F .          (1)
#
# ∇ₛψ_F is taken with the SAME element machinery the vorticity diagnostic uses —
# the contravariant basis a¹, a² of the manifold, assembled and mass-averaged
# into a continuous nodal field — so the forcing lives in the discrete space the
# solution lives in, not in a spectral space beside it. u_F is then exactly
# tangent to the shell, n̂ × g ⊥ n̂ identically, whatever g the discrete gradient
# returned.
#
# The tendency added to the state q = [φ, φu, φv, φw]ᵀ is
#
#   ∂(φu)/∂t += φ u_F  -  ν_r (φu) ,          ∂φ/∂t += -ν_h (φ - φ_ref) .      (2)
#
# The first momentum term is a VELOCITY forcing (∂u/∂t += u_F, which is what a
# vorticity forcing means); the other two are the paper's Eq. (15) large-scale
# dissipation with δ_l = 0 — Rayleigh friction ν_r on the momentum, radiative
# relaxation ν_h of the height towards the resting depth φ_ref = gH, which is
# what initialize() left in qe[:,1]. ν_h is 0 in the shipped decks; note that
# switching it on relaxes the MASS as well, so δmass/mass is no longer a
# conservation check when it is.
#
#---------------------------------------------------------------------------------
# HOW THE COEFFICIENTS EVOLVE
#
# A Markov process of decorrelation time τ (the paper's Eq. 13, c_r = 10
# rotations), advanced once per STEP and held fixed across the RK stages:
#
#   c ← r c + √(1-r²) ξ ,     r = exp(-Δt/τ) ,     ξ ~ N(0,1) i.i.d.           (3)
#
# which is the stationary AR(1) with unit variance and autocorrelation
# exp(-Δt/τ). τ = 0 gives the δ-correlated forcing of Eq. (12): c = ξ, redrawn
# every step. Holding c fixed over the step, rather than redrawing per stage, is
# what makes the five RK stages see one consistent forcing field; redrawing per
# stage would inject a stage-dependent field and the scheme's order would be the
# least of the problems.
#
#---------------------------------------------------------------------------------
# HOW THE AMPLITUDE IS SET  (:forcing_normalize => true, the paper's Eq. 14)
#
# The point of the forcing is to inject energy at a PRESCRIBED rate ε₀, and the
# rate a given field actually injects depends on the flow it is injecting into.
#
# ε₀ IS A SPECIFIC RATE — energy per unit mass per unit time, [m²/s³] — and not
# a total. That is what the deck builds: ε₀ = 2ν_l E/f_up with E = ½U², the mean
# specific kinetic energy, so the quantity it prescribes the growth of is
#
#   Ê = ( ∫ ½φ|u|² dΩ ) / ( ∫ φ dΩ ) ,        dÊ/dt = ε₀ .                     (4)
#
# (φ dΩ is the mass element up to the constant g, which cancels in the ratio.)
# Getting this wrong is not a small error — a rate normalised by the total mass
# instead of per unit mass is out by ~10²⁰ on a giant planet, and the run looks
# like it is still at rest.
#
# Over one step of size Δt the velocity picks up u ← u + α u_F, so
#
#   Δ∫½φ|u|²dΩ = α Δt ∫φ u·u_F dΩ  +  ½ α² Δt² ∫φ|u_F|² dΩ .                  (5)
#
# Dividing by Δt ∫φdΩ and asking for ε₀ gives a quadratic in the amplitude α,
#
#   A α² + B α - ε₀ = 0 ,   A = ½Δt ∫φ|u_F|²dΩ / ∫φdΩ ,  B = ∫φ u·u_F dΩ / ∫φdΩ
#
# with one positive root,
#
#   α = ( -B + √(B² + 4 A ε₀) ) / (2A) ,                                        (6)
#
# which always exists because A > 0 and ε₀ > 0. This is the standard
# construction (Alvelius, Phys. Fluids 11, 1880, 1999) and it is what makes the
# forcing well defined AT REST, where the first-order term B vanishes and the
# whole input is the α² self-term — an ε₀ imposed through the first-order term
# alone would be a division by zero on the first step of this very case.
#
# With :forcing_normalize => false the amplitude is frozen at the α computed on
# the first step and the realised rate is whatever the flow makes of it. The
# deck records the measurement: on Jupiter at τ = 10 rotations that was 0.19 ε₀,
# because Coriolis and geostrophic adjustment decorrelate the flow from a slowly
# varying forcing long before τ. Both rates are reported, so the two can be told
# apart without guessing.
#
#---------------------------------------------------------------------------------
# REPRODUCIBILITY AND PARALLELISM
#
# The coefficients are drawn from a splitmix64 stream with Box-Muller on top,
# written out here rather than taken from Random. Two reasons, both about
# getting the SAME numbers: Base.Random's stream is explicitly not guaranteed
# stable across Julia releases, so a run would stop reproducing on an upgrade;
# and every MPI rank must draw the IDENTICAL coefficient vector, since c_{nm}
# describes one global field. Every rank seeds the same stream and draws the
# same values in the same order, so they stay in lockstep with no communication
# at all — the forcing costs exactly two reduced scalars per step (A and B), not
# a field exchange.
#
# The harmonics themselves are evaluated by the stable fully normalised
# recursion, on the fly, with no npoin × nmode table: at n_f = 24, Δn = 4 the
# table would be 245 columns — 29 MB on the Galewsky grid but 1.2 GB on the
# 64×64 one — while the recursion costs ~380 flops per node per step and works
# at any resolution.
#---------------------------------------------------------------------------------

export St_sphere_forcing
export build_sphere_forcing
export sphere_forcing_step!
export sphere_forcing_apply!
export sphere_forcing_report


#---------------------------------------------------------------------------------
# A reproducible normal stream: splitmix64 (Steele, Lea & Flood 2014) with
# Box-Muller. Deliberately self-contained — see the note above.
#---------------------------------------------------------------------------------
mutable struct St_splitmix64
    s::UInt64
end

@inline function _sm64(rng::St_splitmix64)
    rng.s += 0x9e3779b97f4a7c15
    z = rng.s
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

# uniform on (0,1): 53 mantissa bits, never exactly 0 (log(0) below)
@inline _smrand(rng::St_splitmix64) = (Float64(_sm64(rng) >> 11) + 0.5) * (1.0/9007199254740992.0)

@inline function _smrandn(rng::St_splitmix64)
    u1 = _smrand(rng)
    u2 = _smrand(rng)
    return sqrt(-2.0*log(u1))*cos(2π*u2)
end


#---------------------------------------------------------------------------------
# Everything the forcing needs, allocated once.
#---------------------------------------------------------------------------------
mutable struct St_sphere_forcing{TFloat, TDSS}
    #--- what the deck asked for
    nf      ::Int                 # centre of the forced band
    dn      ::Int                 # its width: |n - n_f| ≤ dn/2
    nlo     ::Int                 # the band, resolved
    nhi     ::Int
    ε₀      ::TFloat              # prescribed energy input rate [m²/s³]
    τ       ::TFloat              # Markov decorrelation time [s] (0 = white)
    Δt      ::TFloat              # the step the amplitude is normalised over
    νr      ::TFloat              # Rayleigh friction on the momentum [1/s]
    νh      ::TFloat              # radiative relaxation of the height [1/s]
    lnorm   ::Bool                # renormalise every step (paper Eq. 14)
    seed    ::UInt64

    #--- the random state
    rng     ::St_splitmix64
    nmode   ::Int
    c       ::Vector{TFloat}      # nmode : the AR(1) coefficients, Eq. (3)
    w       ::Vector{TFloat}      # nmode : -a²/(n(n+1)), folded in once

    #--- recursion tables for the normalised associated Legendre functions
    arec    ::Matrix{TFloat}      # (nhi+1) × (nhi+1) indexed [n+1, m+1]
    brec    ::Matrix{TFloat}
    drec    ::Vector{TFloat}      # nhi+1 : the diagonal step √((2m+1)/(2m))

    #--- fields
    ψ       ::Vector{TFloat}      # npoin : the streamfunction, Eq. (1)
    uF      ::Matrix{TFloat}      # npoin × 3 : the forcing velocity
    acc     ::Matrix{TFloat}      # npoin × neqs : gradient assembly scratch
    dss     ::TDSS               # the neqs-wide assembler cache (nothing in serial)

    #--- diagnostics
    α       ::TFloat              # the amplitude last used, Eq. (5)
    εreal   ::TFloat              # the rate it actually injected [m²/s³]
    εmean   ::TFloat              # its running mean
    nstep   ::Int
    a       ::TFloat              # the shell radius, for the report
end


#---------------------------------------------------------------------------------
# build_sphere_forcing(mesh, metrics, sp, inputs; Δt) -> St_sphere_forcing | nothing
#
# `nothing` when :lsphere_forcing is absent or false, which is what every case
# but Scott & Polvani wants; the caller then pays nothing at all, because the
# apply below is a typed no-op on ::Nothing.
#---------------------------------------------------------------------------------
function build_sphere_forcing(mesh::St_mesh, metrics::St_sphere_metrics,
                              sp::St_sphere_params, inputs;
                              Δt::Real, TF = TFloat, verbose::Bool = true)

    get(inputs, :lsphere_forcing, false) == true || return nothing

    npoin = Int(mesh.npoin)
    neqs  = Int(sp.neqs)
    a     = TF(mesh.radius)

    nf = Int(get(inputs, :forcing_nf, 0))
    dn = Int(get(inputs, :forcing_dn, 0))
    nf > 0 ||
        error(" # ERROR sphere_forcing.jl: :lsphere_forcing => true needs a positive :forcing_nf (the centre of the forced wavenumber band).")
    dn >= 0 ||
        error(" # ERROR sphere_forcing.jl: :forcing_dn must be non-negative.")

    nlo = max(1, nf - dn ÷ 2)          # n = 0 is a constant: no velocity, and 1/(n(n+1)) is a
    nhi = nf + dn ÷ 2                  # division by zero. n = 1 is a solid-body rotation.
    nlo <= nhi ||
        error(" # ERROR sphere_forcing.jl: the forced band is empty; check :forcing_nf and :forcing_dn.")

    ε₀ = TF(get(inputs, :forcing_epsilon, 0.0))
    ε₀ > 0 ||
        error(" # ERROR sphere_forcing.jl: :forcing_epsilon must be positive; a zero energy input rate leaves the fluid at rest, which is not a simulation.")

    τ  = TF(get(inputs, :forcing_tau, 0.0))
    τ >= 0 || error(" # ERROR sphere_forcing.jl: :forcing_tau must be non-negative (0 = δ-correlated).")

    νr = TF(get(inputs, :rayleigh_friction, 0.0))
    νh = TF(get(inputs, :radiative_relaxation, 0.0))
    (νr >= 0 && νh >= 0) ||
        error(" # ERROR sphere_forcing.jl: :rayleigh_friction and :radiative_relaxation must be non-negative; a negative one is an energy source at the largest scales.")

    lnorm = get(inputs, :forcing_normalize, true) == true
    seed  = UInt64(get(inputs, :forcing_seed, 1234))

    #--- the mode list, in the canonical order the kernel walks (m outer, n inner)
    nmode = 0
    for n = nlo:nhi
        nmode += 2n + 1
    end

    #
    # The weight each mode's coefficient is multiplied by, built once. Two
    # factors, and BOTH matter:
    #
    #   -a²/(n(n+1))   the inverse surface Laplacian of Eq. (1), turning a
    #                  vorticity coefficient into a streamfunction one;
    #   √2  (m ≠ 0)    what makes the REAL harmonics orthonormal. The complex
    #                  Y_n^m and Y_n^{-m} carry half the variance each; the real
    #                  cos/sin pair built from them needs √2 to satisfy
    #                  ∫Y²dΩ = 1 over the unit sphere, exactly as the m = 0
    #                  member already does. Leaving it out would give the zonal
    #                  (m = 0) modes twice the variance of every other mode and
    #                  the forcing would not be isotropic — it would have a
    #                  systematic zonal bias, which on a rotating planet is
    #                  precisely the answer this case is asking about.
    #
    w   = zeros(TF, nmode)
    rt2 = sqrt(TF(2))
    idx = 0
    for m = 0:nhi, n = max(nlo, m):nhi
        wk = -a*a/TF(n*(n+1))
        if m == 0
            idx += 1; w[idx] = wk
        else
            idx += 1; w[idx] = rt2*wk               # cos(mλ) member
            idx += 1; w[idx] = rt2*wk               # sin(mλ) member
        end
    end
    idx == nmode ||
        error(string(" # ERROR sphere_forcing.jl: the mode count is inconsistent (", idx, " built, ", nmode, " expected)."))

    #--- the normalised-Legendre recursion coefficients, precomputed
    arec = zeros(TF, nhi+1, nhi+1)
    brec = zeros(TF, nhi+1, nhi+1)
    drec = zeros(TF, nhi+1)
    for m = 1:nhi
        drec[m+1] = sqrt(TF(2m+1)/TF(2m))
    end
    for m = 0:nhi, n = m+2:nhi
        arec[n+1, m+1] = sqrt(TF((2n-1)*(2n+1))/TF((n-m)*(n+m)))
        brec[n+1, m+1] = sqrt(TF((n+m-1)*(n-m-1))/TF((2n-3)*(2n-1)))
    end

    fc = St_sphere_forcing(nf, dn, nlo, nhi, ε₀, τ, TF(Δt), νr, νh, lnorm, seed,
                           St_splitmix64(seed),
                           nmode, zeros(TF, nmode), w,
                           arec, brec, drec,
                           zeros(TF, npoin), zeros(TF, npoin, 3),
                           zeros(TF, npoin, neqs), sp.dss,
                           zero(TF), zero(TF), zero(TF), 0, a)

    #--- the first draw: the stationary distribution, not a ramp from zero
    _forcing_draw!(fc; first = true)

    if verbose && MPI.Comm_rank(get_mpi_comm()) == 0
        T = 2π/Float64(get(inputs, :sp_Omega, get(inputs, :sphere_Omega, 7.292e-5)))
        @printf(" #   forcing: n ∈ [%d, %d] (n_f = %d, Δn = %d), %d modes\n", nlo, nhi, nf, dn, nmode)
        @printf(" #     ε₀ = %.4e m²/s³ ; τ = %.4e s (%.2f rotations, r = %.6f) ; seed = %d\n",
                Float64(ε₀), Float64(τ), Float64(τ)/T, τ > 0 ? exp(-Δt/Float64(τ)) : 0.0, Int(seed))
        @printf(" %s\n", lnorm ? " #     amplitude renormalised every step (paper Eq. 14)" :
                                 " #     amplitude FROZEN at the first step's value; the realised rate will differ from ε₀")
        @printf(" #     dissipation: ν_r = %.4e 1/s (%.3e per rotation) ; ν_h = %.4e 1/s\n",
                Float64(νr), Float64(νr)*T, Float64(νh))
        νh > 0 && println(" #     NOTE ν_h > 0 relaxes the height, so ∫φ dΩ is no longer conserved and δmass/mass stops being a check.")
    end

    return fc
end


#---------------------------------------------------------------------------------
# The AR(1) draw, Eq. (3). Every rank runs this and gets the same numbers.
#---------------------------------------------------------------------------------
function _forcing_draw!(fc::St_sphere_forcing{TF}; first::Bool = false) where {TF}
    c = fc.c
    if fc.τ <= 0 || first
        @inbounds for j in eachindex(c)
            c[j] = TF(_smrandn(fc.rng))
        end
    else
        r = exp(-fc.Δt/fc.τ)
        s = sqrt(max(zero(TF), one(TF) - r*r))
        @inbounds for j in eachindex(c)
            c[j] = r*c[j] + s*TF(_smrandn(fc.rng))
        end
    end
    return fc
end


#---------------------------------------------------------------------------------
# ψ_F at the nodes, Eq. (1). The fully normalised associated Legendre recursion,
# evaluated on the fly; see the header for why there is no table.
#
#   P̄_0^0 = √(1/4π)
#   P̄_m^m = √((2m+1)/2m) sinθ P̄_{m-1}^{m-1}
#   P̄_{m+1}^m = √(2m+3) μ P̄_m^m
#   P̄_n^m = a_n^m ( μ P̄_{n-1}^m - b_n^m P̄_{n-2}^m )
#
# cos(mλ) and sin(mλ) come from the complex-multiply recursion rather than 27
# trig calls per node per step.
#---------------------------------------------------------------------------------
function _forcing_stream_kernel!(ψ::Vector{TF}, crd::AbstractMatrix{TF}, npoin::Int,
                                 nlo::Int, nhi::Int,
                                 c::Vector{TF}, w::Vector{TF},
                                 arec::Matrix{TF}, brec::Matrix{TF},
                                 drec::Vector{TF}) where {TF}

    p00 = sqrt(one(TF)/(TF(4)*TF(π)))

    @inbounds for ip = 1:npoin

        x, y, z = crd[1,ip], crd[2,ip], crd[3,ip]
        r   = sqrt(x*x + y*y + z*z)
        μ   = z/r                                   # sin(latitude)
        sθ  = sqrt(max(zero(TF), one(TF) - μ*μ))    # cos(latitude)

        # cos λ, sin λ from the horizontal projection; at the poles sθ = 0 and
        # every m > 0 harmonic vanishes there anyway, so the value is arbitrary.
        rh = sqrt(x*x + y*y)
        cλ = rh > 0 ? x/rh : one(TF)
        sλ = rh > 0 ? y/rh : zero(TF)

        acc = zero(TF)
        idx = 0

        cm = one(TF)      # cos(0λ)
        sm = zero(TF)     # sin(0λ)
        pmm = p00

        for m = 0:nhi

            if m > 0
                pmm *= drec[m+1]*sθ
                # (cm, sm) ← (cm + i sm)(cλ + i sλ)
                cmn = cm*cλ - sm*sλ
                sm  = sm*cλ + cm*sλ
                cm  = cmn
            end

            # climb n from m to nhi, accumulating once inside the band
            pnm2 = zero(TF)      # P̄_{n-2}^m
            pnm1 = zero(TF)      # P̄_{n-1}^m
            pn   = zero(TF)      # P̄_n^m

            for n = m:nhi
                if n == m
                    pn = pmm
                elseif n == m+1
                    pn = sqrt(TF(2m+3))*μ*pmm
                else
                    pn = arec[n+1,m+1]*(μ*pnm1 - brec[n+1,m+1]*pnm2)
                end

                if n >= nlo
                    idx += 1
                    acc += w[idx]*c[idx]*pn*cm
                    if m > 0
                        idx += 1                       # the sin(mλ) partner shares P̄
                        acc += w[idx]*c[idx]*pn*sm
                    end
                end

                pnm2 = pnm1
                pnm1 = pn
            end
        end

        # The √2 that makes the real m ≠ 0 harmonics orthonormal is already in
        # w[idx] — see build_sphere_forcing — so nothing is owed here.
        ψ[ip] = acc
    end

    return ψ
end


#---------------------------------------------------------------------------------
# ∇ₛψ at the nodes: a¹ ∂ψ/∂ξ + a² ∂ψ/∂η, assembled and mass-averaged exactly the
# way sphere_relative_vorticity! averages the curl, so the result is a
# continuous nodal field rather than a per-element one.
#---------------------------------------------------------------------------------
function _forcing_grad_kernel!(acc::AbstractMatrix{TF}, ψ::Vector{TF},
                               connijk, nelem::Int, ngl::Int,
                               dξdx, dξdy, dξdz, dηdx, dηdy, dηdz,
                               Je, dψ, ω) where {TF}

    fill!(acc, zero(TF))

    @inbounds for iel = 1:nelem
        for j = 1:ngl, i = 1:ngl

            dξ = zero(TF)
            dη = zero(TF)
            for k = 1:ngl
                dξ += dψ[k,i]*ψ[connijk[iel,k,j]]
                dη += dψ[k,j]*ψ[connijk[iel,i,k]]
            end

            a1x, a1y, a1z = dξdx[iel,i,j], dξdy[iel,i,j], dξdz[iel,i,j]
            a2x, a2y, a2z = dηdx[iel,i,j], dηdy[iel,i,j], dηdz[iel,i,j]

            wq = ω[i]*ω[j]*Je[iel,i,j]
            ip = connijk[iel,i,j]

            acc[ip,1] += wq*(a1x*dξ + a2x*dη)
            acc[ip,2] += wq*(a1y*dξ + a2y*dη)
            acc[ip,3] += wq*(a1z*dξ + a2z*dη)
        end
    end

    return acc
end


#
# u_F = n̂ × ∇ₛψ, with n̂ the exact radial unit normal. The cross product is
# identically orthogonal to n̂, so u_F is tangent to the shell to round-off
# whatever the discrete gradient returned — no projection needed.
#
function _forcing_curl_kernel!(uF::Matrix{TF}, acc::AbstractMatrix{TF},
                               crd::AbstractMatrix{TF}, npoin::Int) where {TF}
    @inbounds for ip = 1:npoin
        x, y, z = crd[1,ip], crd[2,ip], crd[3,ip]
        r  = sqrt(x*x + y*y + z*z)
        nx, ny, nz = x/r, y/r, z/r
        gx, gy, gz = acc[ip,1], acc[ip,2], acc[ip,3]
        uF[ip,1] = ny*gz - nz*gy
        uF[ip,2] = nz*gx - nx*gz
        uF[ip,3] = nx*gy - ny*gx
    end
    return uF
end


#
# The three integrals Eq. (5) needs — ∫φ|u_F|²dΩ, ∫φ u·u_F dΩ and the mass
# ∫φ dΩ that turns them into specific rates — over this rank's OWN nodes. Same
# ownership test the conserved integrals use, so a node on a partition seam is
# counted once however many ranks hold it.
#
function _forcing_moments(q::AbstractMatrix{TF}, uF::Matrix{TF}, M::AbstractVector{TF},
                          gip2owner, rank::Int, npoin::Int) where {TF}
    Q = zero(TF)      # ∫φ|u_F|² dΩ
    P = zero(TF)      # ∫φ u·u_F dΩ  =  ∫(φu)·u_F dΩ
    W = zero(TF)      # ∫φ dΩ
    @inbounds for ip = 1:npoin
        gip2owner[ip] == rank || continue
        φ = q[ip,1]
        fx, fy, fz = uF[ip,1], uF[ip,2], uF[ip,3]
        Q += M[ip]*φ*(fx*fx + fy*fy + fz*fz)
        P += M[ip]*(q[ip,2]*fx + q[ip,3]*fy + q[ip,4]*fz)
        W += M[ip]*φ
    end
    return Q, P, W
end


#---------------------------------------------------------------------------------
# sphere_forcing_step!(fc, q, mesh, metrics, sp)
#
# ONCE PER STEP: advance the Markov coefficients, rebuild u_F, and set its
# amplitude so that the step injects ε₀. Called from the step limiter, i.e. on
# the state the step just finished, which is the state the NEXT step starts
# from — so B is evaluated against exactly the u the forcing will act on.
#---------------------------------------------------------------------------------
sphere_forcing_step!(::Nothing, q, mesh, metrics, sp) = nothing

function sphere_forcing_step!(fc::St_sphere_forcing{TF}, q,
                              mesh::St_mesh, metrics::St_sphere_metrics,
                              sp::St_sphere_params) where {TF}

    npoin = Int(mesh.npoin)
    neqs  = Int(sp.neqs)

    _forcing_draw!(fc)

    _forcing_stream_kernel!(fc.ψ, mesh.coords, npoin, fc.nlo, fc.nhi,
                            fc.c, fc.w, fc.arec, fc.brec, fc.drec)

    _forcing_grad_kernel!(fc.acc, fc.ψ, mesh.connijk,
                          Int(mesh.nelem), Int(mesh.ngl),
                          metrics.dξdx, metrics.dξdy, metrics.dξdz,
                          metrics.dηdx, metrics.dηdy, metrics.dηdz,
                          metrics.Je, metrics.dψ, metrics.ω)

    # the gradient is an element-wise sum; complete it across ranks and divide
    # by the mass, exactly as the curl is completed. The cache was built neqs
    # wide, so the whole neqs columns go (4:neqs are zero and stay zero).
    _sphere_dss_scale!(fc.acc, metrics.Minv, npoin, neqs, fc.dss)

    _forcing_curl_kernel!(fc.uF, fc.acc, mesh.coords, npoin)

    #--- the amplitude, Eq. (6)
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    Q, P, W = _forcing_moments(q, fc.uF, metrics.M, mesh.gip2owner, rank, npoin)
    if MPI.Comm_size(comm) > 1
        qpw = TF[Q, P, W]
        MPI.Allreduce!(qpw, MPI.SUM, comm)
        Q, P, W = qpw[1], qpw[2], qpw[3]
    end

    W > 0 ||
        error(" # ERROR sphere_forcing.jl: ∫φ dΩ = 0. The layer has no mass, so a specific energy input rate is meaningless.")

    # per unit mass — see Eq. (4) in the header
    A = TF(0.5)*fc.Δt*Q/W
    B = P/W

    if fc.lnorm || fc.nstep == 0
        A > 0 ||
            error(" # ERROR sphere_forcing.jl: ∫φ|u_F|²dΩ = 0. The forced band produced no velocity — check :forcing_nf against the grid's resolvable wavenumbers.")
        α = (-B + sqrt(B*B + TF(4)*A*fc.ε₀))/(TF(2)*A)
        fc.α = α
    else
        α = fc.α                  # frozen amplitude; the realised rate is measured below
    end

    @inbounds for ip = 1:npoin, k = 1:3
        fc.uF[ip,k] *= α
    end

    # the specific rate this field will actually inject over the step, Eq. (5)
    # with the amplitude folded in. Under renormalisation it is ε₀ by
    # construction, and is kept as a check on the algebra rather than as news;
    # with :forcing_normalize => false it is the number that matters.
    fc.εreal = α*α*A + α*B
    fc.nstep += 1
    fc.εmean += (fc.εreal - fc.εmean)/fc.nstep

    return fc
end


#---------------------------------------------------------------------------------
# sphere_forcing_apply!(RHS, q, qe, fc, npoin)
#
# EVERY STAGE: add the tendency of Eq. (2) to the assembled ∂q/∂t. u_F is held
# fixed across the stages of a step; see the header.
#---------------------------------------------------------------------------------
@inline sphere_forcing_apply!(RHS, q, qe, ::Nothing, npoin::Int) = RHS

function sphere_forcing_apply!(RHS::AbstractMatrix{TF}, q, qe,
                               fc::St_sphere_forcing{TF}, npoin::Int) where {TF}
    uF = fc.uF
    νr = fc.νr
    νh = fc.νh

    @inbounds for ip = 1:npoin
        φ = q[ip,1]
        RHS[ip,2] += φ*uF[ip,1] - νr*q[ip,2]
        RHS[ip,3] += φ*uF[ip,2] - νr*q[ip,3]
        RHS[ip,4] += φ*uF[ip,3] - νr*q[ip,4]
    end

    if νh > 0
        @inbounds for ip = 1:npoin
            RHS[ip,1] -= νh*(q[ip,1] - qe[ip,1])
        end
    end

    return RHS
end


#---------------------------------------------------------------------------------
# One line for the diagnostics stream, and the closing summary.
#---------------------------------------------------------------------------------
sphere_forcing_report(::Nothing) = ""

function sphere_forcing_report(fc::St_sphere_forcing)
    return @sprintf("  ε = %.3e (%.3f ε₀)  α = %.3e",
                    Float64(fc.εreal), Float64(fc.εreal)/Float64(fc.ε₀), Float64(fc.α))
end
