using FFTW: plan_ifft!
using Random: MersenneTwister

# Rogallo synthetic isotropic turbulence on a uniform N³ periodic grid over [-π,π]³.
# Solenoidal velocity with the [41] Yamazaki-Ishihara-Kaneda initial spectrum shape
# E(k) ∝ k⁴ exp(-2(k/kp)²); the field is rescaled to ⟨|u|²⟩ = 1 (3u0²=1, urms=u0=0.577/comp).
#
# The spectral field is built with explicit Hermitian (conjugate) symmetry,
# û(-k) = conj(û(k)), so that a single inverse FFT yields an *exactly* real,
# divergence-free field whose shell-summed energy equals the prescribed E(k).
# This replaces the previous `real(ifft(û))` of a non-Hermitian array, which
# discarded half of the prescribed spectral energy (only masked by the global
# rescale) and corrupted the realized per-mode spectrum.
#
# Allocation: only the three complex spectral arrays and the three real output
# arrays are allocated; the inverse transforms are done in place (plan_ifft!)
# and the variance normalization is a single scalar pass — no N³ temporaries.
function _rogallo_field(N::Int, kp::Float64, seed::Int)
    rng = MersenneTwister(seed)
    uh = zeros(ComplexF64, N, N, N)
    vh = zeros(ComplexF64, N, N, N)
    wh = zeros(ComplexF64, N, N, N)

    Nh = N ÷ 2
    @inline kfreq(i0) = i0 <= Nh ? i0 : i0 - N   # signed wavenumber of 0-based index
    @inline mir(i0)   = (N - i0) % N             # conjugate (mirror) 0-based index
    shapeE(k) = k^4 * exp(-2.0*(k/kp)^2)

    # Visit each (k,-k) conjugate pair once (canonical half via linear-index order),
    # assign the Rogallo-projected coefficient to k and its conjugate to -k.
    @inbounds for c0 in 0:N-1, b0 in 0:N-1, a0 in 0:N-1
        k1 = kfreq(a0); k2 = kfreq(b0); k3 = kfreq(c0)
        ksq = k1*k1 + k2*k2 + k3*k3
        ksq == 0 && continue                     # zero mean flow
        ma0 = mir(a0); mb0 = mir(b0); mc0 = mir(c0)
        L  = a0  + N*(b0  + N*c0)
        Lm = ma0 + N*(mb0 + N*mc0)
        L > Lm && continue                       # conjugate partner already handled

        k     = sqrt(float(ksq))
        amp   = sqrt(shapeE(k)/(4π*ksq))         # per-mode rms: Σ_shell |û|² = E(k)
        θ1 = 2π*rand(rng); θ2 = 2π*rand(rng); φ = 2π*rand(rng)
        α  = amp*cis(θ1)*cos(φ)
        β  = amp*cis(θ2)*sin(φ)
        kperp = sqrt(float(k1*k1 + k2*k2))
        if kperp < 1e-12
            U = α; V = β; W = zero(ComplexF64)
        else
            U = (α*k*k2 + β*k1*k3)/(k*kperp)
            V = (β*k2*k3 - α*k*k1)/(k*kperp)
            W = -β*kperp/k
        end

        a = a0+1; b = b0+1; c = c0+1
        if L == Lm
            # self-conjugate mode (k ≡ -k mod N): coefficient must be real
            uh[a,b,c] = complex(real(U)); vh[a,b,c] = complex(real(V)); wh[a,b,c] = complex(real(W))
        else
            ma = ma0+1; mb = mb0+1; mc = mc0+1
            uh[a,b,c]    = U;       vh[a,b,c]    = V;       wh[a,b,c]    = W
            uh[ma,mb,mc] = conj(U); vh[ma,mb,mc] = conj(V); wh[ma,mb,mc] = conj(W)
        end
    end

    # In-place inverse FFTs (reuse one plan). With Hermitian symmetry the result
    # is real up to round-off, so taking real() loses no physical signal.
    P = plan_ifft!(uh)
    P * uh; P * vh; P * wh

    U = Array{Float64,3}(undef, N, N, N)
    V = Array{Float64,3}(undef, N, N, N)
    W = Array{Float64,3}(undef, N, N, N)
    s2 = 0.0
    @inbounds for i in eachindex(uh)
        ru = real(uh[i]); rv = real(vh[i]); rw = real(wh[i])
        U[i] = ru; V[i] = rv; W[i] = rw
        s2 += ru*ru + rv*rv + rw*rw
    end
    s = sqrt(N^3 / s2)                           # rescale to ⟨|u|²⟩ = 1
    @inbounds for i in eachindex(U)
        U[i] *= s; V[i] *= s; W[i] *= s
    end
    return U, V, W
end

# Periodic trilinear interpolation of three uniform N³ fields on [-π,π]³ at (x,y,z),
# sharing one weight/index computation. NOTE: trilinear interpolation is a low-pass
# filter; for spectrally resolved DNS prefer a generation grid N at least as fine as
# the SEM node spacing so high-k content is not damped before it is sampled.
@inline function _trilin3(U, V, W, N, x, y, z)
    h = 2π/N
    fx = (x + π)/h
    fy = (y + π)/h
    fz = (z + π)/h
    i0 = floor(Int, fx)
    j0 = floor(Int, fy)
    k0 = floor(Int, fz)
    tx = fx - i0
    ty = fy - j0
    tz = fz - k0
    wrap(i) = mod(i, N) + 1
    i1 = wrap(i0); i2 = wrap(i0+1)
    j1 = wrap(j0); j2 = wrap(j0+1)
    k1 = wrap(k0); k2 = wrap(k0+1)
    w000 = (1-tx)*(1-ty)*(1-tz); w100 = tx*(1-ty)*(1-tz)
    w010 = (1-tx)*ty*(1-tz);     w110 = tx*ty*(1-tz)
    w001 = (1-tx)*(1-ty)*tz;     w101 = tx*(1-ty)*tz
    w011 = (1-tx)*ty*tz;         w111 = tx*ty*tz
    @inline interp(F) = @inbounds (
        F[i1,j1,k1]*w000 + F[i2,j1,k1]*w100 + F[i1,j2,k1]*w010 + F[i2,j2,k1]*w110 +
        F[i1,j1,k2]*w001 + F[i2,j1,k2]*w101 + F[i1,j2,k2]*w011 + F[i2,j2,k2]*w111)
    return interp(U), interp(V), interp(W)
end

function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    comm = MPI.COMM_WORLD; rank = MPI.Comm_rank(comm)
    rank == 0 && println(" Initialize DHIT (Rogallo synthetic turbulence) ........................ ")

    qvars    = ["ρ", "ρu", "ρv", "ρw", "ρe"]
    qoutvars = ["ρ", "u", "v", "w", "p"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    PhysConst = PhysicalConst{Float64}(); γ = PhysConst.γ; γm1 = PhysConst.γm1
    Ms = 0.1; p0 = (1.0/Ms)^2 * 1.0/γ; ρe0 = p0/γm1     # c=√(γp0)=10 → Ma_rms≈0.06

    # Same seed on every rank ⇒ identical field; each rank interpolates its own nodes.
    Ngen = Int(get(inputs, :rogallo_N, 64))
    kp   = Float64(get(inputs, :rogallo_kp, 2.0))
    seed = Int(get(inputs, :rogallo_seed, 1))
    U, V, W = _rogallo_field(Ngen, kp, seed)

    for ip = 1:mesh.npoin
        x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]
        u, v, w = _trilin3(U, V, W, Ngen, x, y, z)
        ρ  = 1.0
        ρe = p0/γm1 + 0.5*ρ*(u*u + v*v + w*w)
        q.qn[ip,1] = ρ;  q.qn[ip,2] = ρ*u; q.qn[ip,3] = ρ*v; q.qn[ip,4] = ρ*w; q.qn[ip,5] = ρe; q.qn[ip,end] = p0
        q.qe[ip,1] = 1.0; q.qe[ip,2] = 0.0; q.qe[ip,3] = 0.0; q.qe[ip,4] = 0.0; q.qe[ip,5] = ρe0; q.qe[ip,end] = p0
    end

    rank == 0 && println(" Initialize DHIT (Rogallo synthetic turbulence) ........................ DONE ")
    return q
end
