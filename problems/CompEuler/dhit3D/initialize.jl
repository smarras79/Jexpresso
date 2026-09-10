using FFTW: ifft
using Random: MersenneTwister

# Rogallo synthetic isotropic turbulence on a uniform N³ periodic grid over [-π,π]³.
# Solenoidal velocity with the [41] Yamazaki-Ishihara-Kaneda initial spectrum shape
# E(k) ∝ k⁴ exp(-2(k/kp)²); the field is rescaled to ⟨|u|²⟩ = 1 (3u0²=1, urms=u0=0.577/comp).
function _rogallo_field(N::Int, kp::Float64, seed::Int)
    rng = MersenneTwister(seed)
    uh = zeros(ComplexF64, N, N, N)
    vh = similar(uh)
    wh = similar(uh)
    kf = [ i <= N÷2 ? i : i-N for i in 0:N-1 ]
    shapeE(k) = k == 0 ? 0.0 : k^4 * exp(-2.0*(k/kp)^2)
    @inbounds for a in 1:N, b in 1:N, c in 1:N
        k1 = kf[a]
	k2 = kf[b]
	k3 = kf[c]
        k = sqrt(k1^2 + k2^2 + k3^2); k == 0 && continue
        amp = sqrt(shapeE(k)/(4π*k^2))
        θ1 = 2π*rand(rng); θ2 = 2π*rand(rng); φ = 2π*rand(rng)
        α = amp*cis(θ1)*cos(φ)
	β = amp*cis(θ2)*sin(φ)
	kperp = sqrt(k1^2 + k2^2)
        if kperp < 1e-12
            uh[a,b,c] = α; vh[a,b,c] = β
        else
            uh[a,b,c] = (α*k*k2 + β*k1*k3)/(k*kperp)
            vh[a,b,c] = (β*k2*k3 - α*k*k1)/(k*kperp)
            wh[a,b,c] = -β*kperp/k
        end
    end
    u = real(ifft(uh))*N^3; v = real(ifft(vh))*N^3; w = real(ifft(wh))*N^3
    s = sqrt(1.0 / (sum(u.^2 .+ v.^2 .+ w.^2)/N^3))   # rescale to ⟨|u|²⟩ = 1
    return u.*s, v.*s, w.*s
end

# Periodic trilinear interpolation of a uniform N³ field on [-π,π]³ at (x,y,z).
@inline function _trilin(F, N, x, y, z)
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
    i1 = wrap(i0)
    i2 = wrap(i0+1)
    j1 = wrap(j0)
    j2 = wrap(j0+1)
    k1 = wrap(k0)
    k2 = wrap(k0+1)
    @inbounds return (
        F[i1,j1,k1]*(1-tx)*(1-ty)*(1-tz) + F[i2,j1,k1]*tx*(1-ty)*(1-tz) +
        F[i1,j2,k1]*(1-tx)*ty*(1-tz)     + F[i2,j2,k1]*tx*ty*(1-tz)     +
        F[i1,j1,k2]*(1-tx)*(1-ty)*tz     + F[i2,j1,k2]*tx*(1-ty)*tz     +
        F[i1,j2,k2]*(1-tx)*ty*tz         + F[i2,j2,k2]*tx*ty*tz)
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
        u = _trilin(U, Ngen, x, y, z); v = _trilin(V, Ngen, x, y, z); w = _trilin(W, Ngen, x, y, z)
        ρ  = 1.0
        ρe = p0/γm1 + 0.5*ρ*(u*u + v*v + w*w)
        q.qn[ip,1] = ρ;  q.qn[ip,2] = ρ*u; q.qn[ip,3] = ρ*v; q.qn[ip,4] = ρ*w; q.qn[ip,5] = ρe; q.qn[ip,end] = p0
        q.qe[ip,1] = 1.0; q.qe[ip,2] = 0.0; q.qe[ip,3] = 0.0; q.qe[ip,4] = 0.0; q.qe[ip,5] = ρe0; q.qe[ip,end] = p0
    end

    rank == 0 && println(" Initialize DHIT (Rogallo synthetic turbulence) ........................ DONE ")
    return q
end
