#=============================================================================
 test_vdiffusion.jl -- the implicit VERTICAL DIFFUSION operator (vdiffusion.jl).

     julia --project=<env> test/hevi/test_vdiffusion.jl
     mpiexecjl -n 3 julia --project=<env> test/hevi/test_vdiffusion.jl

 WHAT HAS TO BE TRUE, AND WHY EACH CHECK IS HERE
 -----------------------------------------------
 The implicit half of an ARK split is computed as f_exp = rhs! - f_imp, so a
 WRONG diffusion operator cannot make the scheme inconsistent -- it can only
 fail to remove the stiffness it was added to remove, and that failure looks
 exactly like the blow-up it was supposed to prevent. Nothing downstream would
 point at this file. So the checks below are the only thing standing between a
 sign error here and a week spent looking somewhere else.

   1. it IS d/dz(mu d/dz), to the discretisation's own order
   2. it is self-adjoint and negative semidefinite, which a transposed dpsi
      index or a misplaced metric factor would destroy while leaving (1)
      looking plausible
   3. it agrees TERM FOR TERM with the explicit viscous discretisation --
      the reference implementations below are transcribed verbatim from
      _expansion_visc! in rhs.jl -- because cancellation between a weak form
      and anything else is not cancellation
   4. the 4/3 on the w equation, the one constant here derived on paper
   5. the probed band picks the diffusion up, so HEVI's column LU and the 3D
      IMEX's preconditioner both invert the operator they are actually given
=============================================================================#

using Test, MPI, LinearAlgebra, Printf, OrdinaryDiffEq

MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "mock_sem.jl"))
for f in ("columns.jl", "vdiffusion.jl", "operator.jl", "factorize.jl",
          "acoustic.jl", "ark.jl", "hevi.jl")
    include(joinpath(SRC, f))
end

say(args...) = RANK == 0 && (println(args...); flush(stdout))

const LX, LY, LZ = 400.0, 400.0, 1200.0
const AVMU = [0.0, 3.0, 3.0, 3.0, 5.0]      # AV: these ARE the diffusivities
const SGSMU = [0.0, 1.0, 1.0, 1.0, 2.0]     # SGS: a mask on the closure's mu

say("\n=== implicit vertical diffusion: $NR rank(s) ===")

#-----------------------------------------------------------------------------
# Reference implementations, transcribed from rhs.jl. Kept verbatim -- all
# three sweeps, the full metric contraction, the same dpsi index order -- so
# that "it matches" means it matches what the explicit code actually does and
# not what this file thinks the explicit code does.
#-----------------------------------------------------------------------------

# rhs.jl:2523 -- _expansion_visc!(..., VT::AV, SD::NSD_3D, ::ContGal)
function ref_visc_av!(rξ, rη, rζ, up, vc, ω, ngl, dψ, Je,
                      dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, iel, ieq)
    for m = 1:ngl, l = 1:ngl, k = 1:ngl
        @inbounds begin
            ωJac = ω[k] * ω[l] * ω[m] * Je[iel,k,l,m]
            dqdξ = 0.0; dqdη = 0.0; dqdζ = 0.0
            for ii = 1:ngl
                dqdξ += dψ[ii,k]*up[ii,l,m,ieq]
                dqdη += dψ[ii,l]*up[k,ii,m,ieq]
                dqdζ += dψ[ii,m]*up[k,l,ii,ieq]
            end
            ax=dξdx[iel,k,l,m]; ay=dξdy[iel,k,l,m]; az=dξdz[iel,k,l,m]
            bx=dηdx[iel,k,l,m]; by=dηdy[iel,k,l,m]; bz=dηdz[iel,k,l,m]
            cx=dζdx[iel,k,l,m]; cy=dζdy[iel,k,l,m]; cz=dζdz[iel,k,l,m]
            dqdx = vc[ieq]*(dqdξ*ax + dqdη*bx + dqdζ*cx)
            dqdy = vc[ieq]*(dqdξ*ay + dqdη*by + dqdζ*cy)
            dqdz = vc[ieq]*(dqdξ*az + dqdη*bz + dqdζ*cz)
            ∇ξ = (ax*dqdx + ay*dqdy + az*dqdz)*ωJac
            ∇η = (bx*dqdx + by*dqdy + bz*dqdz)*ωJac
            ∇ζ = (cx*dqdx + cy*dqdy + cz*dqdz)*ωJac
            for i = 1:ngl
                rξ[iel,i,l,m,ieq] -= dψ[i,k]*∇ξ
                rη[iel,k,i,m,ieq] -= dψ[i,l]*∇η
                rζ[iel,k,l,i,ieq] -= dψ[i,m]*∇ζ
            end
        end
    end
end

# The stub closure the SGS reference reads through. `SGS_diffusion` below is
# rhs-side SGS.jl:1413 verbatim; a z-varying mu_turb means a coefficient
# silently treated as constant anywhere would not pass.
struct StubSGS <: AbstractSGSModel
    μ_turb::Vector{Float64}; ltheta_eqn::Bool
    Pr_t::Float64; Sc_t::Float64; μ_mol::Float64; κ_mol::Float64
end
function SGS_diffusion(vc, ieq, ρ, ip, sgs::StubSGS, ltheta_eqn, ::NSD_3D)
    μt = sgs.μ_turb[ip]
    if ieq in (2, 3, 4)
        return (sgs.μ_mol + μt) * vc[ieq]
    elseif ieq == 5
        κ = μt / sgs.Pr_t
        return ltheta_eqn ? κ * vc[ieq] : (ρ*sgs.κ_mol + κ) * vc[ieq]
    else
        return (ρ*sgs.κ_mol + μt/sgs.Sc_t) * vc[ieq]
    end
end

# rhs.jl:2619 -- _expansion_visc!(..., sgs::AbstractSGSModel, ..., NSD_3D)
function ref_visc_sgs!(rξ, rη, rζ, up, vc, ω, ngl, dψ, Je,
                       dξdx,dξdy,dξdz, dηdx,dηdy,dηdz, dζdx,dζdy,dζdz,
                       iel, ieq, conn, sgs, ltheta)
    for m=1:ngl, l=1:ngl, k=1:ngl
        ip = conn[iel,k,l,m]
        @inbounds begin
        ωJac = ω[k]*ω[l]*ω[m]*Je[iel,k,l,m]
        dudξ=0.0;dudη=0.0;dudζ=0.0; dvdξ=0.0;dvdη=0.0;dvdζ=0.0; dwdξ=0.0;dwdη=0.0;dwdζ=0.0
        for ii=1:ngl
            dudξ += dψ[ii,k]*up[ii,l,m,2]; dudη += dψ[ii,l]*up[k,ii,m,2]; dudζ += dψ[ii,m]*up[k,l,ii,2]
            dvdξ += dψ[ii,k]*up[ii,l,m,3]; dvdη += dψ[ii,l]*up[k,ii,m,3]; dvdζ += dψ[ii,m]*up[k,l,ii,3]
            dwdξ += dψ[ii,k]*up[ii,l,m,4]; dwdη += dψ[ii,l]*up[k,ii,m,4]; dwdζ += dψ[ii,m]*up[k,l,ii,4]
        end
        ax=dξdx[iel,k,l,m]; ay=dξdy[iel,k,l,m]; az=dξdz[iel,k,l,m]
        bx=dηdx[iel,k,l,m]; by=dηdy[iel,k,l,m]; bz=dηdz[iel,k,l,m]
        cx=dζdx[iel,k,l,m]; cy=dζdy[iel,k,l,m]; cz=dζdz[iel,k,l,m]
        dudx=dudξ*ax+dudη*bx+dudζ*cx; dudy=dudξ*ay+dudη*by+dudζ*cy; dudz=dudξ*az+dudη*bz+dudζ*cz
        dvdx=dvdξ*ax+dvdη*bx+dvdζ*cx; dvdy=dvdξ*ay+dvdη*by+dvdζ*cy; dvdz=dvdξ*az+dvdη*bz+dvdζ*cz
        dwdx=dwdξ*ax+dwdη*bx+dwdζ*cx; dwdy=dwdξ*ay+dwdη*by+dwdζ*cy; dwdz=dwdξ*az+dwdη*bz+dwdζ*cz
        div_u = dudx + dvdy + dwdz
        ρ = up[k,l,m,1]
        μ = SGS_diffusion(vc, ieq, ρ, ip, sgs, ltheta, NSD_3D())
        local fx, fy, fz
        if ieq == 2
            fx = 2.0*μ*dudx - (2.0/3.0)*μ*div_u; fy = μ*(dudy+dvdx); fz = μ*(dudz+dwdx)
        elseif ieq == 3
            fx = μ*(dudy+dvdx); fy = 2.0*μ*dvdy - (2.0/3.0)*μ*div_u; fz = μ*(dvdz+dwdy)
        elseif ieq == 4
            fx = μ*(dudz+dwdx); fy = μ*(dvdz+dwdy); fz = 2.0*μ*dwdz - (2.0/3.0)*μ*div_u
        else
            dsξ=0.0; dsη=0.0; dsζ=0.0
            for ii=1:ngl
                dsξ += dψ[ii,k]*up[ii,l,m,ieq]; dsη += dψ[ii,l]*up[k,ii,m,ieq]; dsζ += dψ[ii,m]*up[k,l,ii,ieq]
            end
            fx = μ*(dsξ*ax+dsη*bx+dsζ*cx); fy = μ*(dsξ*ay+dsη*by+dsζ*cy); fz = μ*(dsξ*az+dsη*bz+dsζ*cz)
        end
        ∇ξ = (ax*fx+ay*fy+az*fz)*ωJac
        ∇η = (bx*fx+by*fy+bz*fz)*ωJac
        ∇ζ = (cx*fx+cy*fy+cz*fz)*ωJac
        for i=1:ngl
            rξ[iel,i,l,m,ieq] -= dψ[i,k]*∇ξ
            rη[iel,k,i,m,ieq] -= dψ[i,l]*∇η
            rζ[iel,k,l,i,ieq] -= dψ[i,m]*∇ζ
        end
        end
    end
end

# The explicit viscous RHS of a state, assembled exactly as rhs.jl assembles
# it: primitives per element, the loop above, DSS, MPI, then Minv.
function explicit_visc_rhs(pa, u; sgs = nothing)
    mesh = pa.mesh; met = pa.metrics
    nelem = Int(mesh.nelem); ngl = Int(mesh.ngl); npoin = Int(mesh.npoin); neqs = 5
    conn = mesh.connijk
    rξ = zeros(nelem,ngl,ngl,ngl,neqs); rη = zeros(size(rξ)); rζ = zeros(size(rξ))
    up = zeros(ngl,ngl,ngl,neqs)
    for iel = 1:nelem
        for m=1:ngl, l=1:ngl, k=1:ngl
            ip = conn[iel,k,l,m]
            up[k,l,m,1] = u[ip,1]
            for ieq = 2:neqs; up[k,l,m,ieq] = u[ip,ieq]/u[ip,1]; end
        end
        for ieq = 1:neqs
            if sgs === nothing
                ref_visc_av!(rξ, rη, rζ, up, pa.visc_coeff, pa.ω, ngl, pa.basis.dψ, met.Je,
                             met.dξdx, met.dξdy, met.dξdz, met.dηdx, met.dηdy, met.dηdz,
                             met.dζdx, met.dζdy, met.dζdz, iel, ieq)
            else
                ref_visc_sgs!(rξ, rη, rζ, up, pa.visc_coeff, pa.ω, ngl, pa.basis.dψ, met.Je,
                              met.dξdx, met.dξdy, met.dξdz, met.dηdx, met.dηdy, met.dηdz,
                              met.dζdx, met.dζdy, met.dζdz, iel, ieq, conn, sgs, sgs.ltheta_eqn)
            end
        end
    end
    rd = rξ .+ rη .+ rζ
    RHS = zeros(npoin, neqs)
    DSS_rhs!(RHS, rd, conn, nelem, ngl, neqs, pa.SD, pa.AD)
    ca = setup_assembler(mesh.SD, RHS, mesh.ip2gip, mesh.gip2owner)
    ca === nothing || assemble_mpi!(RHS, ca)
    vaux = zeros(npoin)
    for ieq = 1:neqs
        divide_by_mass_matrix!(@view(RHS[:,ieq]), vaux, pa.Minv, neqs, npoin, pa.AD)
    end
    return RHS
end

# The diffusion part of the operator, isolated by differencing against the
# same operator built without it. Nothing else about the two differs.
function diffusion_part(pa, topo, V; mu, sgs = nothing, vars = [1,2,3,4,5], scale43 = 1.0)
    pa2 = sgs === nothing ? pa : merge(pa, (sgs = sgs,))
    o1 = build_hevi_operator(pa2, topo, vars; vdiff = true)
    o0 = build_hevi_operator(pa2, topo, vars; vdiff = false)
    scale43 == 1.0 || (o1.vd.mu[:, findfirst(==(4), vars)] .*= scale43)
    n = Int(pa.mesh.npoin)
    R1 = zeros(n, length(vars)); R0 = zeros(n, length(vars))
    hevi_apply_A!(R1, V, pa2, o1); hevi_apply_A!(R0, V, pa2, o0)
    return R1 .- R0, o1
end

@testset "implicit vertical diffusion" begin

#--- 1. it is d/dz(mu d/dz), and converges like one --------------------------
# Feeding rho_bar * sin(kz) makes the PRIMITIVE exactly sin(kz), so the answer
# is -mu k^2 sin(kz) with no approximation left in it but the discretisation's
# own. Nodes within a quarter-wavelength of the floor and the lid are excluded:
# the weak form's natural zero-flux condition is correct there and is not
# -mu k^2 sin(kz).
function smooth_err(nelz, p)
    pa = build_mock_params(nelx=2, nely=2, nelz=nelz, p=p, Lx=LX, Ly=LY, Lz=LZ,
                           comm=COMM, mu=AVMU)
    np = Int(pa.mesh.npoin); tp = build_column_topology(pa.mesh, COMM)
    k = 4pi / LZ
    V = zeros(np, 5)
    for ip = 1:np; V[ip,5] = pa.qp.qe[ip,1] * sin(k*pa.mesh.coords[3,ip]); end
    D, _ = diffusion_part(pa, tp, V; mu = AVMU)
    e = 0.0; r = 0.0
    for ip = 1:np
        z = pa.mesh.coords[3,ip]
        (z < 250.0 || z > LZ - 250.0) && continue
        ex = -AVMU[5] * k^2 * sin(k*z)
        e = max(e, abs(D[ip,5] - ex)); r = max(r, abs(ex))
    end
    NR > 1 && (e = MPI.Allreduce(e, MPI.MAX, COMM); r = MPI.Allreduce(r, MPI.MAX, COMM))
    return e / r
end
e6  = smooth_err(6, 3); e12 = smooth_err(12, 3); e24 = smooth_err(24, 3)
say(@sprintf("  D vs -mu k^2 sin(kz), p=3:  nelz 6 %.2e -> 12 %.2e -> 24 %.2e", e6, e12, e24))
@test e6 < 5.0e-2
@test e12 < e6 / 4          # h-convergence, not a coincidence at one resolution
@test e24 < e12 / 4
@test smooth_err(12, 5) < 1.0e-4    # and spectral in p

#--- 2. self-adjoint, negative semidefinite ----------------------------------
# The continuous operator d/dz(mu d/dz) is symmetric and dissipative; the weak
# form inherits both exactly. A transposed dpsi index gives a skew operator
# that would still look like a second derivative on a smooth mode.
pa = build_mock_params(nelx=2, nely=2, nelz=6, p=3, Lx=LX, Ly=LY, Lz=LZ,
                       comm=COMM, mu=AVMU)
np = Int(pa.mesh.npoin); topo = build_column_topology(pa.mesh, COMM)
gip = pa.mesh.ip2gip
# Deterministic and GLOBALLY consistent: a node shared by two ranks must get
# the same value on both, so the field is a function of the global index.
X = [sin(0.7*gip[ip] + 1.3*q) for ip = 1:np, q = 1:5]
Y = [cos(0.4*gip[ip] - 0.9*q) for ip = 1:np, q = 1:5]
DX, o1 = diffusion_part(pa, topo, X; mu = AVMU)
DY, _  = diffusion_part(pa, topo, Y; mu = AVMU)
# The multiplicity weight makes the sum a true global inner product rather
# than one that counts shared nodes once per owning rank.
wt = ones(np, 1); ca = setup_assembler(pa.mesh.SD, wt, pa.mesh.ip2gip, pa.mesh.gip2owner)
mult = ones(np); if ca !== nothing; assemble_mpi!(wt, ca); mult .= wt[:,1]; end
function minner(A, B)
    s = 0.0
    for q = 1:size(A,2), ip = 1:size(A,1)
        s += A[ip,q] * B[ip,q] / (pa.Minv[ip] * mult[ip])
    end
    return NR > 1 ? MPI.Allreduce(s, MPI.SUM, COMM) : s
end
Xs = X .* o1.vd.sc; Ys = Y .* o1.vd.sc
a = minner(DX, Ys); b = minner(DY, Xs); d = minner(DX, Xs)
say(@sprintf("  <Dx,y> = %.8g   <Dy,x> = %.8g   (rel %.2e);   <Dx,x> = %.4g",
             a, b, abs(a-b)/max(abs(a),1e-30), d))
@test abs(a - b) <= 1.0e-9 * max(abs(a), abs(b))
@test d <= 0.0

#--- 3. it agrees with the explicit discretisation, term for term ------------
# On a horizontally uniform field the explicit viscous flux divergence IS the
# column-local one for a scalar, so the theta equation must agree everywhere.
# See the header of vdiffusion.jl for why the momentum equations agree only
# in the interior: their stress tensor has horizontal flux components that do
# not vanish, and what those leave behind is a lateral-face term DSS cancels
# between elements and cannot cancel where the domain ends.
onface = [abs(pa.mesh.coords[1,ip]) < 1e-9 || abs(pa.mesh.coords[1,ip]-LX) < 1e-9 ||
          abs(pa.mesh.coords[2,ip]) < 1e-9 || abs(pa.mesh.coords[2,ip]-LY) < 1e-9
          for ip = 1:np]
function hz_uniform_states(pa)
    n = Int(pa.mesh.npoin); qe = pa.qp.qe
    u0 = zeros(n,5); u1 = zeros(n,5)
    for ip = 1:n
        z = pa.mesh.coords[3,ip]
        for ieq = 1:5; u0[ip,ieq] = qe[ip,ieq]; u1[ip,ieq] = qe[ip,ieq]; end
        for ieq in (2,3,4,5)
            u1[ip,ieq] += qe[ip,1]*(0.3*sin(2pi*z/LZ*2) + 0.1*cos(2pi*z/LZ*5) + 0.05*ieq)
        end
    end
    return u0, u1
end
u0, u1 = hz_uniform_states(pa)
V = u1 .- pa.qp.qe

# (a) the AV path: a plain scalar diffusion per equation, so EVERY equation
#     agrees everywhere -- there is no stress tensor to leave a face term.
Aexp = explicit_visc_rhs(pa, u1) .- explicit_visc_rhs(pa, u0)
Dav, _ = diffusion_part(pa, topo, V; mu = AVMU)
worst_av = maximum(abs, Dav .- Aexp) / maximum(abs, Aexp)
say(@sprintf("  AV path vs explicit, all equations, all nodes:      %.2e", worst_av))
@test worst_av < 1.0e-11

# (b) the SGS stress-tensor path
sgs = StubSGS([8.0 + 4.0*sin(2pi*pa.mesh.coords[3,ip]/LZ) for ip=1:np],
              true, 0.7, 0.7, 1.8e-5, 1.4e-5)
pas = merge(pa, (sgs = sgs, visc_coeff = SGSMU))
Sexp = explicit_visc_rhs(pas, u1; sgs=sgs) .- explicit_visc_rhs(pas, u0; sgs=sgs)
Dsgs, _ = diffusion_part(pas, topo, V; mu = SGSMU, sgs = sgs)
for ieq = 2:5
    ref = maximum(abs, @view Sexp[:,ieq])
    din = maximum([abs(Dsgs[ip,ieq]-Sexp[ip,ieq]) for ip=1:np if !onface[ip]]; init=0.0)
    say(@sprintf("  SGS path eq %d vs explicit, interior nodes:          %.2e", ieq, din/ref))
    @test din < 1.0e-11 * ref
end
# theta has no stress tensor, so it must agree on the faces too
@test maximum(abs, Dsgs[:,5] .- Sexp[:,5]) < 1.0e-11 * maximum(abs, @view Sexp[:,5])

#--- 4. the 4/3 on the w equation -------------------------------------------
# tau_zz = 2 mu dw/dz - (2/3) mu div(u), and dw/dz appears in both terms, so
# the column-local coefficient is (2 - 2/3) mu. Dropping it is a 10% error on
# the w equation -- small enough to look like something else, which is exactly
# why it is pinned here.
D43, _ = diffusion_part(pas, topo, V; mu = SGSMU, sgs = sgs, scale43 = 3/4)
ref4 = maximum(abs, @view Sexp[:,4])
bad  = maximum([abs(D43[ip,4]-Sexp[ip,4]) for ip=1:np if !onface[ip]]; init=0.0)
say(@sprintf("  w equation without the 4/3:                         %.2e (must be >> 0)", bad/ref4))
@test bad > 1.0e-2 * ref4

#--- 5. the probed band is the operator, diffusion included -----------------
# This is what makes the whole thing free: assemble_column_band probes the
# operator rather than transcribing a formula, so HEVI's column LU and the 3D
# IMEX's preconditioner absorb the diffusion with no code of their own. If
# they did not, the exact column solve would silently become an approximate
# one and the only symptom would be an iteration count.
owner, own = assign_column_owners(topo, COMM)
cc = build_column_comm(topo, owner, own, COMM, 5)
gdt = 0.02
for vd in (false, true)
    o = build_hevi_operator(pa, topo, [1,2,3,4,5]; vdiff = vd)
    AB, kl, ku, n = assemble_column_band(pa, o, cc, topo, gdt)
    berr = hevi_verify_band(pa, o, cc, topo, AB, kl, ku, n, gdt)
    fac  = build_column_factorization(pa, o, cc, topo, gdt)
    serr = hevi_verify_solve(pa, o, cc, fac)
    say(@sprintf("  vdiff=%-5s  band vs operator %.2e   banded solve residual %.2e",
                 vd, berr, serr))
    @test berr < 1.0e-10
    @test serr < 1.0e-8
end
# and the band really did change -- otherwise the two lines above would both
# pass on an operator that quietly dropped the diffusion.
ovd = build_hevi_operator(pa, topo, [1,2,3,4,5]; vdiff = true)
ono = build_hevi_operator(pa, topo, [1,2,3,4,5]; vdiff = false)
ABv, _, _, _ = assemble_column_band(pa, ovd, cc, topo, gdt)
ABn, _, _, _ = assemble_column_band(pa, ono, cc, topo, gdt)
Δband = cc.nown == 0 ? 0.0 : maximum(maximum(abs, ABv[i] .- ABn[i]) for i = 1:cc.nown)
NR > 1 && (Δband = MPI.Allreduce(Δband, MPI.MAX, COMM))
say(@sprintf("  band entries moved by the diffusion: %.3g", Δband))
@test Δband > 0.0

#--- 6. the deck-level guards ----------------------------------------------
inp_on  = Dict{Symbol,Any}(:implicit_vdiff => true,  :lvisc => true)
inp_off = Dict{Symbol,Any}(:implicit_vdiff => false, :lvisc => true)
inp_nov = Dict{Symbol,Any}(:implicit_vdiff => true,  :lvisc => false)
@test vdiff_enabled(inp_on)
@test !vdiff_enabled(inp_off)
@test !vdiff_enabled(inp_nov)               # nothing to make implicit
# u and v have to join the implicit set: their vertical diffusive rate is the
# same mu/dz^2 as theta's, so leaving them out removes none of the limit.
@test vdiff_vars(pa, inp_on,  [1,4,5]) == [1,2,3,4,5]
@test vdiff_vars(pa, inp_off, [1,4,5]) == [1,4,5]
# :RS under a dynamic closure freezes mu at its molecular floor for the whole
# run, so the operator would carry no diffusion while still costing the band.
@test vdiff_check_linearisation(inp_on, pas, "HEVI", :PS) === nothing
@test_throws ErrorException vdiff_check_linearisation(inp_on, pas, "HEVI", :RS)
@test vdiff_check_linearisation(inp_on, pa, "HEVI", :RS) === nothing   # AV: constant
# and an operator asked for diffusion over too narrow a variable set refuses
@test_throws ErrorException build_hevi_operator(pa, topo, [1,4,5]; vdiff = true)

end # testset

MPI.Barrier(COMM)
say("\n=== implicit vertical diffusion: done ===")
