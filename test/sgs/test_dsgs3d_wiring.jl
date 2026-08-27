#=============================================================================
 test_dsgs3d_wiring.jl -- the REAL SGS_DSGS code, loaded and run.

 test_dsgs3d.jl transcribes the 3D DynSGS kernel and checks the ALGORITHM.
 This file does the complementary half: it includes the actual
 sgsStructs.jl and SGS.jl and exercises the functions as written, so a
 signature that does not match its call site, a field that is not on the
 struct, or an allocator that hands 2D a closure it must not have, fails here
 rather than three hours into a 15.9 M-node run.

 It loads those two files against stubs for the handful of things they reach
 outside themselves -- KernelAbstractions.zeros, MPI (serial: the reductions
 are the identity), the saturation functions, @turbo -- instead of loading
 Jexpresso, for the same reason test_closures.jl does not: the point is to
 test the closure, not to pay for Gridap and MPI to precompile first.

 Run:  julia --project=. test/sgs/test_dsgs3d_wiring.jl
=============================================================================#
module Mini

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# --- stubs ------------------------------------------------------------------
module KernelAbstractions
    zeros(backend, T, dims) = Base.zeros(T, dims)
end
module MPI
    const SUM = :sum
    const MAX = :max
    Allreduce!(x, op, comm) = x            # serial: identity
    Allreduce(x, op, comm)  = x
end
get_mpi_comm() = nothing
qsatw(T,p)   = 0.01
qsati(T,p)   = 0.01
dtqsatw(T,p) = 0.0
dtqsati(T,p) = 0.0
macro turbo(ex); esc(ex); end

using Parameters
include(joinpath(ROOT, "src", "kernel", "abstractTypes.jl"))
include(joinpath(ROOT, "src", "kernel", "physics", "globalConstantsPhysics.jl"))
include(joinpath(ROOT, "src", "kernel", "physics", "sgsStructs.jl"))
include(joinpath(ROOT, "src", "kernel", "physics", "SGS.jl"))
end # module

using Test
using .Mini
const M = Mini

const PC = M.PhysicalConst{Float64}()

# ---- a 2-element, 2x2x2-node synthetic mesh -------------------------------
const NGL   = 2
const NELEM = 2
const NEQS  = 5
const NPOIN = NELEM*NGL^3

conn = zeros(Int64, NELEM, NGL, NGL, NGL)
let ip = 0
    for ie in 1:NELEM, k in 1:NGL, j in 1:NGL, i in 1:NGL
        ip += 1; conn[ie,i,j,k] = ip
    end
end

sgs = M.allocate_SGS(NPOIN, Float64, nothing, PC, M.DSGS();
                     C_s = 0.16, nelem = NELEM, neqs = NEQS, SD = M.NSD_3D(),
                     C1 = 1.0, C2 = 0.5)

@testset "SGS_DSGS wiring" begin

@test sgs isa M.SGS_DSGS
@test sgs isa M.AbstractSGSModel
# 1D/2D must NOT get a closure struct: viscous_rhs_el! there has its own
# DSGS branch and several diagnostics read `params.sgs isa AbstractSGSModel`
# as "there is a closure to ask".
@test M.allocate_SGS(NPOIN, Float64, nothing, PC, M.DSGS(); SD = M.NSD_2D()) === nothing
@test M.allocate_SGS(NPOIN, Float64, nothing, PC, M.DSGS(); SD = M.NSD_1D()) === nothing
@test M.allocate_SGS(NPOIN, Float64, nothing, PC, M.DSGS()) === nothing
@test length(sgs.μ_el) == NELEM && length(sgs.μ_turb) == NPOIN
@test length(sgs.avg) == NEQS && length(sgs.denom) == NEQS && length(sgs.scale) == NEQS
@test sgs.Pr_t == PC.Pr_t && sgs.C_s2 ≈ 0.16^2

# ---- the real compute_dsgs_viscosity! -------------------------------------
ρ0, θ0 = 1.2, 300.0
q  = zeros(NPOIN, NEQS);  q[:,1] .= ρ0;  q[:,5] .= ρ0*θ0
qe = copy(q)
q[1,5] += ρ0*1.0                              # 1 K on one node
q1 = copy(q); q2 = copy(q)
rhs  = zeros(NPOIN, NEQS);  rhs[3,5] = -0.5   # residual in element 1 only
Minv = ones(NPOIN)
Δf   = fill(160.0, NELEM)                     # 160 m element, nop = 4 -> Δ = 40 m
visc = [0.0, 1.0, 1.0, 1.0, 1.0]
μ_dsgs = zeros(NELEM, NEQS)

M.compute_dsgs_viscosity!(sgs, M.NSD_3D(), μ_dsgs,
                          q, q1, q2, qe, rhs, Minv, visc,
                          0.5, conn, Δf, 40.0, 4, PC, nothing,
                          NELEM, NGL, NEQS; lpert = false, lglobal_norms = false)

d5 = max(ρ0*1.0 - ρ0*1.0/NPOIN, 1.0e-3*ρ0*θ0)
@test sgs.ν_el[1] ≈ 1.0*40.0^2*0.5/d5 rtol = 1e-12   # eq. (9), below the cap
@test sgs.ν_el[2] == 0.0                             # no residual -> exactly 0
@test sgs.μ_el[1] ≈ ρ0*sgs.ν_el[1] rtol = 1e-12      # dynamic = rho_bar * nu
# The VTK slots mirror SGS_diffusion's split, and :μ[1] = 0 keeps mass clean.
@test μ_dsgs[1,1] == 0.0
@test μ_dsgs[1,2] ≈ sgs.μ_el[1] rtol = 1e-14
@test μ_dsgs[1,5] ≈ sgs.μ_el[1]/PC.Pr_t rtol = 1e-14

# The global-norm branch must produce the same answer on one rank.
νser = copy(sgs.ν_el)
M.compute_dsgs_viscosity!(sgs, M.NSD_3D(), μ_dsgs,
                          q, q1, q2, qe, rhs, Minv, visc,
                          0.5, conn, Δf, 40.0, 4, PC, nothing,
                          NELEM, NGL, NEQS; lpert = false, lglobal_norms = true)
@test sgs.ν_el == νser

# ---- the real compute_sgs_cache! ------------------------------------------
# uprimitive is (ngl,ngl,ngl,neqs) = (rho,u,v,w,theta); metrics are identity
# over a 1 m element so d/dxi == d/dx.
uprim = zeros(NGL,NGL,NGL,NEQS)
for k in 1:NGL, j in 1:NGL, i in 1:NGL
    uprim[i,j,k,1] = ρ0
    uprim[i,j,k,2] = 1.7*(k-1)      # du/dz = 1.7
    uprim[i,j,k,5] = θ0
end
one4 = ones(NELEM,NGL,NGL,NGL);  zero4 = zeros(NELEM,NGL,NGL,NGL)
dψ   = [-1.0 -1.0; 1.0 1.0]       # ngl = 2 LGL on [-1,1] scaled to unit length
mp   = (Tabs = zeros(1,1), qn = zeros(1), qsatt = zeros(1))

M.compute_sgs_cache!(sgs, uprim, mp, q, NGL, dψ,
                     one4, zero4, zero4,
                     zero4, one4, zero4,
                     zero4, zero4, one4,
                     conn, 1, 40.0^2, 1, M.NSD_3D())

# μ_turb over element 1 is the element's DynSGS coefficient, uniformly.
@test all(sgs.μ_turb[conn[1,i,j,k]] ≈ sgs.μ_el[1] for i in 1:NGL, j in 1:NGL, k in 1:NGL)
# ... and Sij is filled for les_statistics: pure du/dz shear -> S13 = s/2.
@test sgs.S13[conn[1,1,1,1]] ≈ 0.5*1.7 rtol = 1e-12
@test sgs.S11[conn[1,1,1,1]] ≈ 0.0 atol = 1e-14

# With :dsgs_add_smagorinsky the two viscosities ADD, and the Smagorinsky part
# is the textbook rho (Cs Delta)^2 |S| with |S| = |du/dz| for pure shear.
sgs.ladd_smagorinsky = true
M.compute_sgs_cache!(sgs, uprim, mp, q, NGL, dψ,
                     one4, zero4, zero4,
                     zero4, one4, zero4,
                     zero4, zero4, one4,
                     conn, 1, 40.0^2, 1, M.NSD_3D())
@test sgs.μ_turb[conn[1,1,1,1]] ≈ sgs.μ_el[1] + ρ0*0.16^2*40.0^2*1.7 rtol = 1e-12
sgs.ladd_smagorinsky = false

# ---- SGS_diffusion reads it exactly as it reads Smagorinsky's -------------
μt = sgs.μ_turb[conn[1,1,1,1]]
@test M.SGS_diffusion(visc, 2, ρ0, conn[1,1,1,1], sgs, true, M.NSD_3D()) ≈ PC.μ_mol + μt
@test M.SGS_diffusion(visc, 5, ρ0, conn[1,1,1,1], sgs, true, M.NSD_3D()) ≈ μt/PC.Pr_t
@test M.SGS_diffusion(visc, 1, ρ0, conn[1,1,1,1], sgs, true, M.NSD_3D()) == 0.0

end
