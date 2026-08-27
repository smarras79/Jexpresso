#=============================================================================
 test_dsgs_energy_dispatch.jl -- the ENERGY (non-theta) DynSGS path must be
 callable, not just compile.

 WHAT THIS PINS. compute_dsgs_viscosity!(...; ltheta=false) forwards to
 _dsgs_2d_energy!, and the forwarding call omitted `wres` -- the boundary mask
 that keeps the residual from reading the boundary condition rather than a
 discretisation error. Omitting an argument in the middle of a positional list
 is not a quiet mis-binding: dispatch fails outright, so CompEuler/ffs_step
 died with

     MethodError: no method matching _dsgs_2d_energy!(::Matrix{Float64}, ...)

 on its first RHS evaluation, after the mesh, the setup and the initial VTK.

 WHY A TEST AND NOT A RUN. ffs_step's mesh lives under the gitignored meshes/
 directory, so the case cannot run in every checkout. The failure was pure
 dispatch, so a direct call reproduces it exactly and costs no mesh.

 BOTH branches are exercised -- ltheta=false (energy) and ltheta=true (theta) --
 because they take different argument lists to different kernels and only one
 of them was broken.

     julia --project=<env> test/sgs/test_dsgs_energy_dispatch.jl
=============================================================================#

using Test
using Jexpresso: compute_dsgs_viscosity!, DSGS, NSD_2D, PhysicalConst

const TT   = Float64
const NGL  = 4
const NEL  = 2
const NP   = NEL * NGL * NGL
const NEQ  = 4

mk() = (μ = zeros(TT, NEL, NEQ),
        q = ones(TT, NP, NEQ), q1 = ones(TT, NP, NEQ), q2 = ones(TT, NP, NEQ),
        qe = ones(TT, NP, NEQ), rhs = zeros(TT, NP, NEQ),
        Minv = ones(TT, NP), wres = ones(TT, NP),
        vc = ones(TT, NEQ), Δelem = fill(TT(10), NEL),
        conn = reshape(collect(1:NP), NEL, NGL, NGL, 1))

@testset "DynSGS 2D dispatch, both equation sets" begin
    PC = PhysicalConst{TT}()
    for ltheta in (false, true)
        a = mk()
        # A no-flow state: the point is that the CALL resolves and returns, not
        # what viscosity it produces.
        ok = try
            compute_dsgs_viscosity!(a.μ, DSGS(), NSD_2D(), a.q, a.q1, a.q2, a.qe,
                                    a.rhs, a.Minv, a.wres, a.vc, TT(0.1),
                                    a.conn, a.Δelem, PC, TT(0.71), NEL, NGL;
                                    ltheta = ltheta, lglobal_norms = false)
            true
        catch e
            @error "compute_dsgs_viscosity! failed for ltheta=$ltheta" exception=e
            false
        end
        @test ok
        @test all(isfinite, a.μ)          # and it produced numbers, not NaN
        @test all(>=(0), a.μ)             # viscosity is non-negative by construction
    end
end
