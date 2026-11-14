# Shallow-water limiter: depth floor and the :swe_umax velocity cap.
# Limiter numerics only; no mesh, MPI initialization, or GPU is needed.
# Run: julia --startup-file=no test/test_shallow_water_velocity_cap.jl
using Test

module SWLimiterTest
    using Printf
    # Minimal stand-in for the KernelAbstractions backend tag the limiter
    # dispatches on; the device branch is not exercised here.
    struct CPU end
    include(joinpath(@__DIR__, "..", "src", "kernel", "physics", "shallow_water_limiter.jl"))
end
const SWL = SWLimiterTest

function limiter_params(npoin; hmin=1.0e-3, umax=0.0, diagnostic=false)
    inputs = Dict{Symbol,Any}(:swe_hmin => hmin, :swe_umax => umax,
                              :swe_limiter_diagnostic => diagnostic)
    return (mesh=(npoin=npoin,), neqs=3, inputs=inputs, backend=SWL.CPU(),
            swe_limiter_diag=SWL.ShallowWaterLimiterDiagnostics())
end

@testset "shallow-water depth floor and velocity cap" begin
    #          dry      floor-ok  thin+fast  thin+slow  deep+fast   deep+slow
    H  = [-2.0e-3,   1.0e-3,   2.5e-3,    5.0e-3,    150.0,      150.0]
    Hu = [  0.3,     1.0e-3,   0.24,      1.0e-3,    -400.0,     30.0]
    Hv = [ -0.1,     0.0,     -0.32,      2.0e-3,     300.0,    -10.0]
    n  = length(H)
    hmin = 1.0e-3; umax = 2.0

    # cap off: only the depth floor acts (the behaviour before :swe_umax)
    u = vcat(H, Hu, Hv)
    SWL.apply_shallow_water_limiter!(u, limiter_params(n; hmin=hmin, umax=0.0), 0.0)
    @test u[1] == hmin && u[n+1] == 0.0 && u[2n+1] == 0.0
    @test u[2:n] == H[2:n]
    @test u[n+2:2n] == Hu[2:n] && u[2n+2:3n] == Hv[2:n]

    # cap on
    u = vcat(H, Hu, Hv)
    SWL.apply_shallow_water_limiter!(u, limiter_params(n; hmin=hmin, umax=umax), 0.0)
    Hn = u[1:n]; Hun = u[n+1:2n]; Hvn = u[2n+1:3n]
    @test Hn == max.(H, hmin)                       # the cap never touches H
    @test Hun[1] == 0.0 && Hvn[1] == 0.0            # floored node: momentum zeroed
    spd = hypot.(Hun, Hvn) ./ Hn
    @test all(spd .<= umax * (1 + 4eps()))
    for i in (3, 5)                                 # capped exactly to umax, direction kept
        @test spd[i] ≈ umax
        @test atan(Hvn[i], Hun[i]) ≈ atan(Hv[i], Hu[i])
    end
    for i in (2, 4, 6)                              # admissible nodes untouched, bit for bit
        @test Hun[i] == Hu[i] && Hvn[i] == Hv[i]
    end

    # diagnostics count the capped node-visits and the largest |u| removed
    p = limiter_params(n; hmin=hmin, umax=umax, diagnostic=true)
    Hc = max.(H, hmin); Huc = copy(Hu); Hvc = copy(Hv)
    Huc[1] = 0.0; Hvc[1] = 0.0
    SWL._cap_shallow_water_velocity!(p, Hc, Huc, Hvc, umax)
    @test p.swe_limiter_diag.capped_nodes == 2
    @test p.swe_limiter_diag.max_u_before ≈ hypot(0.24, 0.32) / 2.5e-3
    SWL.reset_shallow_water_limiter_diagnostics!(p.swe_limiter_diag)
    @test p.swe_limiter_diag.capped_nodes == 0 && p.swe_limiter_diag.max_u_before == 0.0

    # Float32 state
    u32 = Float32.(vcat(H, Hu, Hv))
    SWL.apply_shallow_water_limiter!(u32, limiter_params(n; hmin=hmin, umax=umax), 0.0)
    @test eltype(u32) == Float32
    @test all(hypot.(u32[n+1:2n], u32[2n+1:3n]) ./ u32[1:n] .<= Float32(umax) * (1 + 4eps(Float32)))
end
