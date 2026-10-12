# Source-term numerics only; no mesh, MPI initialization, or GPU is needed.
# Run: julia --startup-file=no test/test_reduced_gravity_wind.jl
using Test

module ReducedGravitySourceTest
    # Minimal dispatch tags let us load the actual case flux/source files
    # without loading the full solver and its mesh dependencies.
    struct NSD_2D end
    struct St_mesh end
    struct CL end
    struct TOTAL end
    struct PERT end
    const TInt = Int64
    include(joinpath(@__DIR__, "..", "problems", "ShallowWater", "ReducedGravity", "user_flux.jl"))
    include(joinpath(@__DIR__, "..", "problems", "ShallowWater", "ReducedGravity", "user_source.jl"))
end
const RG = ReducedGravitySourceTest

@testset "ReducedGravity wind and dry-layer relaxation" begin
    for T in (Float64, Float32)
        Hdry = T(RG._H_WET_SWE)
        Hfull = T(RG._H_WIND_FULL_SWE)
        Hmid = (Hdry + Hfull) / T(2)
        @test RG._swe_wind_fraction(Hdry) === zero(T)
        @test RG._swe_wind_fraction(-Hdry) === zero(T)
        @test RG._swe_wind_fraction(Hfull) === one(T)
        @test RG._swe_wind_fraction(T(150)) === one(T)
        @test RG._swe_wind_fraction(Hmid) ≈ T(0.5)
        weights = RG._swe_wind_fraction.(range(Hdry, Hfull; length=101))
        @test issorted(weights)
        @test all(w -> zero(T) <= w <= one(T), weights)
        for (H, fraction) in ((Hdry, zero(T)), (Hmid, T(0.5)), (Hfull, one(T)), (T(150), one(T)))
            # y at the basin midpoint gives maximum positive zonal stress.
            q = T[H, 0, 0]
            qe = T[150, 0, 0]
            source = zeros(T, 3)
            RG.user_source!(source, q, qe, Int64(1), RG.CL(), RG.TOTAL();
                            y=T(1e6), ymin=zero(T), ymax=T(2e6))
            @test source[1] == zero(T)
            @test source[2] ≈ fraction * T(RG._TAU0_SWE / RG._RHO_SWE)
            @test source[3] == zero(T)
            gpu = RG.user_source_gpu(q, qe, zero(T), T(1e6), nothing,
                                     T(2e6), zero(T), T(2e6), zero(T), false)
            @test all(x -> x isa T, gpu)
            @test collect(gpu) ≈ source

            # With momentum present, check Coriolis, bottom drag and the
            # complementary dry relaxation independently of wind forcing.
            q[2:3] .= T[2, -3]
            RG.user_source!(source, q, qe, Int64(1), RG.CL(), RG.TOTAL();
                            y=T(1e6), ymin=zero(T), ymax=T(2e6))
            damping = T(RG._KAPPA_SWE) + (one(T) - fraction) * T(RG._SIGMA_DRY_SWE)
            expected = T[0, T(RG._F0_SWE)*q[3] + fraction*T(RG._TAU0_SWE/RG._RHO_SWE) - damping*q[2],
                         -T(RG._F0_SWE)*q[2] - damping*q[3]]
            @test source ≈ expected
            gpu = RG.user_source_gpu(q, qe, zero(T), T(1e6), nothing,
                                     T(2e6), zero(T), T(2e6), zero(T), false)
            @test collect(gpu) ≈ source
            pert = zeros(T, 3)
            RG.user_source!(pert, q, qe, Int64(1), RG.CL(), RG.PERT();
                            y=T(1e6), ymin=zero(T), ymax=T(2e6))
            @test pert == source
        end
    end
end
