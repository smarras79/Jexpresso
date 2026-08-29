#=============================================================================
 test_wall_damping.jl -- sgs_mixing_length2 at, and approaching, the wall.

 THE BUG THIS PINS
 -----------------
 The guard used to read

     (lwall_damping && z > 0.0) || return CsΔ2

 and params_setup.jl builds zwall as max(z - zmin, 0), so EVERY node on the
 lower boundary carries exactly 0.0 and took that early return. The near-wall
 limiter therefore ran backwards at the one place it exists for: it handed the
 wall node the full, undamped (C_sΔ)² while damping its neighbours.

 On the LESICP2 64x64x60 mesh that is ℓ² = 40.96 m² ON the wall against
 6.43 m² at the first node above it (z = 6.91 m) -- a 6.4x spike sitting on the
 node with the smallest h_z and the one the surface heat flux is injected into,
 i.e. the node that sets the explicit vertical viscous limit. It should carry
 ~0.

 So the tests below are, in order: ℓ² -> 0 at the wall; ℓ² monotone increasing
 away from it (the property the old code violated); and the untouched
 behaviours -- ℓ² -> (C_sΔ)² far away, and an exact no-op when the limiter is
 off, which is what every deck that has not asked for it relies on.

     julia --project=<env> test/sgs/test_wall_damping.jl
=============================================================================#

using Test, Printf

# Lift the real definition rather than restating it: a test carrying its own
# copy pins the copy.
const _SGS_SRC = read(joinpath(@__DIR__, "..", "..", "src", "kernel", "physics", "SGS.jl"), String)

function _lift!(needle, label)
    i = findfirst(needle, _SGS_SRC)
    i === nothing && error("$label not found in SGS.jl -- renamed?")
    j = findnext("\nend\n", _SGS_SRC, i[1])
    j === nothing && error("no closing end for $label in SGS.jl")
    include_string(Main, _SGS_SRC[i[1]:j[1]+4], "SGS.jl:$label")
end
_lift!("@inline function sgs_mixing_length2", "sgs_mixing_length2")

say(s) = println(s)

# LESICP2-64x64x60 numbers, so the regression is stated in the units the case
# actually failed in.
const C_S    = 0.16
const C_S2   = C_S * C_S
const Δ      = 40.0            # :les_filter_width => :max gives max(dx,dy,dz)/nop
const Δ2     = Δ * Δ
const KARMAN = 0.4
const CSΔ2   = C_S2 * Δ2       # 40.96 m²
const Z1     = 6.91            # first LGL node above the wall on this mesh

say("\n=== sgs_mixing_length2: the near-wall limit ===")
say(@sprintf("  (C_s*Delta)^2 = %.2f m^2 at C_s = %.2f, Delta = %.1f m", CSΔ2, C_S, Δ))

@testset "sgs_mixing_length2" begin

    @testset "the limiter is OFF: exact no-op at every z" begin
        # This is the promise to every deck that never set :lwall_damping.
        for z in (0.0, 1e-9, Z1, 100.0, 5000.0)
            @test sgs_mixing_length2(C_S2, Δ2, z, KARMAN, false) === CSΔ2
        end
    end

    @testset "ON: the wall node is damped to zero, not to (C_s*Delta)^2" begin
        # The regression. Before the fix this returned CSΔ2.
        ℓ2_wall = sgs_mixing_length2(C_S2, Δ2, 0.0, KARMAN, true)
        @test ℓ2_wall == 0.0
        @test ℓ2_wall != CSΔ2
        say(@sprintf("  z = 0.00 m -> l^2 = %.4f m^2   (was %.2f)", ℓ2_wall, CSΔ2))
    end

    @testset "ON: monotone increasing away from the wall" begin
        # The property the old guard broke: ℓ² jumped DOWN from 40.96 at z = 0
        # to 6.43 at z = 6.91 and then climbed again.
        zs  = [0.0, 0.5, 1.0, 2.0, Z1, 20.0, 50.0, 200.0, 2000.0]
        ℓ2s = [sgs_mixing_length2(C_S2, Δ2, z, KARMAN, true) for z in zs]
        @test issorted(ℓ2s)
        @test all(diff(ℓ2s) .> 0.0)
        for (z, l) in zip(zs, ℓ2s)
            say(@sprintf("  z = %7.2f m -> l^2 = %8.4f m^2", z, l))
        end
    end

    @testset "ON: reduces to (kappa*z)^2 near the wall, (C_s*Delta)^2 far away" begin
        # 1/l^2 = 1/(C_s*Delta)^2 + 1/(kappa*z)^2 -- Mason & Thomson with n = 2.
        z = 0.1                                    # (kappa z)^2 = 0.0016 << 40.96
        @test isapprox(sgs_mixing_length2(C_S2, Δ2, z, KARMAN, true),
                       (KARMAN*z)^2; rtol = 1e-3)
        z = 1.0e5                                  # (kappa z)^2 >> 40.96
        @test isapprox(sgs_mixing_length2(C_S2, Δ2, z, KARMAN, true),
                       CSΔ2; rtol = 1e-3)
        # and the exact Mason & Thomson value in between
        for z in (Z1, 20.0, 100.0)
            κz2 = (KARMAN*z)^2
            @test isapprox(sgs_mixing_length2(C_S2, Δ2, z, KARMAN, true),
                           CSΔ2*κz2/(CSΔ2 + κz2); rtol = 1e-12)
        end
    end

    @testset "never NaN, never negative, never above the unlimited value" begin
        # A NaN here seeds the whole eddy-viscosity field with one.
        for C2 in (0.0, C_S2), z in (-1.0, 0.0, 1e-300, Z1, 1e30)
            l = sgs_mixing_length2(C2, Δ2, z, KARMAN, true)
            @test isfinite(l)
            @test l >= 0.0
            @test l <= C2*Δ2 + 1e-12
        end
        # C_s = 0 and z = 0 is the 0/0 the guard exists for.
        @test sgs_mixing_length2(0.0, Δ2, 0.0, KARMAN, true) == 0.0
        # a negative wall distance (a warped mesh dipping below zmin) is
        # clamped, not squared into a positive damping length
        @test sgs_mixing_length2(C_S2, Δ2, -5.0, KARMAN, true) == 0.0
    end
end

say("\n=== the N^2 denominator floor ===")

@testset "SGS_THETA_FLOOR" begin
    # The floor exists so that a PERT run whose base state never got filled
    # answers N2 = 0 instead of dividing g by an O(0.1 K) perturbation. It has
    # to be unreachable for a real potential temperature and far above the
    # 1e-12 it replaced.
    m = match(r"const SGS_THETA_FLOOR\s*=\s*([0-9.eE+-]+)", _SGS_SRC)
    @test m !== nothing
    floor_val = parse(Float64, m.captures[1])
    @test floor_val >= 1.0            # catches theta' , which is O(0.1 K)
    @test floor_val < 100.0           # but nowhere near a real theta (~300 K)
    say(@sprintf("  SGS_THETA_FLOOR = %.3g K", floor_val))

    # and no site may still be using the old 1e-12 guard on that denominator
    @test !occursin(r"N2_val\s*=\s*abs\(θ_ref\)\s*>\s*1e-12", _SGS_SRC)
end

say("\n=== every N^2 site reads the TOTAL theta ===")

@testset "no compute_sgs_cache! site left on the perturbation" begin
    # The latent PERT bug: theta_ref must be uprimitive[...] + theta_base[ip].
    # A new closure copy-pasted from an old one is exactly how this comes back,
    # so assert on the source rather than on one code path.
    sites = collect(eachmatch(r"θ_ref\s*=\s*uprimitive\[[^\]]*\]([^\n]*)", _SGS_SRC))
    @test length(sites) >= 5
    for s in sites
        @test occursin("theta_base[ip]", s.captures[1])
    end
    say(@sprintf("  %d theta_ref sites, all carrying theta_base[ip]", length(sites)))

    # and the collocation derivative has to pick the base state up too,
    # otherwise dtheta/dz is still the perturbation's
    @test occursin("dθdζ += dψ[ii,m] * theta_base[connijk[iel,k,l,ii]]", _SGS_SRC)
    @test occursin("dθdη += dψ[ii,l] * theta_base[connijk[iel,k,ii]]", _SGS_SRC)
end

say("\n=== done ===\n")
