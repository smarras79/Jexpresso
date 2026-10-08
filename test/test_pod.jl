#---------------------------------------------------------------------------------
# test/test_pod.jl — the POD kernel (src/kernel/rom/pod_core.jl) and the
# rasterizer its modes are drawn on (src/io/plotting/mesh_raster.jl).
#
#   julia test/test_pod.jl
#
# NO `using Jexpresso`, and no package instantiation: both files under test were
# written free of St_mesh, MPI and Plots precisely so that the mathematics can be
# checked on problems whose answer is known in closed form, in a bare Julia with
# Test and LinearAlgebra. The parts that DO need a mesh — the snapshot recorder,
# the field extractors, the writers — are exercised against the real cubed sphere
# in test/test_pod_sphere.jl.
#
# WHAT IS PINNED DOWN HERE, in the order it matters:
#
#   1. POD OF A KNOWN DECOMPOSITION. A field built as f₁(x)g₁(t) + f₂(x)g₂(t)
#      with (f_i,f_j) = δ_ij and Σ_k g₁g₂ = 0 has exactly two modes, and they
#      are f₁ and f₂ with energies ⟨g₁²⟩ and ⟨g₂²⟩. Both algorithms must return
#      exactly that — the modes, the eigenvalues, the coefficients, and nothing
#      in the third mode.
#
#   2. THE INNER PRODUCT IS REALLY THE MASS MATRIX. The test above is run with a
#      deliberately NON-UNIFORM weight, and the modes are checked to be
#      orthonormal under it and NOT under the Euclidean product. This is the
#      mistake the whole file exists to prevent: an off-the-shelf SVD of the
#      snapshot matrix gives the Euclidean answer, which on LGL nodes weights the
#      clustered element-edge nodes several times their share of the domain.
#
#   3. THE TRUNCATION ERROR IS THE ONE ADVERTISED. Σ_{i>r}λ_i/Σλ is the claim POD
#      is used for; here the actual L² error of every rank-r reconstruction is
#      measured and compared against it.
#
#   4. A TRAVELLING WAVE COMES OUT AS A DEGENERATE PAIR. cos(m(λ-ct)) has λ₁ = λ₂
#      and coefficients in quadrature — a circle in the (a₁,a₂) plane. This is
#      the structure the SWsphere deck tells the user to look for in the Galewsky
#      jet, so it is worth pinning down that the code reproduces it in the case
#      where it is exactly true.
#
#   5. THE RASTERIZER IS AN INTERPOLATION, not a smoothing. It is EXACT on fields
#      linear in latitude (barycentric interpolation reproduces linear functions),
#      it never overshoots the nodal data, it leaves no hole in the canvas, and it
#      is continuous across the dateline — the four things a plot of a POD mode
#      depends on if the picture is to be of the mode rather than of the plotter.
#---------------------------------------------------------------------------------

using Test
using LinearAlgebra
using Random
using Printf

include(joinpath(@__DIR__, "..", "src", "kernel", "rom", "pod_core.jl"))
include(joinpath(@__DIR__, "..", "src", "io", "plotting", "mesh_raster.jl"))


#
# A structured longitude-latitude grid dressed as a spectral element mesh: the
# rasterizer only ever asks for lon, lat, connijk, nelem and ngl. Periodic in λ,
# so the last element of each row straddles the dateline and the wrap path is
# exercised by every test that uses it.
#
function latlon_mesh(nelλ::Int, nelφ::Int, ngl::Int; latmax::Float64 = 90.0)
    nx    = nelλ*(ngl-1)                 # λ wraps, so no repeated last column
    ny    = nelφ*(ngl-1) + 1
    npoin = nx*ny
    lon   = zeros(npoin)
    lat   = zeros(npoin)
    node(i, j) = (mod(i-1, nx))*ny + j
    for i = 1:nx, j = 1:ny
        ip = node(i, j)
        lon[ip] = -π + 2π*(i-1)/nx
        lat[ip] = (-latmax + 2*latmax*(j-1)/(ny-1))*π/180
    end
    nelem   = nelλ*nelφ
    connijk = zeros(Int, nelem, ngl, ngl)
    iel = 0
    for eλ = 1:nelλ, eφ = 1:nelφ
        iel += 1
        for i = 1:ngl, j = 1:ngl
            connijk[iel,i,j] = node((eλ-1)*(ngl-1) + i, (eφ-1)*(ngl-1) + j)
        end
    end
    return (lon = lon, lat = lat, connijk = connijk,
            nelem = nelem, ngl = ngl, npoin = npoin)
end

# the discrete L² product the whole file is about
minner(w, u, v) = sum(w[ip]*u[ip]*v[ip] for ip in eachindex(w))


@testset verbose = true "POD" begin

#---------------------------------------------------------------------------------
@testset "a known decomposition, recovered exactly ($method)" for method in (:svd, :snapshot)

    Random.seed!(20260920)
    npoin, K = 400, 16

    w  = 0.1 .+ rand(npoin)                       # a NON-uniform mass matrix
    f1 = randn(npoin); f1 ./= sqrt(minner(w, f1, f1))
    f2 = randn(npoin); f2 .-= minner(w, f1, f2).*f1; f2 ./= sqrt(minner(w, f2, f2))
    q̄  = randn(npoin)

    # a full period, so Σg₁ = Σg₂ = Σg₁g₂ = 0 exactly
    g1 = [cos(2π*k/K)     for k = 0:K-1]
    g2 = [0.5*sin(2π*k/K) for k = 0:K-1]

    X = Array{Float64}(undef, npoin, 1, K)
    for k = 1:K, ip = 1:npoin
        X[ip,1,k] = q̄[ip] + f1[ip]*g1[k] + f2[ip]*g2[k]
    end

    P = pod_from_snapshots(X, w, collect(0.0:K-1); name = "test", method = method)

    @test P.method  == method
    @test P.npoin   == npoin
    @test P.nsnap   == K
    @test length(P.λ) == 2                        # rank 2, and not one mode more

    #--- the eigenvalues are the temporal variances of the two signals
    @test P.λ[1] ≈ sum(abs2, g1)/K   atol = 1e-12
    @test P.λ[2] ≈ sum(abs2, g2)/K   atol = 1e-12
    @test P.total ≈ (sum(abs2, g1) + sum(abs2, g2))/K  atol = 1e-12
    @test sum(P.energy) ≈ 1.0        atol = 1e-12
    @test P.cumenergy[end] ≈ 1.0     atol = 1e-12

    #--- the mean, and the modes (up to the sign POD does not fix)
    @test maximum(abs, P.q̄[:,1] .- q̄) < 1e-12
    @test min(maximum(abs, P.Φ[:,1,1] .- f1), maximum(abs, P.Φ[:,1,1] .+ f1)) < 1e-12
    @test min(maximum(abs, P.Φ[:,1,2] .- f2), maximum(abs, P.Φ[:,1,2] .+ f2)) < 1e-12

    #--- the sign CONVENTION: the largest-magnitude entry of every mode is
    #    positive, so that re-running a case does not invert the colour of a plot
    for i = 1:2
        @test P.Φ[argmax(abs.(P.Φ[:,1,i])), 1, i] > 0
    end

    #--- the coefficients are the projections, and they reproduce the data
    for k = 1:K
        q = pod_reconstruct(P, P.a[k,:])
        @test maximum(abs, q[:,1] .- X[:,1,k]) < 1e-12
    end
    @test maximum(abs, pod_project(P, w, X[:,:,5]) .- P.a[5,:]) < 1e-12

    #--- (2) the inner product really is the weighted one
    @test P.orthoerr < 1e-12
    @test minner(w, P.Φ[:,1,1], P.Φ[:,1,2]) ≈ 0.0 atol = 1e-12
    @test minner(w, P.Φ[:,1,1], P.Φ[:,1,1]) ≈ 1.0 atol = 1e-12
    # …and NOT the Euclidean one, which is what an off-the-shelf SVD would give:
    # if this ever starts passing, the weights have quietly stopped being used.
    @test abs(dot(P.Φ[:,1,1], P.Φ[:,1,1]) - 1.0) > 1e-3

    @test pod_rank_for_energy(P, 0.5)  == 1
    @test pod_rank_for_energy(P, 0.99) == 2
end

#---------------------------------------------------------------------------------
@testset "the truncation error is the advertised one" begin

    Random.seed!(7)
    npoin, K, R = 300, 24, 8
    w = 0.5 .+ rand(npoin)

    X = zeros(npoin, 1, K)
    for m = 1:R                                   # R structures of decaying amplitude
        f = randn(npoin)
        g = randn(K)
        for k = 1:K, ip = 1:npoin
            X[ip,1,k] += (0.5^m)*f[ip]*g[k]
        end
    end

    P = pod_from_snapshots(X, w, collect(1.0:K); name = "decay")
    e = pod_truncation_error(P)

    @test length(e) == length(P.λ) + 1
    @test e[1] ≈ 1.0 atol = 1e-12                 # keeping the mean alone
    @test issorted(e; rev = true)                 # more modes never means more error

    # the MEASURED error of every rank-r reconstruction, against the prediction
    num = zeros(length(P.λ) + 1)
    den = 0.0
    for k = 1:K
        d = X[:,1,k] .- P.q̄[:,1]
        den += minner(w, d, d)
        for r = 0:length(P.λ)
            q = pod_reconstruct(P, P.a[k,1:r])
            res = X[:,1,k] .- q[:,1]
            num[r+1] += minner(w, res, res)
        end
    end
    #
    # Compared as SQUARED energies. Eq. (3) is a statement about energies, and
    # the square root of it is infinitely steep at zero: at full rank the tail
    # is a round-off residual of ~1e-16, whose square root is 1e-8 while the
    # measured error is 1e-15. That gap is the sqrt, not a disagreement.
    #
    for r = 0:length(P.λ)
        @test num[r+1]/den ≈ e[r+1]^2 atol = 1e-13
    end
end

#---------------------------------------------------------------------------------
@testset "a travelling wave is a degenerate pair in quadrature" begin

    # ζ = cos(m(λ - ct)) on a lon-lat grid: an exactly rank-2 field whose two
    # modes must have EQUAL energy — this is the signature the SWsphere deck
    # tells the user to look for in the Galewsky jet.
    msh   = latlon_mesh(8, 4, 4)
    npoin = msh.npoin
    K, m  = 24, 4
    w     = cos.(msh.lat) .+ 0.1                  # a stand-in for the area element

    X = Array{Float64}(undef, npoin, 1, K)
    for k = 1:K
        θ = 2π*(k-1)/K
        for ip = 1:npoin
            X[ip,1,k] = cos(m*msh.lon[ip] - θ)
        end
    end

    P = pod_from_snapshots(X, w, collect(0.0:K-1); name = "wave")

    @test length(P.λ) == 2
    @test P.λ[1] ≈ P.λ[2] rtol = 1e-8             # degenerate, hence a PAIR
    @test P.cumenergy[2] ≈ 1.0 atol = 1e-10

    # in quadrature: the phase portrait is a circle of constant radius, i.e. the
    # pair is one travelling structure and not two standing ones
    rad = [sqrt(P.a[k,1]^2 + P.a[k,2]^2) for k = 1:K]
    @test maximum(rad) - minimum(rad) < 1e-8*maximum(rad)
    @test abs(sum(P.a[:,1].*P.a[:,2])) < 1e-8*sum(abs2, P.a[:,1])
end

#---------------------------------------------------------------------------------
@testset "vector-valued targets are decomposed jointly" begin

    Random.seed!(3)
    npoin, K = 200, 12
    w = fill(1.0/npoin, npoin)

    # one structure per component, sharing ONE temporal coefficient: a single
    # mode of the PAIR, which is the point of a vector target.
    f = randn(npoin, 2)
    g = [sin(2π*k/K) for k = 0:K-1]
    X = Array{Float64}(undef, npoin, 2, K)
    for k = 1:K, c = 1:2, ip = 1:npoin
        X[ip,c,k] = f[ip,c]*g[k]
    end

    P = pod_from_snapshots(X, w, collect(1.0:K); name = "velocity",
                           comps = ["u", "v"], lmean = false)

    @test P.ncomp == 2
    @test P.comps == ["u", "v"]
    @test length(P.λ) == 1                        # ONE mode, with two components
    @test P.energy[1] ≈ 1.0 atol = 1e-12
    @test P.orthoerr < 1e-12
    # the energy of the vector is the sum over its components
    @test P.total ≈ (sum(abs2, g)/K)*(minner(w, f[:,1], f[:,1]) + minner(w, f[:,2], f[:,2])) rtol = 1e-10
    for k = 1:K
        q = pod_reconstruct(P, P.a[k,:])
        @test maximum(abs, q .- X[:,:,k]) < 1e-12
    end
end

#---------------------------------------------------------------------------------
@testset "what POD refuses to invent" begin

    Random.seed!(11)
    npoin, K = 120, 10
    w = rand(npoin) .+ 0.5
    X = randn(npoin, 1, K)

    # K snapshots span at most K dimensions, however fine the grid, and one of
    # them goes to the mean
    P = pod_from_snapshots(X, w, collect(1.0:K))
    @test length(P.λ) <= K - 1

    # the cap is honoured, and the coefficients follow it
    P3 = pod_from_snapshots(X, w, collect(1.0:K); nmodes = 3)
    @test length(P3.λ) == 3
    @test size(P3.a) == (K, 3)
    @test size(P3.Φ) == (npoin, 1, 3)
    @test sum(P3.energy) < 1.0                    # energy is measured against the FULL trace
    @test P3.total ≈ P.total rtol = 1e-12

    # one snapshot is not a decomposition, and a nonsense method is not a method
    @test_throws ErrorException pod_from_snapshots(X[:,:,1:1], w, [0.0])
    @test_throws ErrorException pod_from_snapshots(X, w, collect(1.0:K); method = :nonsense)
    @test_throws ErrorException pod_from_snapshots(X, w[1:10], collect(1.0:K))
    # :svd is serial only — under MPI the weights of the mirrored nodes are zero
    # and M^{1/2} cannot be inverted
    @test_throws ErrorException pod_from_snapshots(X, w, collect(1.0:K);
                                                   method = :svd, lparallel = true)
    # a rank-deficient weight vector is still fine through the correlation matrix
    w0 = copy(w); w0[1:20] .= 0.0
    P0 = pod_from_snapshots(X, w0, collect(1.0:K); method = :auto, lparallel = true)
    @test P0.method == :snapshot
    @test P0.orthoerr < 1e-10
end

#---------------------------------------------------------------------------------
@testset "equirectangular raster" begin

    msh = latlon_mesh(24, 12, 5)
    nlon, nlat = 360, 180

    #--- the canvas: pixel CENTRES, covering the globe exactly once
    λ, φ = equirectangular_grid(nlon, nlat)
    @test length(λ) == nlon && length(φ) == nlat
    @test λ[1] ≈ -179.5 && λ[end] ≈ 179.5
    @test φ[1] ≈  -89.5 && φ[end] ≈  89.5
    @test_throws ErrorException equirectangular_grid(1, 1)

    #--- (5a) EXACT on a field linear in latitude. Barycentric interpolation
    #    reproduces linear functions, so anything above round-off here is a bug
    #    in the geometry, not interpolation error. Latitude and not longitude
    #    because λ itself is discontinuous at the dateline while the field is not.
    flin = [2.0 + 3.0*msh.lat[ip] for ip = 1:msh.npoin]
    λ, φ, F = equirectangular_raster(flin, msh; nlon = nlon, nlat = nlat)
    @test size(F) == (nlon, nlat)
    @test count(isnan, F) == 0                    # no holes anywhere, poles included
    worst = 0.0
    for i = 1:nlon, j = 1:nlat
        worst = max(worst, abs(F[i,j] - (2.0 + 3.0*φ[j]*π/180)))
    end
    @test worst < 1e-12

    #--- (5b) a smooth wave crossing the dateline: second-order accurate, and
    #    NEVER outside the range of the nodal data
    fw(λr, φr) = sin(2λr)*cos(φr)^2 + 0.5*sin(3φr)
    f = [fw(msh.lon[ip], msh.lat[ip]) for ip = 1:msh.npoin]
    λ, φ, F = equirectangular_raster(f, msh; nlon = nlon, nlat = nlat)
    @test count(isnan, F) == 0
    @test minimum(F) >= minimum(f) - 1e-12        # a convex combination cannot
    @test maximum(F) <= maximum(f) + 1e-12        # overshoot the data
    worst = 0.0
    for i = 1:nlon, j = 1:nlat
        worst = max(worst, abs(F[i,j] - fw(λ[i]*π/180, φ[j]*π/180)))
    end
    @test worst < 0.02                            # ~5e-3 at this resolution
    # the dateline: the two end columns are one pixel apart across the seam and
    # must be as accurate as anywhere else. A seam drawn wrong is O(1) here.
    for j = 1:nlat
        @test abs(F[1,j]    - fw(λ[1]*π/180,    φ[j]*π/180)) < 0.02
        @test abs(F[nlon,j] - fw(λ[nlon]*π/180, φ[j]*π/180)) < 0.02
    end

    #--- (5c) the polar caps. A grid that stops at ±85° leaves the caps
    #    uncovered; the fallback fills them from the nearest node in 3-D, and
    #    switching it off leaves them visibly NaN rather than silently wrong.
    capped = latlon_mesh(16, 8, 4; latmax = 85.0)
    fc = [cos(capped.lat[ip]) for ip = 1:capped.npoin]
    _, _, Fgap = equirectangular_raster(fc, capped; nlon = 180, nlat = 90, lfill_gaps = false)
    @test count(isnan, Fgap) > 0
    _, φf, Ffill = equirectangular_raster(fc, capped; nlon = 180, nlat = 90, lfill_gaps = true)
    @test count(isnan, Ffill) == 0
    @test all(isfinite, Ffill)
    @test maximum(Ffill) <= maximum(fc) + 1e-12
    # the cap takes the value of the nearest node, i.e. the one at ±85°
    @test Ffill[1,1] ≈ cos(85*π/180) atol = 1e-6
end

end # testset POD
