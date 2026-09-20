#---------------------------------------------------------------------------------
# test/test_pod_sphere.jl — POD wired to the real thing: the snapshot recorder,
# the field extractors, the mass matrix as inner product, and the writers, on the
# cubed sphere that ships with the SWsphere case.
#
#   julia --project=. test/test_pod_sphere.jl
#
# The mathematics of the decomposition is checked in test/test_pod.jl, which
# needs no package. THIS file checks the wiring, i.e. everything that only means
# something once there is a mesh under it:
#
#   1. THE WEIGHTS ARE THE SURFACE INTEGRAL. Σ_ip w_ip = 4πR². POD is only
#      optimal under the L² product, so if the weights ever stop being the mass
#      matrix — or start double-counting a node — the modes stop being POD modes
#      while still looking perfectly plausible. This is the cheap assertion that
#      catches it.
#
#   2. THE EXTRACTORS READ WHAT THEY SAY THEY READ. A rigid rotation about the
#      polar axis has u_λ = ΩR cos φ and u_φ = 0 exactly, which pins down both
#      the tangent-basis projection and the column order of every target.
#
#   3. RECORD → DECOMPOSE → WRITE, end to end, on an analytic travelling wave
#      whose decomposition is known: exactly two modes, of equal energy, whose
#      span reproduces every snapshot. Then the basis is read back off the disk
#      and has to be the one that was written — that round trip is the interface
#      a reduced-order model built on a previous run depends on.
#
#   4. THE RASTERIZER ON A CUBED SPHERE, where the element edges meet at panel
#      seams and eight nodes sit on cube corners: full coverage, no overshoot,
#      and a spherical harmonic recovered to the resolution of the grid.
#
# The fixture is the .msh shipped with the case, as in test/test_sphere_metrics.jl,
# so the test needs no gmsh binary.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics,
                 pod_settings, pod_recorder, pod_due, pod_record!, pod_finalize!,
                 pod_weights, pod_target, pod_extract!, pod_from_snapshots,
                 pod_project, pod_reconstruct, pod_truncation_error,
                 pod_rank_for_energy, pod_save, pod_load,
                 equirectangular_raster, St_pod
using PartitionedArrays, MPI
using Printf

const CASE_MSH = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere",
                          "cubed_sphere.msh")

# mod_inputs_user_inputs! reads three module globals that run.jl normally sets
# from the command line before it includes a case.
@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphere"
    user_input_file            = "test/test_pod_sphere.jl"
end

function shell_inputs(nop)
    inputs = Dict{Symbol,Any}(
        :lread_gmsh           => true,
        :gmsh_filename        => CASE_MSH,
        :nop                  => nop,
        :interpolation_nodes  => "lgl",
        :backend              => Jexpresso.CPU(),
        :lspherical_shell     => true,
        :sphere_radius        => 6.37122e6,
        :lproject_to_sphere   => true,
        :lgrid_only           => true,
        :sphere_metrics       => :cross_product,
    )
    Jexpresso.mod_inputs_user_inputs!(inputs, 1)   # rank 1 => no banner spam
    inputs[:gmsh_filename] = CASE_MSH
    inputs[:nop]           = nop
    return inputs
end

# the conservative state of a rigid rotation about the polar axis, φ = g·h₀
function rigid_rotation_state(mesh, Ω, φ0)
    npoin = Int(mesh.npoin)
    q = zeros(Float64, npoin, 4)
    for ip = 1:npoin
        x, y = mesh.coords[1,ip], mesh.coords[2,ip]
        q[ip,1] = φ0
        q[ip,2] = φ0*(-Ω*y)          # (Ω ẑ) × x
        q[ip,3] = φ0*( Ω*x)
        q[ip,4] = 0.0
    end
    return q
end


@testset verbose = true "POD on the spherical shell" begin

isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the SWsphere case")

with_mpi() do distribute

    inputs  = shell_inputs(4)
    mesh, _ = mod_mesh_mesh_driver(inputs, 1, distribute)
    metrics = build_sphere_metrics(mesh, inputs; verbose = false)
    npoin   = Int(mesh.npoin)
    R       = mesh.radius

    #-----------------------------------------------------------------------------
    @testset "(1) the weights are the surface integral" begin
        w = pod_weights(mesh, metrics)
        @test length(w) == npoin
        @test all(>(0), w)                                   # serial: every node owned
        @test abs(sum(w) - 4π*R^2)/(4π*R^2) < 1.0e-6
        @test w == Float64.(metrics.M[1:npoin])
    end

    #-----------------------------------------------------------------------------
    @testset "(2) the extractors" begin

        Ω, φ0, g = 1.0e-5, 9.80616*10000.0, 9.80616
        q  = rigid_rotation_state(mesh, Ω, φ0)
        ζ  = [Float64(ip) for ip = 1:npoin]                  # a marker field
        d1 = zeros(npoin, 1)
        d2 = zeros(npoin, 2)
        d4 = zeros(npoin, 4)

        pod_extract!(d1, pod_target(:h), q, ζ, mesh, g)
        @test maximum(abs, d1[:,1] .- φ0/g) < 1e-9

        pod_extract!(d1, pod_target(:phi), q, ζ, mesh, g)
        @test maximum(abs, d1[:,1] .- φ0) < 1e-9

        pod_extract!(d1, pod_target(:vorticity), q, ζ, mesh, g)
        @test d1[:,1] == ζ

        pod_extract!(d4, pod_target(:state), q, ζ, mesh, g)
        @test d4 == q

        # a rigid rotation is purely zonal, with u_λ = ΩR cos φ
        pod_extract!(d1, pod_target(:u), q, ζ, mesh, g)
        @test maximum(abs(d1[ip,1] - Ω*R*cos(mesh.lat[ip])) for ip = 1:npoin) < 1e-6*Ω*R
        pod_extract!(d1, pod_target(:v), q, ζ, mesh, g)
        @test maximum(abs, d1[:,1]) < 1e-6*Ω*R

        # …and the vector target carries the two of them, in that order
        pod_extract!(d2, pod_target(:velocity), q, ζ, mesh, g)
        @test maximum(abs(d2[ip,1] - Ω*R*cos(mesh.lat[ip])) for ip = 1:npoin) < 1e-6*Ω*R
        @test maximum(abs, d2[:,2]) < 1e-6*Ω*R

        @test_throws ErrorException pod_target(:nonsense)
        # the vorticity target cannot be fed a run that never computed one
        @test_throws ErrorException pod_extract!(d1, pod_target(:vorticity), q, nothing, mesh, g)
    end

    #-----------------------------------------------------------------------------
    @testset "(3) record → decompose → write → read back" begin

        outdir = mktempdir()
        nsnap  = 16
        T      = 8*86400.0

        pin = Dict{Symbol,Any}(:lpod             => true,
                               :pod_fields       => [:vorticity],
                               :pod_nsnapshots   => nsnap,
                               :pod_tstart       => 0.0,
                               :pod_nmodes_plot  => 4,
                               :pod_write_vtk    => true,
                               :pod_write_png    => true,
                               :pod_write_data   => true)

        set = pod_settings(pin, 0.0, T)
        @test set.lpod && set.fields == [:vorticity] && set.nsnap == nsnap

        rec = pod_recorder(pin, mesh, 0.0, T; verbose = false)
        @test rec !== nothing
        @test rec.nsnapmax == nsnap + 1
        @test rec.lvort                                       # ζ is needed, and it says so

        #
        # An analytic travelling wave: exactly two modes, in quadrature. The
        # phase advances by one FULL PERIOD over the nsnap+1 samples — not over
        # the window — so that Σcos² = Σsin² and Σcos·sin = 0 exactly over the
        # snapshots. Sampling both ends of one period instead duplicates the
        # zero phase and splits the degenerate pair by ~5 % for that reason
        # alone, which says nothing about the decomposition.
        #
        # Zonal wavenumber 3 and not 4: the cubed sphere is 4-fold symmetric
        # about the polar axis, so a wave-4 pattern sits in a fixed phase
        # relation to the panels — cos(4λ) peaking at panel centres and sin(4λ)
        # at the seams — and the quadrature then treats the two members of the
        # pair differently for a reason that is about the grid and not the flow.
        #
        m    = 3
        q    = rigid_rotation_state(mesh, 1.0e-5, 9.80616*10000.0)
        ζ    = zeros(npoin)
        Δt   = rec.dt/13                                      # a step that does not
        for k = 0:nsnap                                       # divide the sampling interval
            t = k*rec.dt
            θ = 2π*k/(nsnap + 1)
            for ip = 1:npoin
                ζ[ip] = cos(m*mesh.lon[ip] - θ)*cos(mesh.lat[ip])^2
            end
            @test pod_due(rec, t, Δt)
            pod_record!(rec, q, ζ, t, mesh)
            @test !pod_due(rec, t, Δt)                        # …and not twice
        end
        @test rec.nsnap == nsnap + 1
        @test !pod_due(rec, 10T, Δt)                          # the buffer is full

        Ps = pod_finalize!(rec, mesh, metrics, outdir; verbose = false)
        @test length(Ps) == 1
        P = Ps[1]

        @test P.name  == "vorticity"
        @test P.nsnap == nsnap + 1
        @test length(P.λ) == 2                                # rank 2, exactly
        # a degenerate PAIR. The residual split is the cubed sphere's quadrature
        # seeing cos(3λ) and sin(3λ) slightly differently, not the algorithm.
        @test P.λ[1] ≈ P.λ[2] rtol = 0.10
        @test P.cumenergy[2] ≈ 1.0 atol = 1e-10
        @test P.orthoerr < 1e-10
        @test pod_rank_for_energy(P, 0.99) == 2
        @test pod_truncation_error(P)[3] < 1e-6

        # the span reproduces the data it came from
        w = pod_weights(mesh, metrics)
        for k = 1:P.nsnap
            qk = pod_reconstruct(P, P.a[k,:])
            @test maximum(abs, qk[:,1] .- rec.X[1][:,1,k]) < 1e-9
        end
        @test maximum(abs, pod_project(P, w, rec.X[1][:,:,3]) .- P.a[3,:]) < 1e-9

        #--- what was written
        for f in ("pod_vorticity.vtu",
                  "pod_vorticity_spectrum.csv", "pod_vorticity_coefficients.csv",
                  "pod_vorticity.jld2",
                  "pod_vorticity_modes.png", "pod_vorticity_spectrum.png",
                  "pod_vorticity_coefficients.png", "pod_vorticity_mean.png",
                  "pod_vorticity_mode_001.png")
            @test isfile(joinpath(outdir, f))
            @test filesize(joinpath(outdir, f)) > 0
        end
        # one data row per mode, plus the three comment lines and the header
        @test length(readlines(joinpath(outdir, "pod_vorticity_spectrum.csv"))) == 4 + length(P.λ)
        @test length(readlines(joinpath(outdir, "pod_vorticity_coefficients.csv"))) == 2 + P.nsnap

        #--- the round trip a ROM depends on
        Q = pod_load(joinpath(outdir, "pod_vorticity.jld2"))
        @test Q isa St_pod
        @test Q.name == P.name && Q.comps == P.comps
        @test Q.λ == P.λ && Q.a == P.a && Q.Φ == P.Φ && Q.q̄ == P.q̄
        @test Q.total == P.total && Q.method == P.method && Q.lmean == P.lmean

        path = pod_save(P, joinpath(outdir, "again.jld2"))
        @test pod_load(path).λ == P.λ
        @test_throws ErrorException pod_load(joinpath(outdir, "no_such_basis.jld2"))

        #--- and the switch that turns all of it off
        @test pod_recorder(Dict{Symbol,Any}(), mesh, 0.0, T; verbose = false) === nothing
        @test pod_finalize!(nothing, mesh, metrics, outdir; verbose = false) == St_pod{Float64}[]
        @test !pod_due(nothing, 0.0, 1.0)

        rm(outdir; recursive = true, force = true)
    end

    #-----------------------------------------------------------------------------
    @testset "(4) a field that never moves is skipped, not fatal" begin

        # A case whose forcing is off integrates a state at rest, and the POD of
        # a constant is not a small answer — it is undefined. That must cost the
        # field, not the run.
        outdir = mktempdir()
        pin    = Dict{Symbol,Any}(:lpod           => true,
                                  :pod_fields     => [:h],
                                  :pod_nsnapshots => 4,
                                  :pod_write_vtk  => false,
                                  :pod_write_png  => false,
                                  :pod_write_data => false)
        rec = pod_recorder(pin, mesh, 0.0, 100.0; verbose = false)
        q   = rigid_rotation_state(mesh, 0.0, 9.80616*10000.0)
        for k = 0:4
            pod_record!(rec, q, nothing, k*rec.dt, mesh)
        end
        @test rec.nsnap == 5
        @test pod_finalize!(rec, mesh, metrics, outdir; verbose = false) == St_pod{Float64}[]
        rm(outdir; recursive = true, force = true)
    end

    #-----------------------------------------------------------------------------
    @testset "(5) the equirectangular raster on a cubed sphere" begin

        # a degree-4 sectoral harmonic, smooth and crossing every panel seam
        f = [cos(4*mesh.lon[ip])*cos(mesh.lat[ip])^4 for ip = 1:npoin]

        λ, φ, F = equirectangular_raster(f, mesh; nlon = 360, nlat = 180)
        @test size(F) == (360, 180)
        @test count(isnan, F) == 0                            # poles and seams included
        @test minimum(F) >= minimum(f) - 1e-12                # no overshoot: every
        @test maximum(F) <= maximum(f) + 1e-12                # pixel is a convex combination

        worst = 0.0
        for i = 1:360, j = 1:180
            λr, φr = λ[i]*π/180, φ[j]*π/180
            worst  = max(worst, abs(F[i,j] - cos(4λr)*cos(φr)^4))
        end
        @test worst < 0.05                                    # the grid's own resolution
    end
end

end # testset
