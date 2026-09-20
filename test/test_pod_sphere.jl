#---------------------------------------------------------------------------------
# test/test_pod_sphere.jl — POD wired to the real thing: the snapshot recorder,
# the field extractors, the mass matrix as inner product, and the writers, on the
# cubed sphere that ships with the SWsphere case.
#
#   julia --project=. test/test_pod_sphere.jl                 # serial
#   mpiexec -n 4 julia --project=. test/test_pod_sphere.jl     # and under MPI
#
# ANY RANK COUNT WORKS, INCLUDING 1, and the point of running it at more than
# one is that the parallel half of POD has invariants a serial run cannot see.
# The mathematics of the decomposition is checked in test/test_pod.jl and
# test/test_pod_benchmark.jl, neither of which needs a mesh. THIS file checks
# everything that only means something once there is one under it:
#
#   1. THE WEIGHTS ARE THE SURFACE INTEGRAL. Σ w = 4πR², summed over OWNED nodes
#      across all ranks. POD is only optimal under the L² product, so if the
#      weights ever stop being the mass matrix — or start double-counting a node
#      on a partition seam — the modes stop being POD modes while still looking
#      perfectly plausible. The failure is directional and this catches both
#      ends: too high if the ownership test goes missing (a shared node counted
#      once per rank), too low if the cross-rank assembly does.
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
#   4. THE DECOMPOSITION DOES NOT DEPEND ON THE PARTITION (nparts > 1 only, and
#      the reason this file takes an -n):
#        * every rank must come out with the SAME spectrum. The correlation
#          matrix is a sum over nodes completed by one Allreduce, and the K×K
#          eigenproblem is then solved redundantly, so agreement to round-off is
#          the statement that the reduction happened at all;
#        * a node held by several ranks must carry the SAME value of every mode
#          on all of them. That is where the GLOBAL sign convention lives: the
#          sign of a mode is fixed by its largest-magnitude entry, which sits on
#          one rank, and a rank that decided locally would hand back a mode
#          negated with respect to its neighbours — a discontinuity straight
#          through the partition seam, in a field that is supposed to be smooth;
#        * the reduced coordinates (pod_project) must reduce too, or every rank
#          returns its own partial inner product.
#
#   5. THE RASTERIZER ON A CUBED SPHERE, where the element edges meet at panel
#      seams and eight nodes sit on cube corners: full coverage, no overshoot,
#      and a spherical harmonic recovered to the resolution of the grid. Serial
#      only — a rank holds a piece of the sphere, and a map of a piece is not a
#      map (which is also why the PNG writers stand down under MPI and say so).
#
# The fixture is the .msh shipped with the case, as in test/test_sphere_metrics.jl
# and test/test_sphere_parallel.jl, so the test needs no gmsh binary.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics,
                 pod_settings, pod_recorder, pod_due, pod_record!, pod_finalize!,
                 pod_weights, pod_target, pod_targets, pod_extract!, pod_from_snapshots,
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

#
# Gather (global node id, value) from every rank and report the largest
# disagreement between two ranks that hold the SAME node — 0.0 when nothing is
# shared, which is what a one-rank run sees. Same helper, same reason, as
# test/test_sphere_parallel.jl.
#
function worst_shared_disagreement(gips, vals, comm)
    counts = MPI.Allgather(Int32(length(gips)), comm)
    allg   = MPI.Gatherv!(collect(Int64, gips), MPI.VBuffer(
                          Vector{Int64}(undef, sum(counts)), counts), 0, comm)
    allv   = MPI.Gatherv!(collect(Float64, vals), MPI.VBuffer(
                          Vector{Float64}(undef, sum(counts)), counts), 0, comm)
    worst  = 0.0
    if MPI.Comm_rank(comm) == 0
        seen = Dict{Int64,Float64}()
        for k in eachindex(allg)
            g = allg[k]; v = allv[k]
            if haskey(seen, g)
                s = max(abs(v), abs(seen[g]))
                worst = max(worst, s > 0 ? abs(v - seen[g])/s : abs(v - seen[g]))
            else
                seen[g] = v
            end
        end
    end
    return MPI.bcast(worst, 0, comm)
end


isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the SWsphere case")

with_mpi() do distribute

    comm   = MPI.COMM_WORLD
    rank   = MPI.Comm_rank(comm)
    nparts = MPI.Comm_size(comm)

    rank == 0 && @info "POD on the spherical shell, $nparts MPI rank(s)"

    # the global sum every rank needs for an inner product of its own
    reduce_sum!(A) = (MPI.Allreduce!(A, MPI.SUM, comm); A)

@testset verbose = true "POD on the spherical shell ($nparts rank(s))" begin

    inputs  = shell_inputs(4)
    mesh, _ = mod_mesh_mesh_driver(inputs, nparts, distribute)
    metrics = build_sphere_metrics(mesh, inputs; verbose = false)
    npoin   = Int(mesh.npoin)
    R       = mesh.radius

    #-----------------------------------------------------------------------------
    @testset "(1) the weights are the surface integral" begin
        w = pod_weights(metrics.M, mesh)
        @test length(w) == npoin
        @test all(>=(0), w)

        # Summed over OWNED nodes, across every rank. pod_weights has already
        # zeroed the mirrored ones, so this is a plain sum — that is the whole
        # point of it.
        area = MPI.Allreduce(sum(w), MPI.SUM, comm)
        aref = 4π*R^2
        rank == 0 && @printf(" #   area Σw = %.10e , 4πR² = %.10e , rel err %.2e\n",
                             area, aref, abs(area - aref)/aref)
        @test abs(area - aref)/aref < 1.0e-6

        if nparts == 1
            @test all(>(0), w)                                # nothing is mirrored
            @test w == Float64.(metrics.M[1:npoin])
        else
            # a partitioned run must have SOMETHING mirrored, or the test is not
            # exercising the ownership test it claims to
            @test MPI.Allreduce(count(iszero, w), MPI.SUM, comm) > 0
        end
    end

    #-----------------------------------------------------------------------------
    @testset "(2) the extractors" begin

        Ω, φ0, g = 1.0e-5, 9.80616*10000.0, 9.80616
        q  = rigid_rotation_state(mesh, Ω, φ0)
        ζ  = [Float64(ip) for ip = 1:npoin]                  # a marker field
        d1 = zeros(npoin, 1)
        d2 = zeros(npoin, 2)
        d4 = zeros(npoin, 4)

        SWVARS = ["phi", "phiu", "phiv", "phiw"]
        SWOUT  = ["h", "u", "v", "w"]
        tgt(f) = pod_targets([f], SWVARS, SWOUT, 4, true)[1]

        pod_extract!(d1, tgt(:h), q, nothing, ζ, mesh, g)
        @test maximum(abs, d1[:,1] .- φ0/g) < 1e-9

        # "phi" is a SOLUTION variable of this case, resolved by name
        pod_extract!(d1, tgt("phi"), q, nothing, ζ, mesh, g)
        @test maximum(abs, d1[:,1] .- φ0) < 1e-9

        pod_extract!(d1, tgt(:vorticity), q, nothing, ζ, mesh, g)
        @test d1[:,1] == ζ

        # :state stacks every solution variable into one vector target
        st = tgt(:state)
        @test st.ncomp == 4 && st.comps == SWVARS
        pod_extract!(d4, st, q, nothing, ζ, mesh, g)
        @test d4 == q

        # a rigid rotation is purely zonal, with u_λ = ΩR cos φ
        pod_extract!(d1, tgt(:u), q, nothing, ζ, mesh, g)
        @test maximum(abs(d1[ip,1] - Ω*R*cos(mesh.lat[ip])) for ip = 1:npoin) < 1e-6*Ω*R
        pod_extract!(d1, tgt(:v), q, nothing, ζ, mesh, g)
        @test maximum(abs, d1[:,1]) < 1e-6*Ω*R

        # …and the vector target carries the two of them, in that order
        pod_extract!(d2, tgt(:velocity), q, nothing, ζ, mesh, g)
        @test maximum(abs(d2[ip,1] - Ω*R*cos(mesh.lat[ip])) for ip = 1:npoin) < 1e-6*Ω*R
        @test maximum(abs, d2[:,2]) < 1e-6*Ω*R

        # an output variable the case derives, but does not integrate
        @test tgt("h").src === :q                     # :h is the shell's own derived field…
        @test tgt("w").src === :qout                  # …"w" is only an output variable
        @test_throws ErrorException pod_extract!(d1, tgt("w"), q, nothing, ζ, mesh, g)

        @test_throws ErrorException tgt(:nonsense)
        # a shell field is not offered to a flat case
        @test_throws ErrorException pod_targets([:vorticity], SWVARS, SWOUT, 4, false)
        # the vorticity target cannot be fed a run that never computed one
        @test_throws ErrorException pod_extract!(d1, tgt(:vorticity), q, nothing, nothing, mesh, g)
    end

    #-----------------------------------------------------------------------------
    @testset "(3) record → decompose → write → read back" begin

        # ONE directory, made by rank 0 and broadcast: mktempdir() on every rank
        # gives every rank its own, and the pieces of a .pvtu would then land in
        # four different places.
        outdir = MPI.bcast(rank == 0 ? mktempdir() : "", 0, comm)
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

        rec = pod_recorder(pin, mesh, 0.0, T; verbose = false,
                           qvars = ["phi","phiu","phiv","phiw"],
                           qoutvars = ["h","u","v","w"], neqs = 4,
                           nsd = 2, lshell = true)
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
            pod_record!(rec, q, t, mesh; ζ = ζ)
            @test !pod_due(rec, t, Δt)                        # …and not twice
        end
        @test rec.nsnap == nsnap + 1
        @test !pod_due(rec, 10T, Δt)                          # the buffer is full

        Ps = pod_finalize!(rec, mesh, metrics.M, outdir; verbose = false)
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
        # serial takes the SVD of the weighted snapshot matrix; a partitioned run
        # cannot (M^{1/2} is singular where the weights were zeroed) and must
        # fall back to the correlation matrix
        @test P.method === (nparts == 1 ? :svd : :snapshot)

        # the span reproduces the data it came from, node by node, on this rank
        w = pod_weights(metrics.M, mesh)
        for k = 1:P.nsnap
            qk = pod_reconstruct(P, P.a[k,:])
            @test maximum(abs, qk[:,1] .- rec.X[1][:,1,k]) < 1e-9
        end
        # …and the reduced coordinates are the projections, once the inner
        # product has been completed across ranks
        aproj = pod_project(P, w, rec.X[1][:,:,3]; reducer! = reduce_sum!)
        @test maximum(abs, aproj .- P.a[3,:]) < 1e-9

        #-------------------------------------------------------------------------
        # (4) the partition must not show
        #-------------------------------------------------------------------------
        if nparts > 1
            # every rank solved the same K×K eigenproblem, so every rank must
            # have the same spectrum
            λ0 = MPI.bcast(copy(P.λ), 0, comm)
            @test maximum(abs.(P.λ .- λ0)) <= 1e-13*maximum(λ0)
            E0 = MPI.bcast(copy(P.energy), 0, comm)
            @test maximum(abs.(P.energy .- E0)) < 1e-13

            # a node two ranks share must carry the same mode value on both —
            # this is the global sign convention, and a locally-decided sign
            # shows up here as a disagreement of exactly 2
            for i = 1:length(P.λ)
                dΦ = worst_shared_disagreement(@view(mesh.ip2gip[1:npoin]),
                                               @view(P.Φ[1:npoin, 1, i]), comm)
                rank == 0 && @printf(" #   mode %d: worst disagreement on a shared node = %.2e\n", i, dΦ)
                @test dΦ < 1.0e-12
            end
            dM = worst_shared_disagreement(@view(mesh.ip2gip[1:npoin]),
                                           @view(P.q̄[1:npoin, 1]), comm)
            @test dM < 1.0e-12
        end

        #-------------------------------------------------------------------------
        # what was written. Under MPI the .vtu becomes a .pvtu written in pieces,
        # the CSV and the basis are rank 0's, and the PNGs are not written at all
        # — the raster needs the whole sphere on one rank.
        #-------------------------------------------------------------------------
        MPI.Barrier(comm)
        if rank == 0
            # The spectrum and the coefficients are global, so rank 0 writes
            # them once; the modes are partitioned, so the basis is one file per
            # rank (see pod_write_data).
            basis = nparts > 1 ? "pod_vorticity_rank0000.jld2" : "pod_vorticity.jld2"
            files = String["pod_vorticity" * (nparts > 1 ? ".pvtu" : ".vtu"),
                           "pod_vorticity_spectrum.csv",
                           "pod_vorticity_coefficients.csv",
                           basis]
            nparts == 1 && append!(files, ["pod_vorticity_modes.png",
                                           "pod_vorticity_spectrum.png",
                                           "pod_vorticity_coefficients.png",
                                           "pod_vorticity_mean.png",
                                           "pod_vorticity_mode_001.png"])
            for f in files
                @test isfile(joinpath(outdir, f))
                @test filesize(joinpath(outdir, f)) > 0
            end
            nparts > 1 && @test !isfile(joinpath(outdir, "pod_vorticity_modes.png"))
            # one data row per mode, plus the three comment lines and the header
            @test length(readlines(joinpath(outdir, "pod_vorticity_spectrum.csv"))) == 4 + length(P.λ)
            @test length(readlines(joinpath(outdir, "pod_vorticity_coefficients.csv"))) == 2 + P.nsnap

            #--- the round trip a ROM depends on
            Q = pod_load(joinpath(outdir, basis); verbose = false)
            @test Q isa St_pod
            @test Q.name == P.name && Q.comps == P.comps
            @test Q.λ == P.λ && Q.a == P.a && Q.total == P.total
            @test Q.method == P.method && Q.lmean == P.lmean
            # rank 0's slice of the modes, which is what rank 0 wrote
            @test Q.Φ == P.Φ && Q.q̄ == P.q̄

            path = pod_save(P, joinpath(outdir, "again.jld2"))
            @test pod_load(path).λ == P.λ
            @test_throws ErrorException pod_load(joinpath(outdir, "no_such_basis.jld2"))
        end

        #--- and the switch that turns all of it off
        @test pod_recorder(Dict{Symbol,Any}(), mesh, 0.0, T; verbose = false) === nothing
        @test pod_finalize!(nothing, mesh, metrics.M, outdir; verbose = false) == St_pod{Float64}[]
        # …and POD refuses a run whose grid changes under it
        @test_throws ErrorException pod_recorder(Dict{Symbol,Any}(:lpod => true, :lamr => true),
                                                 mesh, 0.0, T; verbose = false,
                                                 qvars = ["phi"], neqs = 1, nsd = 2, lshell = true)

        MPI.Barrier(comm)
        rank == 0 && rm(outdir; recursive = true, force = true)
    end

    #-----------------------------------------------------------------------------
    @testset "(4b) a field that never moves is skipped, not fatal" begin

        # A case whose forcing is off integrates a state at rest, and the POD of
        # a constant is not a small answer — it is undefined. That must cost the
        # field, not the run. The test that it is skipped is collective: the
        # decision is reduced, so every rank must reach the same one.
        outdir = MPI.bcast(rank == 0 ? mktempdir() : "", 0, comm)
        pin    = Dict{Symbol,Any}(:lpod           => true,
                                  :pod_fields     => [:h],
                                  :pod_nsnapshots => 4,
                                  :pod_write_vtk  => false,
                                  :pod_write_png  => false,
                                  :pod_write_data => false)
        rec = pod_recorder(pin, mesh, 0.0, 100.0; verbose = false,
                           qvars = ["phi","phiu","phiv","phiw"],
                           qoutvars = ["h","u","v","w"], neqs = 4,
                           nsd = 2, lshell = true)
        q   = rigid_rotation_state(mesh, 0.0, 9.80616*10000.0)
        for k = 0:4
            pod_record!(rec, q, k*rec.dt, mesh)
        end
        @test rec.nsnap == 5
        @test pod_finalize!(rec, mesh, metrics.M, outdir; verbose = false) == St_pod{Float64}[]
        MPI.Barrier(comm)
        rank == 0 && rm(outdir; recursive = true, force = true)
    end

    #-----------------------------------------------------------------------------
    # Serial only: a rank holds a piece of the sphere, and a map of a piece is
    # not a map. This is the same reason the PNG writers stand down under MPI.
    #-----------------------------------------------------------------------------
    if nparts == 1
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
    else
        rank == 0 && @info "raster tests skipped: they are serial by nature (a rank holds a piece of the sphere)"
    end

end # testset

end # with_mpi
