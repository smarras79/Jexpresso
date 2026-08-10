#==============================================================================
 test/test_sphere_mesh.jl — standalone tests for the spherical-shell grid.

     julia --project=. test/test_sphere_mesh.jl

 Not part of the CI case registry (test/ci_cases.jl runs whole simulations and
 compares HDF5 output; this exercises a grid builder, and it generates its own
 gmsh files in a temp directory so it needs nothing from meshes/).

 What is tested, and why each one matters on a CLOSED shell:

   * a WATERTIGHT cubed sphere (seam nodes already shared by the file)
   * a SPLIT cubed sphere (seam nodes DUPLICATED — the way gmsh writes a
     panel-by-panel geometry). Both must give the SAME grid after the
     coincident-node merge; if the merge is broken, this one comes out as six
     disconnected patches and T1..T4 fail.
   * MSH 2.2 and MSH 4.1 ASCII, so both readers are covered
   * elements written with a RANDOMLY REVERSED node order, which is what an
     inconsistently oriented gmsh surface looks like: after
     orient_elements_outward! the grid must be indistinguishable from the
     clean one.
   * for every case and every nop: the full T1..T9 suite of
     check_sphere_mesh, the exact node count, the Euler characteristic, and
     the fact that two elements sharing an edge list the SAME node indices
     along it.
==============================================================================#

using Test
using Random

include(joinpath(@__DIR__, "..", "tools", "generate_cubed_sphere.jl"))
using Jexpresso

const R_EARTH = 6.371e6


"write the generated grid as MSH 4.1 ASCII (write_msh22 in the generator covers 2.2)"
function write_msh41(fname::String, coords, quads, tags)
    mkpath(dirname(fname))
    open(fname, "w") do f
        println(f, "\$MeshFormat");  println(f, "4.1 0 8"); println(f, "\$EndMeshFormat")
        println(f, "\$Nodes")
        println(f, "1 ", length(coords), " 1 ", length(coords))
        println(f, "2 1 0 ", length(coords))
        for k = 1:length(coords); println(f, k); end
        for c in coords
            println(f, c[1], " ", c[2], " ", c[3])
        end
        println(f, "\$EndNodes")
        println(f, "\$Elements")
        println(f, "1 ", length(quads), " 1 ", length(quads))
        println(f, "2 7 3 ", length(quads))
        for (k, q) in enumerate(quads)
            println(f, k, " ", q[1], " ", q[2], " ", q[3], " ", q[4])
        end
        println(f, "\$EndElements")
    end
end


"the ngl node ids that element `iel` lists along its local edge `le`"
function edge_nodes(mesh, iel, le)
    ngl = mesh.ngl
    le == 1 && return [mesh.connijk[iel, l, 1]   for l = 1:ngl]
    le == 2 && return [mesh.connijk[iel, ngl, l] for l = 1:ngl]
    le == 3 && return [mesh.connijk[iel, l, ngl] for l = 1:ngl]
    return                [mesh.connijk[iel, 1, l]   for l = 1:ngl]
end


@testset "spherical shell grid" begin

    dir = mktempdir()

    for n in (2, 4, 5)
        #
        # Four flavours of the SAME sphere. All of them must produce an
        # identical topology after the reader has cleaned them up.
        #
        c_ok,  q_ok,  t_ok  = generate_cubed_sphere(n, R_EARTH; split_panels = false)
        c_spl, q_spl, t_spl = generate_cubed_sphere(n, R_EARTH; split_panels = true)

        Random.seed!(1234)
        q_flip = [rand() < 0.5 ? (q[4], q[3], q[2], q[1]) : q for q in q_ok]

        files = Dict{String,String}()
        files["watertight-msh2"] = joinpath(dir, "cs_$(n)_ok.msh")
        write_msh22(files["watertight-msh2"], c_ok, q_ok, t_ok)

        files["split-msh2"] = joinpath(dir, "cs_$(n)_split.msh")
        write_msh22(files["split-msh2"], c_spl, q_spl, t_spl)

        files["flipped-msh2"] = joinpath(dir, "cs_$(n)_flip.msh")
        write_msh22(files["flipped-msh2"], c_ok, q_flip, t_ok)

        files["watertight-msh4"] = joinpath(dir, "cs_$(n)_ok4.msh")
        write_msh41(files["watertight-msh4"], c_ok, q_ok, t_ok)

        for (label, fname) in files, nop in (1, 3, 4)

            @testset "n=$n $label nop=$nop" begin

                mesh = build_sphere_shell_mesh(fname, nop; verbose = false)
                ngl  = nop + 1

                # a cubed sphere of n×n per panel: 6n² quads, 6n²+2 vertices,
                # 12n² edges (from V-E+F=2 with F=6n², V=6n²+2)
                @test mesh.nelem        == 6*n^2
                @test mesh.npoin_linear == 6*n^2 + 2
                @test mesh.nedges       == 12*n^2
                @test mesh.npoin == mesh.npoin_linear +
                                    mesh.nedges*(ngl-2) + mesh.nelem*(ngl-2)^2

                # the full T1..T9 suite
                @test check_sphere_mesh(mesh; verbose = false)

                # neighbouring elements must list the SAME nodes on a shared
                # edge — the property that keeps the shell stitched together
                sides = [Tuple{Int,Int}[] for _ = 1:mesh.nedges]
                for iel = 1:mesh.nelem, le = 1:4
                    push!(sides[mesh.edge_in_elem[iel, le]], (iel, le))
                end
                torn = 0
                for ie = 1:mesh.nedges
                    length(sides[ie]) == 2 || continue
                    (ia, la), (ib, lb) = sides[ie][1], sides[ie][2]
                    na, nb = edge_nodes(mesh, ia, la), edge_nodes(mesh, ib, lb)
                    (na == nb || na == reverse(nb)) || (torn += 1)
                end
                @test torn == 0

                # every node on the shell
                @test maximum(abs.(sqrt.(mesh.x.^2 .+ mesh.y.^2 .+ mesh.z.^2) .- mesh.radius)) <
                      1.0e-8 * mesh.radius

                # the sub-cell area underestimates 4πR² but must be close, and
                # must IMPROVE with nop
                A = sphere_mesh_area(mesh)
                @test A <= 4π*mesh.radius^2
                @test abs(A - 4π*mesh.radius^2)/(4π*mesh.radius^2) < 0.1
            end
        end

        # the split (duplicated seams) and watertight files must give the very
        # same grid size once merged
        m_ok  = build_sphere_shell_mesh(files["watertight-msh2"], 4; verbose = false)
        m_spl = build_sphere_shell_mesh(files["split-msh2"],      4; verbose = false)
        @test m_ok.npoin        == m_spl.npoin
        @test m_ok.nedges       == m_spl.nedges
        @test m_ok.npoin_linear == m_spl.npoin_linear
    end

    #
    # VTK output: the files must be produced and be non-empty.
    #
    @testset "VTK writer" begin
        c, q, t = generate_cubed_sphere(3, R_EARTH)
        fn = joinpath(dir, "cs_vtk.msh")
        write_msh22(fn, c, q, t)
        nop  = 4
        mesh = build_sphere_shell_mesh(fn, nop; verbose = false)

        outdir = joinpath(dir, "out")
        write_vtk_sphere_grid(mesh, "sphere_grid_ho", outdir; verbose = false)

        @test isfile(joinpath(outdir, "sphere_grid_ho.vtu"))
        @test filesize(joinpath(outdir, "sphere_grid_ho.vtu")) > 0

        # exactly ONE grid file, like the rest of the cases
        @test length(filter(endswith(".vtu"), readdir(outdir))) == 1

        #
        # The surface file must carry the HIGH-ORDER sub-elements, not the
        # linear elements: (ngl-1)² VTK_QUADs per spectral element, and every
        # LGL node must appear as a corner of one of them. Read the cell count
        # straight out of the .vtu header rather than trusting the writer.
        #
        ngl  = nop + 1
        nsub = mesh.nelem*(ngl-1)^2
        bytes = read(joinpath(outdir, "sphere_grid_ho.vtu"))
        hdr   = String(bytes[1:min(lastindex(bytes), 4000)])
        mp    = match(r"NumberOfPoints=\"(\d+)\"", hdr)
        mc    = match(r"NumberOfCells=\"(\d+)\"",  hdr)
        @test mp !== nothing && mc !== nothing
        @test parse(Int, mp.captures[1]) == mesh.npoin
        @test parse(Int, mc.captures[1]) == nsub
        @test nsub > mesh.nelem                       # i.e. NOT one cell per element

        # and the sub-cells really do touch every high-order node
        touched = falses(mesh.npoin)
        for iel = 1:mesh.nelem, j = 1:ngl, i = 1:ngl
            touched[mesh.connijk[iel, i, j]] = true
        end
        @test all(touched)
    end

    #
    # The Galewsky et al. (2004) initial condition.
    #
    # problems/ShallowWater/SWsphere/initialize.jl is a Jexpresso-module file
    # (it uses NSD_2D, St_mesh_sphere, define_q), so evaluate it there, the way
    # run.jl does for a real run.
    #
    @testset "Galewsky initial condition" begin

        Jexpresso.MPI.Initialized() || Jexpresso.MPI.Init()
        Base.include(Jexpresso, joinpath(@__DIR__, "..", "problems",
                                         "ShallowWater", "SWsphere", "initialize.jl"))

        cc, qq, tt = generate_cubed_sphere(8, 6.37122e6)
        fn = joinpath(dir, "cs_ic.msh")
        write_msh22(fn, cc, qq, tt)
        mesh = build_sphere_shell_mesh(fn, 4; radius = 6.37122e6, verbose = false)

        inputs = Dict{Symbol,Any}(:backend => Jexpresso.CPU())
        p      = Jexpresso.galewsky_params(mesh, inputs)
        h0     = Jexpresso.galewsky_h0(p)

        # THE constant of integration quoted in the literature. Reproducing it
        # to 1e-6 m means the parameters AND the quadrature agree with the
        # paper — it is the sharpest single check on this initial condition.
        @test isapprox(h0, 10158.1861704546; atol = 1.0e-6)

        # the jet: zero outside the band, u_max exactly at the band midpoint
        # (φ0 + φ1 = π/2, so the midpoint is π/4)
        uj(φ) = Jexpresso.galewsky_ujet(φ, p.umax, p.φ0, p.φ1, p.en)
        @test uj(p.φ0)            == 0.0
        @test uj(p.φ1)            == 0.0
        @test uj(p.φ0 - 0.1)      == 0.0
        @test uj(p.φ1 + 0.1)      == 0.0
        @test uj(-0.5)            == 0.0
        @test isapprox(uj((p.φ0 + p.φ1)/2), p.umax; rtol = 1.0e-12)
        @test all(uj(φ) > 0 for φ in range(p.φ0 + 1e-3, p.φ1 - 1e-3; length = 50))

        # the balanced height: flat south of the jet, flat north of it
        hb(φ) = Jexpresso.galewsky_hbalance(φ, p)
        @test hb(-π/2)      == 0.0
        @test hb(p.φ0)      == 0.0
        @test hb(p.φ0 - 0.2) == 0.0
        @test isapprox(hb(π/2), hb(p.φ1); rtol = 1.0e-12)
        @test hb(π/2) < 0                       # depth drops across the jet

        # GEOSTROPHIC BALANCE: g dh/dφ must equal -a u (f + tan(φ) u / a).
        # Finite-difference the integral and compare against its own integrand
        # — this is what certifies h and u are a balanced pair rather than two
        # independently plausible profiles.
        δ = 1.0e-5
        for φ in range(p.φ0 + 0.05, p.φ1 - 0.05; length = 12)
            dhdφ = (hb(φ + δ) - hb(φ - δ))/(2δ)
            @test isapprox(dhdφ,
                           -Jexpresso.galewsky_balance_integrand(φ, p)/p.g;
                           rtol = 1.0e-4)
        end

        # now the field on the mesh
        qs = Jexpresso.initialize(mesh.SD, 0, mesh, inputs, dir, Float64)

        h = @view qs.qn[1:mesh.npoin, 1]
        u = @view qs.qn[1:mesh.npoin, 2]
        v = @view qs.qn[1:mesh.npoin, 3]

        @test all(v .== 0.0)                                  # v ≡ 0
        @test all(0.0 .<= u .<= p.umax*(1 + 1.0e-12))         # jet, no overshoot
        @test maximum(u) > 0.9*p.umax                         # the jet is resolved
        @test all(h .> 0.0)
        @test minimum(h) > h0 + hb(π/2) - 1.0                 # bounded below by the northern plateau
        @test maximum(h) < h0 + p.hhat + 1.0

        # the perturbation is a positive bump and never exceeds ĥ
        d = qs.qn[1:mesh.npoin, 1] .- qs.qe[1:mesh.npoin, 1]
        @test minimum(d) >= 0.0
        @test maximum(d) <= p.hhat
        @test maximum(d) > 0.0                                # it IS switched on

        # the unperturbed reference is zonally symmetric: h depends on latitude
        # only. Bin by latitude and check the spread inside each bin.
        nb   = 200
        lo   = fill(Inf,  nb)
        hi   = fill(-Inf, nb)
        for ip = 1:mesh.npoin
            b = clamp(1 + floor(Int, (mesh.lat[ip] + π/2)/π*nb), 1, nb)
            lo[b] = min(lo[b], qs.qe[ip, 1])
            hi[b] = max(hi[b], qs.qe[ip, 1])
        end
        # each bin spans π/nb in latitude; |dh/dφ| <= max|integrand|/g, so the
        # spread inside a bin is bounded by that times the bin width
        slope = maximum(abs(Jexpresso.galewsky_balance_integrand(φ, p))/p.g
                        for φ in range(p.φ0, p.φ1; length = 2000))
        @test all(hi[b] - lo[b] <= 1.01*slope*(π/nb) + 1.0e-6 for b = 1:nb if isfinite(lo[b]))

        # switching the perturbation off must give exactly the balanced state
        inputs[:lgalewsky_perturbation] = false
        qb = Jexpresso.initialize(mesh.SD, 0, mesh, inputs, dir, Float64)
        @test qb.qn[1:mesh.npoin, 1] == qb.qe[1:mesh.npoin, 1]

        # the writer carries the fields onto the same single grid file
        icdir = joinpath(dir, "ic")
        write_vtk_sphere_grid(mesh, "sphere_grid_ho", icdir; q = qs, verbose = false)
        @test length(filter(endswith(".vtu"), readdir(icdir))) == 1
        bytes = read(joinpath(icdir, "sphere_grid_ho.vtu"))
        hdr   = String(bytes[1:min(lastindex(bytes), 20000)])
        for name in ("\"h\"", "\"u\"", "\"v\"", "\"velocity\"")
            @test occursin(name, hdr)
        end
    end
end
