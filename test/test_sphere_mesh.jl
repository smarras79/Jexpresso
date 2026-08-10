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
using MPI

include(joinpath(@__DIR__, "..", "tools", "generate_cubed_sphere.jl"))
using Jexpresso

# build_sphere_shell_mesh -> println_rank -> get_mpi_comm needs a live MPI, and
# `using Jexpresso` from the precompile cache never evaluates the module body
# that would have called MPI.Init().
MPI.Initialized() || MPI.Init()

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

                # sphere_mesh_area is a crude POLYGONAL sum (the SEM quadrature
                # of the area is Σ M[ip], checked in the metrics as M4), so it
                # always underestimates 4πR² and is 12% low at the coarsest grid
                # here (n=2, nop=1). The bound only has to exclude nonsense.
                A = sphere_mesh_area(mesh)
                @test A <= 4π*mesh.radius^2
                @test abs(A - 4π*mesh.radius^2)/(4π*mesh.radius^2) < 0.25
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
        # u is C^∞ with every derivative vanishing at the edges, so it
        # UNDERFLOWS to exactly 0 close to them (exp(-1493) at φ0+1e-3, 1.8e-60
        # at φ0+0.01). Positivity is only meaningful in the core.
        @test all(uj(φ) > 0 for φ in range(p.φ0 + 0.15, p.φ1 - 0.15; length = 50))

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
        δ     = 1.0e-5
        slope = maximum(abs(Jexpresso.galewsky_balance_integrand(x, p))/p.g
                        for x in range(p.φ0, p.φ1; length = 2000))
        for φ in range(p.φ0 + 0.05, p.φ1 - 0.05; length = 12)
            dhdφ = (hb(φ + δ) - hb(φ - δ))/(2δ)
            # |dh/dφ| runs from 6.0e3 in the core to 2.9e-7 near the edges, so
            # the check is relative to the PEAK slope: a pointwise rtol would be
            # comparing quadrature noise against a vanishing quantity.
            @test isapprox(dhdφ,
                           -Jexpresso.galewsky_balance_integrand(φ, p)/p.g;
                           rtol = 1.0e-4, atol = 1.0e-6*slope)
        end

        # now the field on the mesh
        qs = Jexpresso.initialize(mesh.SD, 0, mesh, inputs, dir, Float64)

        # q.qout carries the primitives [h, u, v, w]; q.qn the conservative
        # state [φ, φu, φv, φw] with φ = g h
        h    = @view qs.qout[1:mesh.npoin, 1]
        umag = [sqrt(qs.qout[ip,2]^2 + qs.qout[ip,3]^2 + qs.qout[ip,4]^2)
                for ip = 1:mesh.npoin]

        @test all(0.0 .<= umag .<= p.umax*(1 + 1.0e-12))      # jet, no overshoot
        @test maximum(umag) > 0.9*p.umax                      # the jet is resolved
        @test all(h .> 0.0)
        @test minimum(h) > h0 + hb(π/2) - 1.0                 # bounded below by the northern plateau
        @test maximum(h) < h0 + p.hhat + 1.0

        # φ = g h, exactly
        @test all(isapprox(qs.qn[ip,1], p.g*qs.qout[ip,1]; rtol=1e-14) for ip = 1:mesh.npoin)

        # THE constraint: the initial momentum is tangential to the shell
        @test Jexpresso.sphere_normal_momentum(qs.qn, mesh; ivar=2) <
              1.0e-9 * maximum(abs, @view qs.qn[1:mesh.npoin, 2:4])

        # the perturbation is a positive bump on the geopotential, bounded by g·ĥ
        d = qs.qn[1:mesh.npoin, 1] .- qs.qe[1:mesh.npoin, 1]
        @test minimum(d) >= 0.0
        @test maximum(d) <= p.g*p.hhat
        @test maximum(d) > 0.0                                # it IS switched on

        # the unperturbed reference is zonally symmetric: h depends on latitude
        # only. Bin by latitude and check the spread inside each bin.
        nb   = 200
        lo   = fill(Inf,  nb)
        hi   = fill(-Inf, nb)
        for ip = 1:mesh.npoin
            b = clamp(1 + floor(Int, (mesh.lat[ip] + π/2)/π*nb), 1, nb)
            lo[b] = min(lo[b], qs.qe[ip, 1]/p.g)
            hi[b] = max(hi[b], qs.qe[ip, 1]/p.g)
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
        inputs[:lgalewsky_perturbation] = true

        # the writer carries the fields onto the same single grid file
        icdir = joinpath(dir, "ic")
        write_vtk_sphere_grid(mesh, "sphere_grid_ho", icdir; q = qs, verbose = false)
        @test length(filter(endswith(".vtu"), readdir(icdir))) == 1
        bytes = read(joinpath(icdir, "sphere_grid_ho.vtu"))
        hdr   = String(bytes[1:min(lastindex(bytes), 20000)])
        for name in ("\"phi\"", "\"phiu\"", "\"h\"", "\"u\"", "\"velocity\"",
                     "\"momentum_normal\"")
            @test occursin(name, hdr)
        end
    end

    #
    # The RHS: fluxes, sources and the Lagrange multiplier of
    # Marras, Kopera & Giraldo (2015), QJRMS 141: 1727-1739, Eq. (8)-(11).
    #
    @testset "SWE fluxes, sources and the Lagrange multiplier" begin

        Jexpresso.MPI.Initialized() || Jexpresso.MPI.Init()
        for f in ("initialize.jl", "user_flux.jl", "user_source.jl", "user_primitives.jl")
            Base.include(Jexpresso, joinpath(@__DIR__, "..", "problems",
                                             "ShallowWater", "SWsphere", f))
        end

        cc, qq, tt = generate_cubed_sphere(6, 6.37122e6)
        fn = joinpath(dir, "cs_rhs.msh")
        write_msh22(fn, cc, qq, tt)
        mesh = build_sphere_shell_mesh(fn, 4; radius = 6.37122e6, verbose = false)

        inputs = Dict{Symbol,Any}(:backend => Jexpresso.CPU())
        qs = Jexpresso.initialize(mesh.SD, 0, mesh, inputs, dir, Float64)
        r  = mesh.radius
        ω  = 7.292e-5
        g  = 9.80616

        F = zeros(4); G = zeros(4); H = zeros(4); S = zeros(4)

        for ip in (1, 7, 123, mesh.npoin ÷ 3, mesh.npoin)

            qv = qs.qn[ip, 1:4]
            φ, φu, φv, φw = qv
            u  = (φu/φ, φv/φ, φw/φ)
            x  = (mesh.x[ip], mesh.y[ip], mesh.z[ip])

            #----------------------------------------------------------- fluxes
            Jexpresso.user_flux!(F, G, H, mesh.SD, qv, qv, mesh,
                                 Jexpresso.CL(), Jexpresso.TOTAL(); neqs=4, ip=ip)

            # mass flux is the momentum itself
            @test (F[1], G[1], H[1]) == (φu, φv, φw)

            # momentum flux tensor T_ij = φ u_i u_j + (φ²/2) δ_ij: symmetric,
            # with the pressure exactly on the diagonal
            T = [F[2] G[2] H[2]; F[3] G[3] H[3]; F[4] G[4] H[4]]
            @test T ≈ T'
            for i = 1:3, j = 1:3
                ref = φ*u[i]*u[j] + (i == j ? 0.5*φ*φ : 0.0)
                @test isapprox(T[i,j], ref; rtol=1e-14, atol=1e-9)
            end

            #---------------------------------------------------------- sources
            Jexpresso.user_source!(S, qv, qv, mesh.npoin,
                                   Jexpresso.CL(), Jexpresso.TOTAL();
                                   neqs=4, x=x[1], y=x[2], z=x[3])

            @test S[1] == 0.0                                   # no mass source

            Smom = (S[2], S[3], S[4])
            f0   = 2.0*ω*x[3]/r                                 # 2ω sin(lat)

            # Coriolis part alone: -f (x × φu), f = 2ωz/r²
            fpar = 2.0*ω*x[3]/(r*r)
            cor  = (-fpar*(x[2]*φw - x[3]*φv),
                    -fpar*(x[3]*φu - x[1]*φw),
                    -fpar*(x[1]*φv - x[2]*φu))

            # ... is tangential, and perpendicular to the momentum
            @test abs(cor[1]*x[1] + cor[2]*x[2] + cor[3]*x[3]) <
                  1e-10*max(1.0, sum(abs, cor))*r
            @test abs(cor[1]*φu + cor[2]*φv + cor[3]*φw) <
                  1e-10*max(1.0, sum(abs, cor))*max(1.0, sqrt(φu^2+φv^2+φw^2))
            # ... and has the magnitude of the classical f₀ |φu|
            @test isapprox(sqrt(sum(c^2 for c in cor)),
                           abs(f0)*sqrt(φu^2 + φv^2 + φw^2); rtol=1e-12)

            # THE LAGRANGE MULTIPLIER. Whatever is left of the momentum source
            # after removing Coriolis must be μx with μ = -φ|u|²/r².
            lag  = (Smom[1]-cor[1], Smom[2]-cor[2], Smom[3]-cor[3])
            u2   = u[1]^2 + u[2]^2 + u[3]^2
            μref = -φ*u2/(r*r)
            for k = 1:3
                @test isapprox(lag[k], μref*x[k]; rtol=1e-12, atol=1e-12)
            end

            # it is purely NORMAL (parallel to x) and CENTRIPETAL (inward),
            # of magnitude φ|u|²/r — the force that bends a parcel of "mass" φ
            # moving at |u| around a circle of radius r
            @test μref <= 0.0
            @test isapprox(sqrt(sum(c^2 for c in lag)), φ*u2/r; rtol=1e-9, atol=1e-9)
            if u2 > 0
                @test lag[1]*x[1] + lag[2]*x[2] + lag[3]*x[3] < 0    # points inward
            end
        end

        #-------------------------------------------------- the projection P
        # P = I - xxᵀ/r² (paper Eq. 10-11) must: kill the normal component,
        # leave a tangential field untouched, and be idempotent.
        w = copy(qs.qn)
        @test Jexpresso.sphere_normal_momentum(w, mesh; ivar=2) <
              1.0e-9*maximum(abs, @view w[1:mesh.npoin, 2:4])

        d0 = Jexpresso.project_momentum_to_sphere!(w, mesh; ivar=2)
        @test d0 < 1.0e-9*maximum(abs, @view w[1:mesh.npoin, 2:4])   # nothing to remove
        @test w[1:mesh.npoin, 2:4] ≈ qs.qn[1:mesh.npoin, 2:4]        # ... so unchanged

        # now push a deliberate normal component in and check it is removed
        bad = copy(qs.qn)
        for ip = 1:mesh.npoin
            s = 0.37*sqrt(bad[ip,2]^2 + bad[ip,3]^2 + bad[ip,4]^2 + 1.0)/mesh.radius
            bad[ip,2] += s*mesh.x[ip]
            bad[ip,3] += s*mesh.y[ip]
            bad[ip,4] += s*mesh.z[ip]
        end
        @test Jexpresso.sphere_normal_momentum(bad, mesh; ivar=2) > 0.3

        removed = Jexpresso.project_momentum_to_sphere!(bad, mesh; ivar=2)
        @test removed > 0.3
        @test Jexpresso.sphere_normal_momentum(bad, mesh; ivar=2) <
              1.0e-9*maximum(abs, @view bad[1:mesh.npoin, 2:4])
        # the tangential part survived untouched
        @test bad[1:mesh.npoin, 2:4] ≈ qs.qn[1:mesh.npoin, 2:4]
        # idempotent
        again = copy(bad)
        Jexpresso.project_momentum_to_sphere!(again, mesh; ivar=2)
        @test again ≈ bad

        #-------------------------------------------- THE steady-state check
        #
        # The UNPERTURBED Galewsky jet is an exact steady solution, so the whole
        # right-hand side must vanish on it:
        #
        #     dq/dt = -∇ₛ·F(q) + S(q) = 0 .
        #
        # This is the test that certifies the formulation as a whole rather than
        # term by term: get the pressure term, the Coriolis sign, the Lagrange
        # multiplier or the manifold metrics wrong and the cancellation is
        # destroyed.
        #
        # The divergence here is the REAL SEM operator — sphere_rhs! with the
        # metric terms of sphere_metrics.jl — not a finite-difference stand-in.
        # It cannot be zero (the initial state is only represented to the
        # polynomial order of the grid), so what is checked is that it CONVERGES
        # under refinement, which is the statement that the operator is
        # consistent with the equations.
        #
        inputs2 = Dict{Symbol,Any}(:backend => Jexpresso.CPU(),
                                   :interpolation_nodes => Jexpresso.LGL(),
                                   :lgalewsky_perturbation => false)
        resid = Float64[]
        for nopr in (3, 5)
            meshr = build_sphere_shell_mesh(fn, nopr; radius = 6.37122e6, verbose = false)
            metr  = build_sphere_metrics(meshr, inputs2; verbose = false)
            @test check_sphere_metrics(meshr, metr; verbose = false)

            qr  = Jexpresso.initialize(meshr.SD, 0, meshr, inputs2, dir, Float64)
            spr = build_sphere_params(meshr, metr, inputs2; neqs = 4)
            RHS = zeros(meshr.npoin, 4)
            sphere_rhs!(RHS, qr.qn, qr.qe, meshr, metr, spr, Jexpresso.TOTAL())

            A = sum(metr.M)
            push!(resid, sqrt(sum(metr.M[ip]*RHS[ip,1]^2 for ip = 1:meshr.npoin)/A))
        end
        # raising the order must cut the residual hard (measured: ~20x from
        # nop=3 to nop=5 on this grid); a factor of 2 is a very loose floor that
        # a wrong equation could not clear
        @test resid[2] < 0.5*resid[1]

        #------------------------------------------------------- primitives
        up = zeros(4)
        for ip in (1, 55, mesh.npoin)
            qv = qs.qn[ip, 1:4]
            Jexpresso.user_primitives!(qv, qv, up, Jexpresso.TOTAL())
            @test up[1] == qv[1]
            @test isapprox(up[2], qv[2]/qv[1]; rtol=1e-15)

            uo = zeros(4)
            Jexpresso.user_uout!(ip, Jexpresso.TOTAL(), uo, qv, qv)
            @test isapprox(uo[1], qv[1]/g;      rtol=1e-14)   # h
            @test isapprox(uo[2], qv[2]/qv[1];  rtol=1e-15)   # u
        end
    end

    #
    # Advancing in time: the SSP-RK3 loop, the conservation properties it must
    # have, and the constraint it must hold.
    #
    @testset "time integration" begin

        Jexpresso.MPI.Initialized() || Jexpresso.MPI.Init()

        cc, qq, tt = generate_cubed_sphere(6, 6.37122e6)
        fn = joinpath(dir, "cs_tl.msh")
        write_msh22(fn, cc, qq, tt)
        mesh = build_sphere_shell_mesh(fn, 4; radius = 6.37122e6, verbose = false)

        base = Dict{Symbol,Any}(:backend => Jexpresso.CPU(),
                                :interpolation_nodes => Jexpresso.LGL(),
                                :SOL_VARS_TYPE => Jexpresso.TOTAL(),
                                :lgalewsky_perturbation => false,   # the STEADY state
                                :tinit => 0.0, :tend => 3600.0,
                                :ndiagnostics_outputs => 0,
                                :lwrite_initial => false,
                                :cfl => 0.35)

        metrics = build_sphere_metrics(mesh, base; verbose = false)

        # the mass matrix IS the SEM quadrature of the shell: Σ M = 4πR²
        @test isapprox(sum(metrics.M), 4π*mesh.radius^2; rtol = 1.0e-8)

        q  = Jexpresso.initialize(mesh.SD, 0, mesh, base, dir, Float64)
        sp = build_sphere_params(mesh, metrics, base; neqs = 4)

        # the CFL step is positive and finite
        Δt = sphere_cfl_dt(q.qn, mesh, metrics; cfl = 0.35)
        @test Δt > 0 && isfinite(Δt)

        mass0, ener0, drift0 = sphere_diagnostics(q.qn, mesh, metrics)
        @test drift0 < 1.0e-9*maximum(abs, @view q.qn[1:mesh.npoin, 2:4])

        h0 = copy(q.qn[1:mesh.npoin, 1])
        t  = sphere_time_loop!(mesh, metrics, sp, q, base, dir; verbose = false)

        @test isapprox(t, base[:tend]; rtol = 1.0e-12)

        mass1, ener1, drift1 = sphere_diagnostics(q.qn, mesh, metrics)

        # MASS is conserved to round-off: the scheme is conservative and a
        # closed manifold has no boundary flux to lose it through
        @test abs(mass1 - mass0)/mass0 < 1.0e-8

        # the flow stayed ON the shell — the whole point of the multiplier and
        # the per-stage projection
        @test drift1 < 1.0e-8*maximum(abs, @view q.qn[1:mesh.npoin, 2:4])

        # nothing went unphysical
        @test all(isfinite, q.qn[1:mesh.npoin, 1:4])
        @test minimum(@view q.qn[1:mesh.npoin, 1]) > 0

        # the balanced jet is a STEADY solution, so the drift away from it is
        # pure discretization error and must SHRINK when the order is raised
        err4 = maximum(abs.(q.qn[1:mesh.npoin,1] .- h0))

        mesh6 = build_sphere_shell_mesh(fn, 6; radius = 6.37122e6, verbose = false)
        met6  = build_sphere_metrics(mesh6, base; verbose = false)
        q6    = Jexpresso.initialize(mesh6.SD, 0, mesh6, base, dir, Float64)
        sp6   = build_sphere_params(mesh6, met6, base; neqs = 4)
        h06   = copy(q6.qn[1:mesh6.npoin, 1])
        sphere_time_loop!(mesh6, met6, sp6, q6, base, dir; verbose = false)
        err6  = maximum(abs.(q6.qn[1:mesh6.npoin,1] .- h06))

        @test err6 < err4

        # the filter must conserve mass exactly (mass-weighted DSS average)
        finp = merge(base, Dict{Symbol,Any}(:lfilter => true))
        spf  = build_sphere_params(mesh, metrics, finp; neqs = 4)
        @test spf.lfilter
        qf   = Jexpresso.initialize(mesh.SD, 0, mesh, finp, dir, Float64)
        mf0, _, _ = sphere_diagnostics(qf.qn, mesh, metrics)
        sphere_filter!(qf.qn, mesh, metrics, spf)
        mf1, _, _ = sphere_diagnostics(qf.qn, mesh, metrics)
        @test abs(mf1 - mf0)/mf0 < 1.0e-12
    end
end
