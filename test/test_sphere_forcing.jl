#---------------------------------------------------------------------------------
# test/test_sphere_forcing.jl — the Scott & Polvani stochastic forcing operator,
# src/kernel/operators/sphere_forcing.jl, on the cubed sphere that ships with
# the SWsphere case.
#
#   julia --project=. test/test_sphere_forcing.jl                 # serial
#   mpiexec -n 4 julia --project=. test/test_sphere_forcing.jl     # and under MPI
#
# ANY RANK COUNT WORKS, INCLUDING 1.
#
# WHAT IS CHECKED, and why each one is worth a test:
#
#   1. THE BASIS IS ORTHONORMAL ON THIS GRID. ∫Y_j Y_k dΩ = a² δ_jk, measured
#      with the SEM mass matrix. This is two statements at once: that the
#      normalised associated Legendre recursion and the √2 on the m ≠ 0 real
#      harmonics are right, and that the GRID RESOLVES THE FORCED BAND at all.
#      The second is not decoration — force at an n the grid cannot represent
#      and the "forcing" is an aliasing pattern, which looks like turbulence and
#      is not. The check is run at the band the decks actually use, n_f = 24.
#
#      The √2 is the one that fails silently: without it the zonal (m = 0) modes
#      carry twice the variance of every other mode, so the forcing acquires a
#      systematic zonal bias — on a rotating planet, a bias towards precisely
#      the jets the case exists to measure.
#
#   2. THE ROUND TRIP. ζ(u_F) = ζ_F: build a known random vorticity field in the
#      band, let the operator turn it into u_F through ψ = ∇ₛ⁻²ζ and u = n̂ × ∇ₛψ,
#      then take the DISCRETE curl of u_F with the solver's own diagnostic and
#      compare. One number that exercises the harmonics, the -n(n+1)/a²
#      eigenvalue of the surface Laplacian, the assembled gradient, the cross
#      product and every sign in the chain. Any of them wrong and this misses.
#
#   3. TANGENCY. n̂ × g is orthogonal to n̂ identically, so u_F must be on the
#      shell to round-off whatever the discrete gradient returned. If this ever
#      drifts, the Lagrange projection is silently eating part of the forcing.
#
#   4. THE INJECTION RATE IS ε₀, AND IT IS SPECIFIC. The amplitude is the root
#      of a quadratic (Eq. 6 of the operator) chosen so the step injects exactly
#      the prescribed rate. ε₀ is energy per unit MASS per unit time — the decks
#      build it as 2ν_l·½U²/f_up — so the quadratic has to be normalised by
#      ∫φ dΩ. Getting that wrong is out by ~10²⁰ on a giant planet and the run
#      looks like it is still at rest, which is exactly the symptom this case
#      was stuck on before the operator existed. Checked twice: against the
#      operator's own bookkeeping, and by integrating from rest and watching
#      Ê(t) = ∫½φ|u|²dΩ / ∫φdΩ track ε₀·t.
#
#   5. EVERY RANK DRAWS THE SAME FIELD (nparts > 1 only). c_{nm} describes ONE
#      global field; the ranks stay in lockstep by replaying the same stream
#      rather than by communicating. If that ever diverges each rank forces its
#      own partition with its own field and the seams become the answer.
#---------------------------------------------------------------------------------
using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics, build_sphere_params,
                 build_sphere_forcing, sphere_forcing_step!, sphere_forcing_apply!,
                 sphere_relative_vorticity!, _forcing_stream_kernel!
using PartitionedArrays, MPI, Printf

@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphere"
    user_input_file            = "test_sphere_forcing"
end

const MSH = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere", "cubed_sphere.msh")
const NF  = 24            # the band the decks force at
const DN  = 4
const EPS = 2.5e-11       # a specific energy input rate [m²/s³]
const DT  = 100.0

inputs = Dict{Symbol,Any}(
    :lspherical_shell => true, :lread_gmsh => true, :gmsh_filename => MSH,
    :nop => 5, :interpolation_nodes => "lgl",
    :sphere_radius => 6.37122e6, :cubed_sphere_map => :none,
    :sphere_metrics => :curl_invariant, :lproject_to_sphere => true,
    :lcheck_grid => false, :lstop_on_bad_grid => false,
    :backend => Jexpresso.CPU(), :nvars => 4, :neqs => 4,
    :lfilter => false, :lvisc => false,
    :lsphere_forcing => true, :forcing_epsilon => EPS,
    :forcing_nf => NF, :forcing_dn => DN, :forcing_tau => 0.0,
    :forcing_seed => 7, :forcing_normalize => true,
    :rayleigh_friction => 0.0, :radiative_relaxation => 0.0,
)
Jexpresso.mod_inputs_user_inputs!(inputs, 1)

const R = Ref{Any}()
with_mpi() do distribute
    nparts  = MPI.Comm_size(MPI.COMM_WORLD)
    m, _    = mod_mesh_mesh_driver(inputs, nparts, distribute)
    mt      = build_sphere_metrics(m, inputs; verbose = false)
    R[]     = (m, mt, build_sphere_params(m, mt, inputs; neqs = 4))
end
const mesh, metrics, sp = R[]

const comm   = MPI.COMM_WORLD
const rank   = MPI.Comm_rank(comm)
const nparts = MPI.Comm_size(comm)
const npoin  = Int(mesh.npoin)
const a      = Float64(mesh.radius)
const M      = metrics.M

# the mass matrix with non-owned nodes zeroed, so an integral is counted once
# however many ranks hold the node
const w = let v = copy(M)
    if nparts > 1
        for ip = 1:npoin
            mesh.gip2owner[ip] == rank || (v[ip] = 0.0)
        end
    end
    v
end
gsum(x) = nparts > 1 ? MPI.Allreduce(x, MPI.SUM, comm) : x
gmax(x) = nparts > 1 ? MPI.Allreduce(x, MPI.MAX, comm) : x

const fc   = build_sphere_forcing(mesh, metrics, sp, inputs; Δt = DT, verbose = false)
const φref = 9.80616*10000.0
const q    = let z = zeros(Float64, npoin, 5)
    for ip = 1:npoin; z[ip,1] = φref; end
    z
end

rank == 0 && @printf("\n# test_sphere_forcing.jl : npoin = %d, nparts = %d, n ∈ [%d,%d], %d modes\n",
                     npoin, nparts, fc.nlo, fc.nhi, fc.nmode)

@testset "sphere_forcing" begin

#--- 1. the real harmonic basis is orthonormal on this grid --------------------
@testset "harmonic basis" begin
    nm = fc.nmode
    # unit-amplitude weights: the same 1 / √2 convention build_sphere_forcing
    # folds into fc.w, with the inverse Laplacian left out
    wunit = let v = zeros(Float64, nm), idx = 0
        for m = 0:fc.nhi, n = max(fc.nlo, m):fc.nhi
            if m == 0
                idx += 1; v[idx] = 1.0
            else
                idx += 1; v[idx] = sqrt(2.0)
                idx += 1; v[idx] = sqrt(2.0)
            end
        end
        @test idx == nm
        v
    end

    Y    = zeros(Float64, npoin, nm)
    ctmp = zeros(Float64, nm)
    ψtmp = zeros(Float64, npoin)
    for j = 1:nm
        fill!(ctmp, 0.0); ctmp[j] = 1.0
        _forcing_stream_kernel!(ψtmp, mesh.coords, npoin, fc.nlo, fc.nhi,
                                ctmp, wunit, fc.arec, fc.brec, fc.drec)
        Y[:,j] .= ψtmp
    end

    dmax, omax = 0.0, 0.0
    for j = 1:nm, k = j:nm
        g = gsum(sum(w[ip]*Y[ip,j]*Y[ip,k] for ip = 1:npoin))/(a*a)
        j == k ? (dmax = max(dmax, abs(g - 1))) : (omax = max(omax, abs(g)))
    end
    rank == 0 && @printf("#   orthonormality: max|diag-1| = %.3e , max|offdiag| = %.3e\n", dmax, omax)

    @test dmax < 1.0e-3            # resolved, and normalised — the √2 lives here
    @test omax < 1.0e-3            # and mutually orthogonal
end

#--- 2. the round trip: the discrete curl of u_F is the prescribed ζ_F ---------
@testset "curl of u_F" begin
    sphere_forcing_step!(fc, q, mesh, metrics, sp)

    wunit = let v = zeros(Float64, fc.nmode), idx = 0
        for m = 0:fc.nhi, n = max(fc.nlo, m):fc.nhi
            if m == 0; idx += 1; v[idx] = 1.0
            else;      idx += 1; v[idx] = sqrt(2.0); idx += 1; v[idx] = sqrt(2.0); end
        end
        v
    end
    ζt = zeros(Float64, npoin)
    _forcing_stream_kernel!(ζt, mesh.coords, npoin, fc.nlo, fc.nhi, fc.c, wunit,
                            fc.arec, fc.brec, fc.drec)
    ζt .*= fc.α                                  # the amplitude the operator chose

    qf = zeros(Float64, npoin, 5)
    for ip = 1:npoin
        qf[ip,1] = φref
        qf[ip,2] = φref*fc.uF[ip,1]
        qf[ip,3] = φref*fc.uF[ip,2]
        qf[ip,4] = φref*fc.uF[ip,3]
    end
    ζg = zeros(Float64, npoin)
    sphere_relative_vorticity!(ζg, qf, mesh, metrics, sp)

    num = gsum(sum(w[ip]*(ζg[ip]-ζt[ip])^2 for ip = 1:npoin))
    den = gsum(sum(w[ip]*ζt[ip]^2          for ip = 1:npoin))
    err = sqrt(num/den)
    rank == 0 && @printf("#   ||ζ(u_F) - ζ_F|| / ||ζ_F|| = %.4e\n", err)

    @test den > 0
    @test err < 0.02          # the SEM curl of the SEM gradient of a band-limited field
end

#--- 3. u_F is on the shell ----------------------------------------------------
@testset "tangency" begin
    d = 0.0
    for ip = 1:npoin
        x, y, z = mesh.coords[1,ip], mesh.coords[2,ip], mesh.coords[3,ip]
        r = sqrt(x*x + y*y + z*z)
        d = max(d, abs(fc.uF[ip,1]*x + fc.uF[ip,2]*y + fc.uF[ip,3]*z)/r)
    end
    u = gmax(maximum(abs, fc.uF))
    d = gmax(d)
    rank == 0 && @printf("#   max|u_F·x̂| / max|u_F| = %.3e\n", d/u)
    @test u > 0
    @test d/u < 1.0e-12
end

#--- 4. the injection rate is ε₀, per unit mass --------------------------------
@testset "injection rate" begin
    @test isapprox(fc.εreal, fc.ε₀; rtol = 1.0e-10)

    # and the integrated statement: from rest, Ê(t) = ε₀ t exactly, with the
    # forcing integrated by hand so the check is on the operator and not on the
    # RK scheme
    qq  = copy(q)
    nst = 40
    for k = 1:nst
        sphere_forcing_step!(fc, qq, mesh, metrics, sp)
        for ip = 1:npoin
            qq[ip,2] += DT*qq[ip,1]*fc.uF[ip,1]
            qq[ip,3] += DT*qq[ip,1]*fc.uF[ip,2]
            qq[ip,4] += DT*qq[ip,1]*fc.uF[ip,3]
        end
    end
    E    = gsum(sum(w[ip]*0.5*(qq[ip,2]^2+qq[ip,3]^2+qq[ip,4]^2)/qq[ip,1] for ip = 1:npoin))
    Mtot = gsum(sum(w[ip]*qq[ip,1] for ip = 1:npoin))
    rank == 0 && @printf("#   Ê(t)/(ε₀ t) after %d steps = %.8f\n", nst, (E/Mtot)/(fc.ε₀*nst*DT))
    @test isapprox(E/Mtot, fc.ε₀*nst*DT; rtol = 1.0e-8)
end

#--- 5. the tendency lands where it should ------------------------------------
@testset "apply" begin
    RHS = zeros(Float64, npoin, 5)
    qe  = copy(q)
    fc.νr = 1.0e-6
    sphere_forcing_apply!(RHS, q, qe, fc, npoin)
    # at rest the friction contributes nothing and the momentum tendency is φu_F
    @test maximum(abs, view(RHS, :, 1)) == 0.0                     # continuity untouched (ν_h = 0)
    @test all(RHS[ip,2] ≈ φref*fc.uF[ip,1] for ip = 1:npoin)
    # ...and with momentum present the friction is -ν_r φu
    q2 = copy(q); q2[:,2] .= 1.0e3
    fill!(RHS, 0.0)
    sphere_forcing_apply!(RHS, q2, qe, fc, npoin)
    @test all(RHS[ip,2] ≈ φref*fc.uF[ip,1] - 1.0e-6*1.0e3 for ip = 1:npoin)
    fc.νr = 0.0
end

#--- 6. every rank drew the same field -----------------------------------------
if nparts > 1
    @testset "same field on every rank" begin
        c0 = copy(fc.c)
        MPI.Bcast!(c0, 0, comm)
        @test maximum(abs, c0 .- fc.c) == 0.0          # bit-identical, not just close
        @test gmax(Float64(fc.α)) == Float64(fc.α)     # so the amplitude agrees too
    end
end

end # testset
rank == 0 && println("# test_sphere_forcing.jl ... DONE\n")
