# =============================================================================
# test_adaptive_angular_manufactured  (fresh rewrite, NOT compiled or run)
#
# What it checks, with one regrid BETWEEN two solves (production pattern):
#   1. Solve on the base angular mesh               -> error e0, fingerprint
#   2. Regrid (lobe-based adaptive refinement), invalidate every angular cache,
#      rebuild infrastructure
#   3. Solve on the new mesh, warm-started from the remapped solution -> e1, matvecs_warm
#   4. Solve on the new mesh from zero                                 -> e1c, matvecs_cold
#   5. PASS/FAIL lines:
#        - warm and cold converge to the same solution
#        - warm start needs fewer matvecs than cold start
#        - e1 <= e0 (refining near the sun direction reduces the error)
#
# Manufactured solution (narrow angular peak at the sun direction, so adaptivity matters):
#   I(x,Ω) = I0 + A sin(k·x) f(Ω),  f = 1 + b Ωx + c Ωy + d Ωz + B exp(-(1-Ω·Ωs)/w)
# The scattering integral S(Ω) = ∫ P(Ω,Ω') I(x,Ω') dΩ' is evaluated by a FINE reference
# quadrature, so the error includes the angular quadrature error of the coupling matrix C
# that adaptivity is supposed to reduce:
#   S = I0 s0(Ω) + A sin(k·x) F(Ω),  s0 = ∫P,  F = ∫P f
# No assumption is made about the normalisation of phase_function_hg.
#
# ASSUMPTIONS to verify (marked "ASSUMPTION" below):
#   A1  (now checked against your source) adaptive_angular_regrid! returns (infra_new, snapshots, info);
#   A2  infra_new has `metrics` and `n_dirs` like `setup` (as test_adaptive_angular_checks assumes)
#   A3  g_vary mode: build_phase_coupling_matrices_elementwise uses the plain mean of g_field over the
#       element's nodes. In that mode the error is only indicative (element-wise g makes the source
#       discontinuous at element faces); only the warm/cold consistency check is meaningful.
#
# Run one `g` per job (0.2, 0.6, 0.85), plus one job with g_vary=true.
# =============================================================================

# ---------- fine reference quadrature on the sphere (midpoint rule) ----------
struct RefQuad
    θ::Vector{Float64}; ϕ::Vector{Float64}; w::Vector{Float64}
    Ωx::Vector{Float64}; Ωy::Vector{Float64}; Ωz::Vector{Float64}
end
function make_ref_quad(nθ::Int, nϕ::Int)
    θ = Float64[]; ϕ = Float64[]; w = Float64[]
    for i = 1:nθ, j = 1:nϕ
        t = (i - 0.5) * π / nθ; p = (j - 0.5) * 2π / nϕ
        push!(θ, t); push!(ϕ, p); push!(w, sin(t) * (π / nθ) * (2π / nϕ))
    end
    RefQuad(θ, ϕ, w, sin.(θ) .* cos.(ϕ), sin.(θ) .* sin.(ϕ), cos.(θ))
end

lobe_f(Ωx, Ωy, Ωz, Ωs, b, c, d, B, w) =
    1 + b*Ωx + c*Ωy + d*Ωz + B*exp(-(1 - (Ωx*Ωs[1] + Ωy*Ωs[2] + Ωz*Ωs[3])) / w)

# ---------- remap of the warm start: nearest direction within the same spatial element ----------
snapshot_dirs(mesh, I_all) = Dict(iel => (copy(mesh.extra_mesh[iel].extra_coords[1,:]),
                                         copy(mesh.extra_mesh[iel].extra_coords[2,:]),
                                         [copy(I_all[(iel,dd)]) for dd = 1:mesh.extra_mesh[iel].extra_npoin])
                                  for iel = 1:mesh.nelem)

function remap_nearest(mesh_new, snap)
    out = Dict{Tuple{Int,Int},Vector{Float64}}()
    for iel = 1:mesh_new.nelem
        em = mesh_new.extra_mesh[iel]; θo, ϕo, vals = snap[iel]
        Ωo = [(sin(θo[m])*cos(ϕo[m]), sin(θo[m])*sin(ϕo[m]), cos(θo[m])) for m = 1:length(θo)]
        for dd = 1:em.extra_npoin
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]
            Ω = (sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ))
            best = 1; bd = Inf
            for m = 1:length(Ωo)
                dist = 1 - (Ω[1]*Ωo[m][1] + Ω[2]*Ωo[m][2] + Ω[3]*Ωo[m][3])
                dist < bd && (bd = dist; best = m)
            end
            out[(iel,dd)] = copy(vals[best])
        end
    end
    return out
end

# ---------- everything that keys on the angular mesh must go ----------
function invalidate_all_angular_caches!()
    clear_face_plans!()
    for nm in (:ASYNC_LAYOUTS, :ASYNC_DEPS, :GS_CYCLE_CACHE)         # no-ops if clear_face_plans! already did it
        try
            isdefined(@__MODULE__, nm) && empty!(getfield(@__MODULE__, nm))
        catch
        end
    end
    return nothing
end

# ---------- manufactured case; mesh-angular-independent, so it is valid before and after a regrid ----------
function lobe_manufactured_case(mesh, rq::RefQuad; g_const::Float64, g_vary::Bool, g_hi::Float64,
                                I0, A, b, c, d, B, lobe_w, Ωs, κ_val, σ_val)
    Lx = mesh.xmax-mesh.xmin; Ly = mesh.ymax-mesh.ymin; Lz = mesh.zmax-mesh.zmin
    kx = 2π/Lx; ky = 2π/Ly; kz = 2π/Lz
    g_field = (x,y,z,ip) -> g_vary ? g_const + (g_hi - g_const) * (z - mesh.zmin) / Lz : g_const

    iel_of = zeros(Int, length(mesh.x))
    for iel = 1:mesh.nelem, k=1:mesh.ngl, j=1:mesh.ngl, i=1:mesh.ngl
        ip = mesh.connijk[iel,i,j,k]; iel_of[ip] == 0 && (iel_of[ip] = iel)
    end
    gmean = zeros(mesh.nelem)                                               # ASSUMPTION A3
    for iel = 1:mesh.nelem
        s = 0.0
        for k=1:mesh.ngl, j=1:mesh.ngl, i=1:mesh.ngl
            ip = mesh.connijk[iel,i,j,k]; s += g_field(mesh.x[ip], mesh.y[ip], mesh.z[ip], ip)
        end
        gmean[iel] = s / mesh.ngl^3
    end

    fq = [lobe_f(rq.Ωx[m], rq.Ωy[m], rq.Ωz[m], Ωs, b, c, d, B, lobe_w) for m = 1:length(rq.w)]
    rows = Dict{NTuple{3,Float64},NTuple{2,Float64}}()
    function row(θ::Float64, ϕ::Float64, g::Float64)
        get!(rows, (θ, ϕ, round(g; digits=10))) do
            s0 = 0.0; F = 0.0
            for m = 1:length(rq.w)
                Pw = phase_function_hg(θ, rq.θ[m], ϕ, rq.ϕ[m], g) * rq.w[m]
                s0 += Pw; F += Pw * fq[m]
            end
            (s0, F)
        end
    end

    I_manuf(x,y,z,θ,ϕ) = begin
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        I0 + A*sin(kx*x + ky*y + kz*z) * lobe_f(Ωx, Ωy, Ωz, Ωs, b, c, d, B, lobe_w)
    end
    source_fn = (x,y,z,ip,θ,ϕ) -> begin
        Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
        f = lobe_f(Ωx, Ωy, Ωz, Ωs, b, c, d, B, lobe_w)
        φ = kx*x + ky*y + kz*z
        Ω_dot_grad = A*cos(φ) * (Ωx*kx + Ωy*ky + Ωz*kz) * f
        s0, F = row(Float64(θ), Float64(ϕ), gmean[iel_of[ip]])
        Ω_dot_grad + κ_val*I_manuf(x,y,z,θ,ϕ) - σ_val * (I0*s0 + A*sin(φ)*F)
    end
    case = RadiativeCase(
        kappa_fn = (x,y,z,ip) -> κ_val,
        sigma_fn = (x,y,z,ip) -> σ_val,
        phase_fn = phase_function_hg, g_const = g_const, g_field = g_field,
        bc_fn = (x,y,z,ip,nx,ny,nz,θ,ϕ) -> I_manuf(x,y,z,θ,ϕ),
        source_fn = source_fn)
    return case, I_manuf
end

function max_error_vs_manufactured(mesh, I_all, I_manuf, comm)
    e = 0.0
    for iel = 1:mesh.nelem
        em = mesh.extra_mesh[iel]
        for dd = 1:em.extra_npoin
            θ = em.extra_coords[1,dd]; ϕ = em.extra_coords[2,dd]; v = I_all[(iel,dd)]
            for k=1:mesh.ngl, j=1:mesh.ngl, i=1:mesh.ngl
                gip = mesh.connijk[iel,i,j,k]
                e = max(e, abs(v[lidx(i,j,k,mesh.ngl)] - I_manuf(mesh.x[gip], mesh.y[gip], mesh.z[gip], θ, ϕ)))
            end
        end
    end
    return MPI.Allreduce(e, MPI.MAX, comm)
end

function max_diff(mesh, Ia, Ib, comm)
    e = 0.0
    for iel = 1:mesh.nelem, dd = 1:mesh.extra_mesh[iel].extra_npoin
        e = max(e, maximum(abs.(Ia[(iel,dd)] .- Ib[(iel,dd)])))
    end
    return MPI.Allreduce(e, MPI.MAX, comm)
end

# ---------- regrid wrapper: the only place that touches the adaptivity API (A1, A2) ----------
function regrid_between_solves(mesh, sem, setup, comm, case, builder, Ωs;
                               tol_quad, max_level, strength_tol, lobe_k)
    κ, σ, _, _ = build_case_fields(mesh, case)
    npoin = length(mesh.x)
    g_node = [case.g_field(mesh.x[ip], mesh.y[ip], mesh.z[ip], ip) for ip = 1:npoin]
    F_node = ones(npoin)                                                   # unit sun flux
    infra_new, _, info = adaptive_angular_regrid!(mesh, sem, setup, comm, case.phase_fn, g_node, σ, F_node,
                                                  Ωs_to_angles(Ωs)..., builder;
                                                  tol_quad=tol_quad, max_level=max_level,
                                                  strength_tol=strength_tol, lobe_k=lobe_k)
    return infra_new, info
end
Ωs_to_angles(Ωs) = (acos(Ωs[3]), mod(atan(Ωs[2], Ωs[1]), 2π))

# =============================================================================
function test_adaptive_angular_manufactured(inputs, nparts, distribute, sem, partitioned_model;
        g::Float64=0.6, g_vary::Bool=false, g_hi::Float64=0.85,
        θsun::Float64=0.6, ϕsun::Float64=0.8,
        B::Float64=3.0, lobe_w::Float64=0.05,
        I0=4.2, A=1.0, b=0.1, c=0.1, d=0.1,
        tau_kappa::Float64=0.3, tau_sigma::Float64=1.0,        # optical depths across the box
        tol_quad::Float64=1e-3, max_level::Int=3, strength_tol::Float64=1e-3, lobe_k::Float64=2.0,
        solver::Symbol=:gmres, tol::Float64=1e-10, max_iter::Int=300,
        gmres_restart::Int=40, gmres_max_restarts::Int=10,
        ref_nθ::Int=240, ref_nϕ::Int=480,
        opts::SweepOptions=SweepOptions())
    comm = get_mpi_comm(); rank = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
 
    setup = setup_unified_ncf_test(inputs, nparts, distribute, sem, partitioned_model;
                                    spatial_target_gels=Int[], angular_targets=Tuple{Int,Int}[])
    mesh = setup.mesh; metrics = setup.metrics          # spatial metrics do not change under angular regrid
    L = max(mesh.xmax-mesh.xmin, mesh.ymax-mesh.ymin, mesh.zmax-mesh.zmin)
    κ_val = tau_kappa / L; σ_val = tau_sigma / L
    Ωs = (sin(θsun)*cos(ϕsun), sin(θsun)*sin(ϕsun), cos(θsun))
    rq = make_ref_quad(ref_nθ, ref_nϕ)
    builder = make_angular_mesh_builder(inputs, inputs[:backend])
    case, I_manuf = lobe_manufactured_case(mesh, rq; g_const=g, g_vary=g_vary, g_hi=g_hi,
                                           I0=I0, A=A, b=b, c=c, d=d, B=B, lobe_w=lobe_w, Ωs=Ωs,
                                           κ_val=κ_val, σ_val=σ_val)
    kw = (solver=solver, tol=tol, max_iter=max_iter, gmres_restart=gmres_restart,
          gmres_max_restarts=gmres_max_restarts, opts=opts)
    nmv(rep) = length(rep[:matvec_times])
    rank == 0 && println("adaptive manufactured: g=$g vary=$g_vary sun=($θsun,$ϕsun) lobe(B=$B,w=$lobe_w) " *
                         "tau_kappa=$tau_kappa tau_sigma=$tau_sigma solver=$solver ranks=$np")
 
    # ---- 1. base mesh ----
    rep0 = Dict{Symbol,Any}()
    I0_all, _, _ = run_unified_case_solve(setup, metrics, comm, case; report=rep0, kw...)
    e0 = max_error_vs_manufactured(mesh, I0_all, I_manuf, comm)
    nd0 = MPI.Allreduce(sum(setup.n_dirs), +, comm)
    rank == 0 && println("[base]   pairs=$nd0 matvecs=$(nmv(rep0)) max error=$e0")
    fp0 = solution_fingerprint(mesh, I0_all, comm); rank == 0 && foreach(println, fp0)
 
    # ---- 2. regrid between solves ----
    snap = snapshot_dirs(mesh, I0_all)
    invalidate_all_angular_caches!()
    setup1, rinfo = regrid_between_solves(mesh, sem, setup, comm, case, builder, Ωs;
                                   tol_quad=tol_quad, max_level=max_level,
                                   strength_tol=strength_tol, lobe_k=lobe_k)
    invalidate_all_angular_caches!()                                       # the rebuild may have repopulated stale tables
    mesh1 = setup1.mesh
    nd1 = MPI.Allreduce(sum(setup1.n_dirs), +, comm)
    rank == 0 && println("[regrid] pairs $nd0 -> $nd1 (x$(round(nd1/nd0; digits=2))), levels min/max = $(minimum(values(rinfo.levels)))/$(maximum(values(rinfo.levels)))")
    I_init = remap_nearest(mesh1, snap)
 
    # ---- 3. warm-started solve on the new mesh ----
    rep1 = Dict{Symbol,Any}()
    I1_all, _, _ = run_unified_case_solve(setup1, metrics, comm, case;
                                          I_cg_all_init=I_init, report=rep1, kw...)
    e1 = max_error_vs_manufactured(mesh1, I1_all, I_manuf, comm)
    rank == 0 && println("[warm]   matvecs=$(nmv(rep1)) max error=$e1")
 
    # ---- 4. cold-started solve on the new mesh ----
    rep2 = Dict{Symbol,Any}()
    I2_all, _, _ = run_unified_case_solve(setup1, metrics, comm, case; report=rep2, kw...)
    e2 = max_error_vs_manufactured(mesh1, I2_all, I_manuf, comm)
    dwc = max_diff(mesh1, I1_all, I2_all, comm)
    rank == 0 && println("[cold]   matvecs=$(nmv(rep2)) max error=$e2")
    fp1 = solution_fingerprint(mesh1, I1_all, comm); rank == 0 && foreach(println, fp1)
 
    # ---- 5. verdicts ----
    ok_same = dwc < 1e-6
    ok_fewer = nmv(rep1) < nmv(rep2)
    ok_err = e1 <= e0
    if rank == 0
        println("CHECK warm vs cold max diff = $dwc  -> $(ok_same ? "PASS" : "FAIL") (want < 1e-6)")
        println("CHECK warm matvecs $(nmv(rep1)) < cold matvecs $(nmv(rep2))  -> $(ok_fewer ? "PASS" : "FAIL")")
        println("CHECK error after regrid $e1 <= before $e0  -> $(ok_err ? "PASS" : (g_vary ? "n/a (g_vary: indicative only)" : "FAIL"))")
        println("ADAPTIVE-MANUFACTURED g=$g vary=$g_vary ranks=$np pairs $nd0->$nd1 e0=$e0 e1=$e1 warm=$(nmv(rep1)) cold=$(nmv(rep2))")
    end
    return (e0=e0, e1=e1, warm_cold_diff=dwc, matvecs_warm=nmv(rep1), matvecs_cold=nmv(rep2), pairs=(nd0, nd1))
end
