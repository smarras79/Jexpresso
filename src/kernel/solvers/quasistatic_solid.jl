#---------------------------------------------------------------------------------
# quasistatic_solid.jl — quasi-static, small-strain solid mechanics on
# Jexpresso's spectral-element discretisation, for an isotropic linear
# VISCOELASTIC solid.
#
# Selected by `:lquasistatic_solid => true` in user_inputs.jl. The driver then
# skips params_setup and the time integrator: there is no inertia and no
# M dq/dt = RHS. Time enters only through the material's memory and through
# time-dependent loads, and every time step is ONE linear solve.
#
# First case using it: problems/Elasticity/viscoelastic_cylinder2d (Artioli,
# Beirão da Veiga, Lovadina & Sacco, Comput. Mech. 60 (2017), Sec. 4.2).
#
# NOTHING HERE IS A NEW DISCRETISATION. Every spatial operator is Jexpresso's
# own, built on the grid, basis, metric terms and LGL quadrature of sem_setup:
#
#   stiffness       build_gradient_product_matrices (element_matrices.jl — the
#                   directional components of build_laplace_matrix, same layout
#                   and quadrature) assembled by DSS_laplace_sparse
#   tractions       compute_segment_integral! + DSS_segment_integral!
#                   (surface_integral.jl), with the tangent-based edge Jacobian
#                   and normal of build_metric_terms! — curved-edge aware, so an
#                   :exact_geometry boundary is integrated on the circle
#   body force      the lumped LGL mass matrix sem.matrix.M, as standard_linsolve!
#   stress output   _∇f! (operators.jl) + DSS_rhs! + sem.matrix.M: the SEM nodal
#                   gradient
#   field output    write_output(…, VTK()) through the case's user_uout!, with
#                   simulation.pvd kept by init_pvd_file / append_pvd_entry
#
# What this file adds is the material model and the time recursion, which act
# on assembled matrices and nodal vectors only.
#
# THE MATERIAL  (Artioli et al. Sec. 3.1, eqs. 3.2-3.7)
# -----------------------------------------------------
# Generalized Maxwell in the DEVIATOR, elastic in VOLUME:
#
#     σ = K tr(ε) 1 + s,        s(t) = ∫ 2G(t−s) ė(s) ds,
#     G(t) = G [ g₀ + Σₘ gₘ exp(−t/τₘ) ]                       (Prony series)
#
# with e = dev ε, and one internal variable per Prony term — the deviatoric
# strain qₘ in the spring of the m-th spring-dashpot branch:
#
#     q̇ₘ + qₘ/τₘ = ė,          s = 2G [ g₀ e + Σₘ gₘ qₘ ].
#
# NOTATION. τₘ is the paper's λ⁽ᵐ⁾, the RELAXATION TIME. It is written τ on
# purpose: the Lamé constant λ = K − 2G/3 is also around, and the two have
# nothing to do with each other.
#
# PLANE STRAIN. ε_zz = 0, but the material point is three-dimensional:
# e_zz = −tr(ε)/3 ≠ 0, so q_zz and s_zz are not zero either and σ_zz evolves in
# time. Only xx, yy, xy enter the weak form (ε(v)_zz = 0); σ_zz is output.
#
# TIME INTEGRATION — the exponential recursion
# --------------------------------------------
# Over a step Δt the internal-variable ODE is integrated exactly for a given
# strain history, which leaves the recursion
#
#     qₘⁿ⁺¹ = dₘ qₘⁿ + γₘ (eⁿ⁺¹ − eⁿ),        dₘ = exp(−Δt/τₘ)                (*)
#
# with γₘ chosen by :viscoelastic_integrator:
#
#   "exponential"     γ = (τ/Δt)(1 − e^{−Δt/τ})   strain LINEAR in time over the
#      (default)        step: the recursion of Zienkiewicz, Taylor & Fox (2014)
#                       and FEAP, which is what Artioli et al. ran. 2nd order.
#   "midpoint"        γ = e^{−Δt/2τ}               the hereditary integral by the
#                       midpoint rule: Simo & Hughes, Computational
#                       Inelasticity (1998), Sec. 10.3, Box 10.1. 2nd order.
#   "backward_euler"  d = γ = 1/(1 + Δt/τ)         1st order — the control that
#                       lets a Δt sweep tell the orders apart.
#
# THE HISTORY IS A NODAL FIELD. The discrete strain is linear in the nodal
# displacement and (*) is linear and starts from zero, so qₘ at every
# quadrature point is exactly dev ε(wₘ) of a nodal displacement-like field wₘ
# that obeys the same recursion:
#
#     wₘⁿ⁺¹ = dₘ wₘⁿ + γₘ (uⁿ⁺¹ − uⁿ).
#
# The stress is σ = K tr ε(u) 1 + 2G dev ε(g₀u + Σₘ gₘwₘ), and with the
# assembled operators
#
#     v·A_vol u = ∫ div v div u,         v·A_dev z = ∫ ε(v) : 2 dev ε(z)
#
# (A_dev in Voigt form is D_dev = [4/3 −2/3 0; −2/3 4/3 0; 0 0 1]) each step is
#
#     (K A_vol + G̃ A_dev) uⁿ⁺¹ = f(tⁿ⁺¹) − G A_dev Σₘ gₘ (dₘ wₘⁿ − γₘ uⁿ),
#      G̃ = G (g₀ + Σₘ gₘ γₘ)                          the algorithmic modulus.
#
# This is the quadrature-point recursion to round-off — not an approximation of
# it — and it needs no per-quadrature-point storage. (A NONLINEAR material,
# plasticity say, would: there the history is not a linear function of the
# displacement.) G̃ depends on Δt only, so for a fixed step one factorisation
# serves the whole run; factorisations are cached by G̃.
#
# In Jexpresso's flat layout u = [ux(1:npoin); uy(1:npoin)], and with Sαβ the
# matrices DSS_laplace_sparse assembles from build_gradient_product_matrices,
#
#     A_vol = [ Sxx  Sxy ;  Syx  Syy ]
#     A_dev = [ (4/3)Sxx + Syy     Syx − (2/3)Sxy ;
#               Sxy − (2/3)Syx     Sxx + (4/3)Syy ],      Syx = Sxyᵀ.
#
# THE INITIAL STEP. With :linitial_elastic_solve => true (default) the loads
# are on at t = tinit: a Δt = 0 step of (*) (d = γ = 1, G̃ = G·Σg) gives the
# instantaneous, glassy response before the time steps — a true Heaviside load,
# the problem the correspondence-principle solution describes. With false the
# run starts unloaded and the loads enter during the FIRST increment, which,
# through (*), the integrator reads as a ramp over [tinit, tinit+Δt]: what
# FEAP does when a load is switched on at step 1, and what Artioli et al.
# Fig. 5.2 shows.
#
# BOUNDARY CONDITIONS come from the case through one hook (user_bc.jl):
#
#     user_solid_bc(tag::String, x, y, t, nx, ny) -> NamedTuple
#
# called at every node of every boundary edge with that edge's gmsh tag and
# OUTWARD unit normal. Fields: `ux = value` / `uy = value` prescribe that
# displacement component (a symmetry plane or roller is ONE component);
# `traction = (tx, ty)` prescribes σ·n; an empty NamedTuple, `(;)`, is a
# traction-free surface. Constrained DOFs are eliminated — the solve is on the
# free block, which stays symmetric positive definite. The set of constrained
# DOFs is fixed at tinit; the values may depend on t.
#
# THE OTHER HOOKS
#
#     user_solid_material()           (user_flux.jl)    (; E, ν, g0, g, τ)
#     user_solid_body_force(x, y, t)  (user_source.jl)  (bx, by) per unit volume
#     user_solid_exact(x, y, t)       (user_analytic.jl, OPTIONAL)
#                                     (; ux, uy[, sxx, syy, sxy, szz])
#     initialize(...)                 qvars = SOLID_FIELDS, displacement zero
#     user_uout!(...)                 the VTK fields, from the SOLID_FIELDS and
#                                     the exact solution at each node
#
# OUTPUT
#   * VTK (:outformat => "vtk") at every :diagnostics_at_times time, through the
#     standard writer: iter_<k>.pvtu + simulation.pvd.
#   * solid_history.csv — always: every time step, the displacement at the
#     :solid_probes points (and the exact one there), and the error norms.
#
# LIMITS: 2D (NSD_2D), serial, CPU, inexact (LGL-collocated) quadrature, no
# Laguerre/semi-infinite elements. Each is checked up front, with a message.
#---------------------------------------------------------------------------------

const SOLID_INTEGRATORS = (:exponential, :midpoint, :backward_euler)

# The fields the solver hands to write_output: the case's initialize() must
# declare exactly these as qvars (its user_uout! decides what is written).
const SOLID_FIELDS = ["ux", "uy", "sxx", "syy", "sxy", "szz"]

"""
    maxwell_step_coefficients(Δt, τ, integrator) -> (d, γ)

Coefficients of one step of `qⁿ⁺¹ = d qⁿ + γ (eⁿ⁺¹ − eⁿ)` for the Maxwell
internal variable `q̇ + q/τ = ė`. `Δt == 0` is the instantaneous (glassy) jump,
`d = γ = 1`, whatever the integrator.
"""
function maxwell_step_coefficients(Δt::Real, τ::Real, integrator::Symbol)
    Δt == 0 && return (1.0, 1.0)
    x = Float64(Δt)/Float64(τ)
    if integrator === :exponential
        return (exp(-x), -expm1(-x)/x)          # (1 − e^{−x})/x, accurate at small x
    elseif integrator === :midpoint
        return (exp(-x), exp(-x/2))
    elseif integrator === :backward_euler
        return (1.0/(1.0 + x), 1.0/(1.0 + x))
    end
    error(" # ERROR quasistatic_solid.jl: unknown :viscoelastic_integrator ", integrator,
          ". Use one of ", join(string.(SOLID_INTEGRATORS), ", "), ".")
end

#
# The case's material, checked and completed with the elastic constants.
#
function solid_material_from_case()

    isdefined(@__MODULE__, :user_solid_material) ||
        error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid => true needs ",
              "user_solid_material() in the case's user_flux.jl, returning (; E, ν, g0, g, τ).")

    m  = user_solid_material()
    E  = Float64(m.E)
    ν  = Float64(m.ν)
    g0 = Float64(get(m, :g0, 1.0))
    g  = Float64[x for x in get(m, :g, Float64[])]
    τ  = Float64[x for x in get(m, :τ, Float64[])]

    E > 0               || error(" # ERROR quasistatic_solid.jl: E = $E must be positive.")
    -1.0 < ν < 0.5      || error(" # ERROR quasistatic_solid.jl: ν = $ν must lie in (-1, 1/2).")
    length(g) == length(τ) ||
        error(" # ERROR quasistatic_solid.jl: the Prony series has ", length(g),
              " weights g but ", length(τ), " relaxation times τ.")
    all(>(0), τ)        || error(" # ERROR quasistatic_solid.jl: every relaxation time τ must be positive, got $τ.")
    (g0 >= 0 && all(>=(0), g)) ||
        error(" # ERROR quasistatic_solid.jl: Prony weights must be non-negative, got g0 = $g0, g = $g.")
    g0 + sum(g; init=0.0) > 0 || error(" # ERROR quasistatic_solid.jl: g0 + Σg must be positive.")

    K = E/(3.0*(1.0 - 2.0*ν))
    G = E/(2.0*(1.0 + ν))

    return (; E, ν, K, G, g0, g, τ, nprony = length(g))
end

#
# DSS_laplace_sparse sizes its result from the largest index it stored, and it
# drops entries below eps: an off-diagonal block such as Sxy can come out
# smaller than npoin x npoin. Square it up before the blocks are stacked.
#
function _solid_npoin_square(S::SparseMatrixCSC, n::Integer)
    size(S) == (n, n) && return S
    I, J, V = findnz(S)
    return sparse(I, J, V, n, n)
end

"""
    solid_sem_operators(sem) -> (A_vol, A_dev)

The two plane-strain elasticity operators, from Jexpresso's element matrices
(build_gradient_product_matrices) and sparse DSS (DSS_laplace_sparse), in the
flat layout u = [ux; uy]: `v·A_vol u = ∫ div v div u`, `v·A_dev z = ∫ ε(v):2dev ε(z)`.
"""
function solid_sem_operators(sem)
    mesh = sem.mesh
    N    = mesh.nop
    Lxx, Lxy, Lyy = build_gradient_product_matrices(mesh.SD, sem.basis.ψ, sem.basis.dψ, sem.ω,
                                                    mesh.nelem, mesh, sem.metrics, N, N, Float64)
    n   = Int(mesh.npoin)
    Sxx = _solid_npoin_square(DSS_laplace_sparse(mesh, Lxx), n)
    Sxy = _solid_npoin_square(DSS_laplace_sparse(mesh, Lxy), n)
    Syy = _solid_npoin_square(DSS_laplace_sparse(mesh, Lyy), n)
    Syx = sparse(transpose(Sxy))

    A_vol = [Sxx  Sxy;
             Syx  Syy]
    A_dev = [(4/3)*Sxx + Syy    Syx - (2/3)*Sxy;
             Sxy - (2/3)*Syx    Sxx + (4/3)*Syy]
    return A_vol, A_dev
end

#
# Boundary loads and constraints at time t. Tractions are integrated with
# Jexpresso's segment integral (ω times the curved-edge Jacobian Jef, then the
# DSS onto the nodes); Dirichlet components are collected as dof => value, in
# the flat layout (ux of node ip is dof ip, uy is npoin + ip).
#
function solid_boundary_conditions!(f, dirichlet::Dict{Int,Float64}, t, sem, S_edge, F_edge, S_flux)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    npoin   = Int(mesh.npoin)

    empty!(dirichlet)
    fill!(S_edge, 0.0)
    fill!(S_flux, 0.0)
    for iedge = 1:mesh.nedges_bdy
        tag = mesh.bdy_edge_type[iedge]
        tag === nothing && continue
        tag = String(tag)
        fill!(F_edge, 0.0)
        for k = 1:ngl
            ip = mesh.poin_in_bdy_edge[iedge, k]
            x  = mesh.coords[1, ip]; y = mesh.coords[2, ip]
            bc = user_solid_bc(tag, x, y, t, metrics.nx[iedge, k], metrics.ny[iedge, k])
            if haskey(bc, :traction)
                F_edge[k, 1], F_edge[k, 2] = bc.traction
            end
            haskey(bc, :ux) && _solid_set_dirichlet!(dirichlet, ip,         Float64(bc.ux), tag, x, y)
            haskey(bc, :uy) && _solid_set_dirichlet!(dirichlet, npoin + ip, Float64(bc.uy), tag, x, y)
        end
        compute_segment_integral!(S_edge, F_edge, sem.ω, metrics.Jef, iedge, ngl)
    end
    DSS_segment_integral!(S_flux, S_edge, nothing, mesh.nedges_bdy, ngl, mesh.connijk,
                          mesh.poin_in_bdy_edge, mesh.bdy_edge_in_elem, 2)
    @views f[1:npoin]         .+= S_flux[:, 1]
    @views f[npoin+1:2*npoin] .+= S_flux[:, 2]
    return nothing
end

function _solid_set_dirichlet!(dirichlet, dof, value, tag, x, y)
    old = get(dirichlet, dof, nothing)
    if old !== nothing && abs(old - value) > 1.0e-12*max(1.0, abs(value))
        error(" # ERROR quasistatic_solid.jl: two boundary edges prescribe different values ",
              "(", old, " and ", value, ") on the same displacement component at (", x, ", ", y,
              "), tag \"", tag, "\".")
    end
    dirichlet[dof] = value
    return nothing
end

#
# f += ∫ φ b dΩ with the LGL-collocated (lumped) mass, as standard_linsolve!.
#
function solid_add_body_force!(f, t, sem)
    mesh  = sem.mesh
    npoin = Int(mesh.npoin)
    M     = sem.matrix.M
    for ip = 1:npoin
        bx, by = user_solid_body_force(mesh.coords[1, ip], mesh.coords[2, ip], t)
        f[ip]         += M[ip]*bx
        f[npoin + ip] += M[ip]*by
    end
    return nothing
end

"""
    solid_nodal_gradient!(grad, field, sem, f_el, ∇f_el)

The SEM nodal gradient of a nodal field: Jexpresso's element gradient _∇f!
(LGL-collocated, weighted by ωJ), DSS_rhs!, and the lumped mass matrix — i.e.
the element derivatives at a shared node averaged with the mass weights.
`grad` is npoin x 2.
"""
function solid_nodal_gradient!(grad, field, sem, f_el, ∇f_el)
    mesh    = sem.mesh
    metrics = sem.metrics
    ngl     = mesh.ngl
    for iel = 1:mesh.nelem
        for j = 1:ngl, i = 1:ngl
            f_el[i, j] = field[mesh.connijk[iel, i, j]]
        end
        _∇f!(∇f_el, f_el, ngl, sem.basis.dψ, sem.ω, metrics.Je,
             metrics.dξdx, metrics.dξdy, metrics.dηdx, metrics.dηdy,
             iel, CL(), sem.QT, mesh.SD, ContGal())
    end
    fill!(grad, 0.0)
    DSS_rhs!(grad, ∇f_el, mesh.connijk, mesh.nelem, ngl, 2, mesh.SD, ContGal())
    grad ./= sem.matrix.M
    return grad
end

function _solid_check_supported(sem, qp, inputs)
    mesh = sem.mesh
    mesh.SD isa NSD_2D ||
        error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid is implemented for 2D (plane strain) only.")
    sem.QT isa Inexact ||
        error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid needs the LGL-collocated (inexact) ",
              "quadrature; remove :lexact_integration from the deck.")
    inputs[:backend] == CPU() ||
        error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid runs on the CPU backend only.")
    MPI.Comm_size(get_mpi_comm()) == 1 ||
        error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid assembles one global sparse ",
              "matrix and is serial; run it on one MPI rank.")
    mesh.nelem_semi_inf == 0 ||
        error(" # ERROR quasistatic_solid.jl: Laguerre/semi-infinite elements are not supported.")
    for f in (:user_solid_bc, :user_solid_body_force)
        isdefined(@__MODULE__, f) ||
            error(" # ERROR quasistatic_solid.jl: :lquasistatic_solid needs `", f,
                  "` from the case (user_bc.jl / user_source.jl). See the header of ",
                  "src/kernel/solvers/quasistatic_solid.jl.")
    end
    String.(qp.qvars) == SOLID_FIELDS ||
        error(" # ERROR quasistatic_solid.jl: the case's initialize() must declare qvars = ",
              SOLID_FIELDS, ", got ", qp.qvars, ".")
    return nothing
end

"""
    quasistatic_solid_solve!(sem, qp, inputs, OUTPUT_DIR, TFloat)

March the quasi-static viscoelastic problem from `:tinit` to `:tend` in steps of
`:Δt`. See the header of this file for the model, the hooks and the output.
"""
function quasistatic_solid_solve!(sem, qp, inputs, OUTPUT_DIR::String, TFloat)

    _solid_check_supported(sem, qp, inputs)

    mesh  = sem.mesh
    ngl   = mesh.ngl
    npoin = Int(mesh.npoin)
    ndof  = 2*npoin
    M     = sem.matrix.M

    mat        = solid_material_from_case()
    integrator = Symbol(lowercase(string(get(inputs, :viscoelastic_integrator, "exponential"))))
    integrator in SOLID_INTEGRATORS ||
        error(" # ERROR quasistatic_solid.jl: :viscoelastic_integrator => \"", integrator,
              "\" is not one of ", join(string.(SOLID_INTEGRATORS), ", "), ".")
    linitial   = get(inputs, :linitial_elastic_solve, true) == true
    tinit      = Float64(inputs[:tinit])
    tend       = Float64(inputs[:tend])
    Δt         = Float64(inputs[:Δt])
    (Δt > 0 && tend > tinit) ||
        error(" # ERROR quasistatic_solid.jl: need Δt > 0 and tend > tinit (got Δt = $Δt, [$tinit, $tend]).")

    has_exact = get(inputs, :_has_analytic, false) == true && isdefined(@__MODULE__, :user_solid_exact)
    compare   = has_exact && linitial

    # The step times: uniform Δt, the last step shortened to land on tend.
    tsteps = Float64[tinit]
    while tsteps[end] < tend - 1.0e-10*max(1.0, abs(tend))
        push!(tsteps, min(tsteps[end] + Δt, tend))
    end
    abs(tsteps[end] - tend) < 1.0e-9*max(1.0, abs(tend)) && (tsteps[end] = tend)
    nsteps = length(tsteps) - 1

    # Output times: the deck's :diagnostics_at_times (a number or a list).
    _dt = get(inputs, :diagnostics_at_times, tend)
    out_times = _dt isa Number ? Float64[_dt] : Float64[x for x in _dt]
    is_output(t) = any(to -> abs(to - t) <= 1.0e-9*max(1.0, abs(t)), out_times)
    lvtk = inputs[:outformat] isa VTK

    # Probes: the mesh node nearest each requested point.
    probes = Tuple{String,Int}[]
    for p in get(inputs, :solid_probes, ())
        name, px, py = String(p[1]), Float64(p[2]), Float64(p[3])
        ipbest = argmin([(mesh.coords[1, ip] - px)^2 + (mesh.coords[2, ip] - py)^2 for ip = 1:npoin])
        push!(probes, (name, ipbest))
    end

    #-----------------------------------------------------------------------
    # Report the set-up.
    #-----------------------------------------------------------------------
    println(" # Quasi-static viscoelastic solid (plane strain) ....................")
    @printf("   SEM: nop = %d, %d elements, %d nodes, %d displacement DOFs\n", mesh.nop, mesh.nelem, npoin, ndof)
    @printf("   material: E = %g, ν = %g  →  K = %.6g, G = %.6g\n", mat.E, mat.ν, mat.K, mat.G)
    if mat.nprony == 0
        println("   no Prony terms: linear ELASTIC solid")
    else
        @printf("   Prony series G(t)/G = %g", mat.g0)
        for m = 1:mat.nprony
            @printf(" + %g exp(-t/%g)", mat.g[m], mat.τ[m])
        end
        println("      (τ = relaxation time, the paper's λ; NOT the Lamé constant)")
        abs(mat.g0 + sum(mat.g) - 1.0) > 1.0e-12 &&
            @printf("   NOTE g0 + Σg = %.6g ≠ 1: the instantaneous shear modulus is G·%.6g, not G\n",
                    mat.g0 + sum(mat.g), mat.g0 + sum(mat.g))
    end
    @printf("   time: [%g, %g], Δt = %g, %d steps, integrator \"%s\"\n", tinit, tend, Δt, nsteps, integrator)
    println(linitial ?
            "   loads ON at t = tinit: instantaneous (glassy) elastic solve first (Heaviside load)" :
            "   loads enter during the FIRST increment (FEAP convention: a ramp over [tinit, tinit+Δt])")
    for (name, ip) in probes
        @printf("   probe %-4s → node %d at (%.6g, %.6g)\n", name, ip, mesh.coords[1, ip], mesh.coords[2, ip])
    end
    has_exact && println("   exact solution: user_solid_exact (user_analytic.jl)",
                         linitial ? "" : " — NOT compared: it describes a Heaviside load at tinit")

    #-----------------------------------------------------------------------
    # Operators (Jexpresso's element matrices + DSS) and the constraints.
    #-----------------------------------------------------------------------
    A_vol, A_dev = solid_sem_operators(sem)

    nedges = Int(mesh.nedges_bdy)
    S_edge = zeros(nedges, ngl, 2)
    F_edge = zeros(ngl, 2)
    S_flux = zeros(npoin, 2)

    u = zeros(ndof)
    u[1:npoin] .= @view qp.qn[1:npoin, 1]
    u[npoin+1:ndof] .= @view qp.qn[1:npoin, 2]
    maximum(abs, u; init=0.0) == 0.0 ||
        error(" # ERROR quasistatic_solid.jl: the initial displacement from initialize() is not zero. ",
              "The quasi-static path starts from the undeformed, unstressed state: the history ",
              "that produced a prior deformation is not known.")

    f         = zeros(ndof)
    dirichlet = Dict{Int,Float64}()
    solid_boundary_conditions!(f, dirichlet, tinit, sem, S_edge, F_edge, S_flux)
    cdofs  = sort!(collect(keys(dirichlet)))
    isfree = trues(ndof); isfree[cdofs] .= false
    fdofs  = findall(isfree)
    println("   ", length(cdofs), " constrained DOFs, ", length(fdofs), " free")

    factor_cache = Dict{Float64,Any}()
    A_fc_cache   = Dict{Float64,Any}()
    function system_for(Gt)
        get!(factor_cache, Gt) do
            A = mat.K*A_vol + Gt*A_dev
            A_fc_cache[Gt] = A[fdofs, cdofs]
            try
                cholesky(Symmetric(A[fdofs, fdofs]))
            catch e
                e isa PosDefException || rethrow()
                error(" # ERROR quasistatic_solid.jl: the stiffness on the free DOFs is not positive ",
                      "definite. The Dirichlet conditions in user_solid_bc probably leave a rigid-body ",
                      "motion (translation or rotation) unconstrained.")
            end
        end
    end

    #-----------------------------------------------------------------------
    # One load/time step tⁿ → tⁿ⁺¹ = t over δt (δt = 0: the glassy jump).
    # w[m] is the nodal history field of Prony term m (see the header).
    #-----------------------------------------------------------------------
    nM = mat.nprony
    w  = [zeros(ndof) for _ = 1:nM]
    d  = zeros(nM); γ = zeros(nM)
    hist = zeros(ndof)
    function advance!(t, δt)
        Gt = mat.G*mat.g0
        fill!(hist, 0.0)
        for m = 1:nM
            d[m], γ[m] = maxwell_step_coefficients(δt, mat.τ[m], integrator)
            Gt += mat.G*mat.g[m]*γ[m]
            hist .+= mat.g[m] .* (d[m] .* w[m] .- γ[m] .* u)       # Σ gₘ (dₘ wₘⁿ − γₘ uⁿ)
        end

        fill!(f, 0.0)
        solid_boundary_conditions!(f, dirichlet, t, sem, S_edge, F_edge, S_flux)
        sort!(collect(keys(dirichlet))) == cdofs ||
            error(" # ERROR quasistatic_solid.jl: the set of Dirichlet DOFs changed at t = $t; ",
                  "user_solid_bc may change the VALUES it prescribes in time, not WHICH components.")
        solid_add_body_force!(f, t, sem)
        nM > 0 && (f .-= mat.G .* (A_dev*hist))

        F    = system_for(Gt)
        uc   = Float64[dirichlet[i] for i in cdofs]
        rhs  = f[fdofs] .- A_fc_cache[Gt]*uc
        unew = similar(u)
        unew[fdofs] = F \ rhs
        unew[cdofs] = uc

        for m = 1:nM
            w[m] .= d[m] .* w[m] .+ γ[m] .* (unew .- u)
        end
        u .= unew
        return nothing
    end

    #-----------------------------------------------------------------------
    # Nodal stresses σ = K div u 1 + 2G dev ε(g₀u + Σ gₘwₘ), by the SEM nodal
    # gradient; the six SOLID_FIELDS in the flat layout write_output reads.
    #-----------------------------------------------------------------------
    f_el  = zeros(ngl, ngl)
    ∇f_el = zeros(mesh.nelem, ngl, ngl, 2)
    gux = zeros(npoin, 2); guy = zeros(npoin, 2); gzx = zeros(npoin, 2); gzy = zeros(npoin, 2)
    z   = zeros(ndof)
    fields = zeros(npoin*length(SOLID_FIELDS))
    function compute_fields!()
        z .= mat.g0 .* u
        for m = 1:nM
            z .+= mat.g[m] .* w[m]
        end
        solid_nodal_gradient!(gux, @view(u[1:npoin]),      sem, f_el, ∇f_el)
        solid_nodal_gradient!(guy, @view(u[npoin+1:ndof]), sem, f_el, ∇f_el)
        solid_nodal_gradient!(gzx, @view(z[1:npoin]),      sem, f_el, ∇f_el)
        solid_nodal_gradient!(gzy, @view(z[npoin+1:ndof]), sem, f_el, ∇f_el)
        for ip = 1:npoin
            p   = mat.K*(gux[ip, 1] + guy[ip, 2])               # K tr ε(u)
            exx = gzx[ip, 1]; eyy = gzy[ip, 2]
            exy = 0.5*(gzx[ip, 2] + gzy[ip, 1])
            trz = exx + eyy                                     # ε_zz = 0
            fields[ip]           = u[ip]
            fields[npoin + ip]   = u[npoin + ip]
            fields[2npoin + ip]  = p + 2.0*mat.G*(exx - trz/3.0)
            fields[3npoin + ip]  = p + 2.0*mat.G*(eyy - trz/3.0)
            fields[4npoin + ip]  = 2.0*mat.G*exy
            fields[5npoin + ip]  = p + 2.0*mat.G*(-trz/3.0)
        end
        return fields
    end

    #-----------------------------------------------------------------------
    # Diagnostics: error norms (LGL mass-matrix norms), history file, VTK.
    #-----------------------------------------------------------------------
    nvar   = length(SOLID_FIELDS)
    qexact = zeros(npoin, nvar + 1)
    function exact_fields!(t)
        has_exact || return qexact
        for ip = 1:npoin
            ex = user_solid_exact(mesh.coords[1, ip], mesh.coords[2, ip], t)
            qexact[ip, 1] = ex.ux; qexact[ip, 2] = ex.uy
            if haskey(ex, :sxx)
                qexact[ip, 3] = ex.sxx; qexact[ip, 4] = ex.syy
                qexact[ip, 5] = ex.sxy; qexact[ip, 6] = ex.szz
            end
        end
        return qexact
    end
    function error_norms()
        compare || return (NaN, NaN, NaN, NaN)
        emax = 0.0; e2 = 0.0; n2 = 0.0; es2 = 0.0; ns2 = 0.0
        for ip = 1:npoin
            dx = fields[ip] - qexact[ip, 1]; dy = fields[npoin + ip] - qexact[ip, 2]
            emax = max(emax, hypot(dx, dy))
            e2  += M[ip]*(dx^2 + dy^2)
            n2  += M[ip]*(qexact[ip, 1]^2 + qexact[ip, 2]^2)
            ds = (fields[2npoin + ip] - qexact[ip, 3], fields[3npoin + ip] - qexact[ip, 4],
                  fields[4npoin + ip] - qexact[ip, 5], fields[5npoin + ip] - qexact[ip, 6])
            es2 += M[ip]*(ds[1]^2 + ds[2]^2 + 2.0*ds[3]^2 + ds[4]^2)
            ns2 += M[ip]*(qexact[ip, 3]^2 + qexact[ip, 4]^2 + 2.0*qexact[ip, 5]^2 + qexact[ip, 6]^2)
        end
        return (emax, sqrt(e2), sqrt(e2/max(n2, eps())), ns2 > 0 ? sqrt(es2/ns2) : NaN)
    end

    hist_path = joinpath(OUTPUT_DIR, "solid_history.csv")
    hio = open(hist_path, "w")
    header = String["t"]
    for (name, _) in probes
        push!(header, string(name, "_ux"), string(name, "_uy"))
        compare && push!(header, string(name, "_ux_exact"), string(name, "_uy_exact"))
    end
    compare && append!(header, ["err_u_max", "err_u_L2", "err_u_L2_rel", "err_sigma_L2_rel"])
    println(hio, join(header, ","))

    pvd_path = joinpath(OUTPUT_DIR, "simulation.pvd")
    lvtk && init_pvd_file(pvd_path)
    uaux = zeros(npoin, nvar + 1)
    iout = 0

    function record!(t, n)
        compute_fields!()
        exact_fields!(t)
        errs = error_norms()
        row = Float64[t]
        for (_, ip) in probes
            push!(row, u[ip], u[npoin + ip])
            compare && push!(row, qexact[ip, 1], qexact[ip, 2])
        end
        compare && append!(row, collect(errs))
        println(hio, join((@sprintf("%.15e", v) for v in row), ","))

        (is_output(t) || n == 0 || n == nsteps) || return nothing

        s = @sprintf("   t = %-9.5g", t)
        for (name, ip) in probes
            s *= @sprintf("  %s: ux = %.8f", name, u[ip])
            compare && (s *= @sprintf(" (exact %.8f)", qexact[ip, 1]))
        end
        compare && (s *= @sprintf("  |e|max = %.3e  e_L2,rel = %.3e", errs[1], errs[3]))
        println(s)

        if lvtk && is_output(t)
            iout += 1
            write_output(mesh.SD, fields, uaux, t, iout, mesh, nothing,
                         nothing, nothing, 0.0, 0.0, 0.0,
                         OUTPUT_DIR, inputs, qp.qvars, qp.qoutvars, inputs[:outformat];
                         nvar = nvar, qexact = qexact)
            append_pvd_entry(pvd_path, t, "iter_$(iout).pvtu")
        end
        return nothing
    end

    #-----------------------------------------------------------------------
    # March.
    #-----------------------------------------------------------------------
    try
        linitial && advance!(tinit, 0.0)          # glassy response to the load at tinit
        record!(tinit, 0)
        for n = 1:nsteps
            advance!(tsteps[n + 1], tsteps[n + 1] - tsteps[n])
            record!(tsteps[n + 1], n)
        end
    finally
        close(hio)
    end

    # Hand the final state back through the usual state array.
    for ivar = 1:nvar, ip = 1:npoin
        qp.qn[ip, ivar] = fields[(ivar - 1)*npoin + ip]
    end

    println("   history: ", hist_path)
    lvtk && println("   fields:  ", pvd_path, " (", iout, " frames)")
    println(" # Quasi-static viscoelastic solid (plane strain) .................... DONE")

    return (; u, fields, mat, tsteps, probes)
end
