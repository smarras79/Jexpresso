#---------------------------------------------------------------------------------
# Two-dimensional flux emergence in a two-temperature solar atmosphere for the
# ideal GLM-MHD equations — the Parker (undular magnetic buoyancy) instability
# of an isolated horizontal flux sheet after Shibata et al. (1989a), as set up
# in
#
#   D. Son, Y. Jang, T. Magara,
#   "A Comparative Analysis of High-resolution Shock-capturing Schemes for
#    Two-dimensional Magnetohydrodynamic Simulation of Flux Emergence in the
#    Solar Atmosphere", ApJS 277:46 (2025), Section 2.1, Eqs. (1)-(5).
#   https://doi.org/10.3847/1538-4365/adb617
#
# Nondimensional units (paper Table 1): length H₀ (photospheric pressure scale
# height), velocity C_s (photospheric adiabatic sound speed), time τ₀ = H₀/C_s,
# density ρ₀, pressure p₀ = ρ₀C_s², temperature T₀ = mC_s²/(γk_B), magnetic
# field B₀ = (ρ₀C_s²)^½ (Heaviside-Lorentz: magnetic pressure ½B²), gravity
# g₀ = C_s²/(γH₀). With H₀ = C_s = ρ₀ = 1 the perfect-gas law is
#
#     p = ρ T / γ,          g₀ = 1/γ,          γ = 1.05.
#
# Jexpresso's y is the paper's vertical coordinate z. The domain is
# [0, X_max] × [0, Z_max] = [0, 80] × [0, 35].
#
# 1. Two-temperature atmosphere (Eq. 1):
#
#     T(z) = T_ch + (T_cor - T_ch)/2 [tanh((z - z_cor)/w_tr) + 1]
#     T_ch = T₀,  T_cor = 25 T₀,  z_cor = 18 H₀,  w_tr = 0.6 H₀.
#
# 2. Magnetic flux sheet (Eqs. 2-3), parallel to x, inside the cold layer:
#
#     B(z) = [2 p(z)/β(z)]^½,      β(z) = β_* / f(z),
#     f(z) = ¼ [1 + tanh((z - z₀)/w₀)] [1 - tanh((z - z₁)/w₁)],
#     z₀ = 4 H₀,  z₁ = z₀ + D = 8 H₀,  w₀ = w₁ = 0.5 H₀,
#
#   i.e. the magnetic pressure is ½B² = p f/β_*. β_* is the plasma beta at
#   the sheet center. THE PAPER DOES NOT PRINT ITS VALUE. It is inferred here
#   from the paper's Fig. 1(b): integrating Eq. (4) for β_* = 1 gives
#   max B_x/B₀ = 0.120 at z = 4.3 H₀ and log₁₀ρ = -2.8, -7.9, -8.2 at
#   z = 8.5, 20, 35 H₀ against ≈ 0.113, -2.8, -7.9, -8.2 read off the figure
#   (β_* = 2 would give 0.091 and -3.1, -8.2, -8.5). β_* = 1 is also the
#   standard case of Shibata et al. (1989a), whose expansion-law constants
#   a₁ = 0.062, a₂ = 0.3 the paper adopts (its Eqs. 40-43).
#
# 3. Magnetostatic equilibrium (Eq. 4):
#
#     d p_total/dz = -ρ g₀,      p_total = p + ½B² = p (1 + f/β_*),
#
#   with ρ = γ p/T. Writing P = p_total this is dP/dz = -P/[T(z)(1 + f/β_*)],
#   integrated numerically from P(0) = p(0)(1 + f(0)/β_*), p(0) = ρ₀T₀/γ = 1/γ
#   (the paper normalizes ρ to its value at z = 0, Fig. 1(b)).
#
# 4. Perturbation (Eq. 5), inside the sheet over the central region
#    X_max/2 - λ/2 < x < X_max/2 + λ/2:
#
#     V_x = f(z) A C_s sin[2π (x - X_max/2)/λ],   A = 0.05,   λ = 20 H₀.
#
# The out-of-plane components w and Bz are identically zero for this problem
# (as is ψ), but they are carried along because the equation set implemented
# here is the full nine-field GLM-MHD system.
#
# 5. Discrete (well-balanced) equilibrium. The continuous state above is in
#    exact magnetostatic balance; its LGL interpolant is not, exactly. The
#    discrete residual of the vertical momentum balance of the initial state
#    is evaluated once (fe_well_balanced_table) and its negative is added as
#    a static source by user_source.jl, so that the initial state is an exact
#    equilibrium of the DISCRETE operator and no spurious settling flow, and
#    no spurious DynSGS residual, is seeded. On the shipped 80×35 mesh the
#    residual is small — at most 4·10⁻⁴ of ρg₀, at z = 0 — so the correction
#    is a refinement, not a rescue (the 0.2 C_s x-uniform settling seen in
#    early runs was resistive erosion of the sheet by an over-sensitive
#    DynSGS floor, see :dsgs_local_rel in user_inputs.jl). The residual is
#    x-independent, so it is tabulated against height, on a virtual column of
#    elements with the mesh's element height and LGL nodes — the same table
#    on every MPI rank, whatever the partition. Switch off with
#    fe_well_balanced[] = false.
#
# The divergence-cleaning speed c_h is set to the maximum wave speed
# |v| + c_f of the initial condition (the coronal sound speed √25 = 5 C_s, plus
# the perturbation) and kept constant throughout the simulation. The paper
# re-evaluates c_h every step (its Eq. 11); the initial value already
# bounds the coronal sound speed, and the loop's Alfvén speed of 4-7 C_s at
# late times (paper Fig. 5) is of the same order.
#---------------------------------------------------------------------------------

# Case parameters (paper Section 2.1). Refs so that they can be changed from
# the REPL without re-including the case.
if !@isdefined(fe_beta_star)
    const fe_beta_star = Ref{Float64}(1.0)    # plasma β at the sheet center (inferred, see header)
end
const fe_Tch   = 1.0      # T_ch/T₀
const fe_Tcor  = 25.0     # T_cor/T₀
const fe_zcor  = 18.0     # z_cor/H₀
const fe_wtr   = 0.6      # w_tr/H₀
const fe_z0    = 4.0      # z₀/H₀   lower boundary of the sheet
const fe_D     = 4.0      # D/H₀    sheet thickness
const fe_z1    = fe_z0 + fe_D
const fe_w0    = 0.5      # w₀/H₀
const fe_w1    = 0.5      # w₁/H₀
const fe_A     = 0.05     # perturbation amplitude (in C_s)
const fe_λ     = 20.0     # perturbation wavelength (in H₀)
const fe_Xmax  = 80.0     # X_max/H₀ (also the horizontal size of FE_80x35.geo)
const fe_Zmax  = 35.0     # Z_max/H₀

# Paper Eq. (1)
@inline fe_temperature(z) = fe_Tch + 0.5*(fe_Tcor - fe_Tch)*(tanh((z - fe_zcor)/fe_wtr) + 1.0)

# Paper Eq. (3)
@inline fe_fsheet(z) = 0.25*(1.0 + tanh((z - fe_z0)/fe_w0))*(1.0 - tanh((z - fe_z1)/fe_w1))

#
# Magnetostatic profile p(z), ρ(z), B_x(z) tabulated on a fine uniform grid
# 0 ≤ z ≤ zmax (paper Eq. 4, see header item 3) and linearly interpolated at
# the mesh points. Trapezoidal integration on a 1e-4 H₀ grid: the relative
# error of the table is O(1e-9), far below the SEM discretization error.
#
function fe_hydrostatic_table(zmax; nz=350001)

    γ    = γ_mhd
    β    = fe_beta_star[]
    z    = collect(range(0.0, zmax, length=nz))
    dz   = z[2] - z[1]

    # dP/dz = -P/[T(z)(1 + f/β)]  ->  ln P(z) = ln P(0) - ∫₀ᶻ dz'/[T(1 + f/β)]
    integrand = [1.0/(fe_temperature(zz)*(1.0 + fe_fsheet(zz)/β)) for zz in z]
    lnP = zeros(nz)
    for k = 2:nz
        lnP[k] = lnP[k-1] - 0.5*(integrand[k] + integrand[k-1])*dz
    end

    p0   = 1.0/γ                                  # ρ(0) = 1, T(0) = 1
    P0   = p0*(1.0 + fe_fsheet(0.0)/β)
    p    = similar(z)
    ρ    = similar(z)
    B    = similar(z)
    for k = 1:nz
        P    = P0*exp(lnP[k])
        f    = fe_fsheet(z[k])
        p[k] = P/(1.0 + f/β)
        ρ[k] = γ*p[k]/fe_temperature(z[k])
        B[k] = sqrt(2.0*p[k]*f/β)                 # ½B² = p f/β
    end
    return z, p, ρ, B
end

@inline function fe_interp(ztab, vtab, z)
    nz = length(ztab)
    dz = ztab[2] - ztab[1]
    k  = clamp(1 + floor(Int, (z - ztab[1])/dz), 1, nz - 1)
    w  = clamp((z - ztab[k])/dz, 0.0, 1.0)
    return (1.0 - w)*vtab[k] + w*vtab[k+1]
end

#
# Well-balanced correction (header item 5). Refs shared with user_source.jl.
#
if !@isdefined(fe_well_balanced)
    const fe_well_balanced = Ref{Bool}(true)
end
if !@isdefined(fe_wb_y)
    const fe_wb_y    = Ref{Vector{Float64}}(Float64[])   # sorted node heights
    const fe_wb_corr = Ref{Vector{Float64}}(Float64[])   # -residual of the ρv equation at rest
end

# Lagrange derivative matrix on the nodes ξ (barycentric form)
function fe_derivative_matrix(ξ::AbstractVector)
    n = length(ξ)
    w = ones(n)
    for j = 1:n, k = 1:n
        k == j && continue
        w[j] /= (ξ[j] - ξ[k])
    end
    D = zeros(n, n)
    for i = 1:n
        for j = 1:n
            i == j && continue
            D[i,j] = (w[j]/w[i])/(ξ[i] - ξ[j])
        end
        D[i,i] = -sum(D[i,:])
    end
    return D
end

#
# Discrete residual R(z) = -d(p + ½Bx²)/dz|_SEM - ρ g₀ of the vertical
# momentum equation at rest, on a virtual column of uniform elements of
# height Δy_el with the LGL node layout ξ (in [-1,1]) between 0 and ymax.
# At an element interface the CG assembly with lumped mass averages the two
# one-sided derivatives (equal weights for equal elements), which is what is
# tabulated there. Returns the sorted node heights and the CORRECTION -R.
#
function fe_well_balanced_table(Δy_el, ξ, ymax, ztab, ptab, ρtab, Btab)
    ngl = length(ξ)
    ne  = max(1, round(Int, ymax/Δy_el))
    D   = fe_derivative_matrix(ξ)
    ys   = Float64[]
    corr = Float64[]
    prev_top = nothing               # residual at the top node of the element below
    for ie = 1:ne
        y0 = (ie - 1)*Δy_el
        yn = [y0 + 0.5*(ξ[j] + 1.0)*Δy_el for j = 1:ngl]
        P  = [fe_interp(ztab, ptab, y) + 0.5*fe_interp(ztab, Btab, y)^2 for y in yn]
        ρ  = [fe_interp(ztab, ρtab, y) for y in yn]
        dPdy = (2.0/Δy_el)*(D*P)
        R = -dPdy .- ρ .* g_mhd
        for j = 1:ngl
            if j == 1 && prev_top !== nothing
                corr[end] = -0.5*(prev_top + R[1])      # shared interface node
            else
                push!(ys, yn[j]); push!(corr, -R[j])
            end
        end
        prev_top = R[ngl]
    end
    return ys, corr
end

# Correction at height y (nearest tabulated node; the table holds every
# node height of the column, so the lookup is exact on the mesh).
@inline function fe_wb_lookup(y)
    ys = fe_wb_y[]
    n  = length(ys)
    n == 0 && return 0.0
    k = searchsortedfirst(ys, y)
    k > n && return fe_wb_corr[][n]
    k == 1 && return fe_wb_corr[][1]
    return (y - ys[k-1] < ys[k] - y) ? fe_wb_corr[][k-1] : fe_wb_corr[][k]
end

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 2D ideal GLM-MHD (flux emergence, Son et al. 2025) ........... "
    end

    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: the length of qvars defines neqs. Slot 4 MUST carry the total
    # energy ρE (not ρw) because Jexpresso's shared 2D kernels assume the
    # energy lives in slot 4 — see the header of user_flux.jl.
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T", "vA", "β"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------

    if (inputs[:backend] != CPU())
        error(" problems/MHD/fluxEmergenceSon2025: only the CPU backend is supported for now.")
    end
    if (inputs[:SOL_VARS_TYPE] != TOTAL())
        error(" problems/MHD/fluxEmergenceSon2025: only SOL_VARS_TYPE = TOTAL() is supported.")
    end

    γ   = γ_mhd
    γm1 = γ - 1.0

    #
    # Global extent of the mesh (the rank-local extrema are not the domain
    # under MPI). Used for the sanity check below, the hydrostatic table and
    # the absorbing layer.
    #
    npoin = mesh.npoin
    xmin_g = MPI.Allreduce(minimum(view(mesh.x, 1:npoin)), MPI.MIN, comm)
    xmax_g = MPI.Allreduce(maximum(view(mesh.x, 1:npoin)), MPI.MAX, comm)
    ymin_g = MPI.Allreduce(minimum(view(mesh.y, 1:npoin)), MPI.MIN, comm)
    ymax_g = MPI.Allreduce(maximum(view(mesh.y, 1:npoin)), MPI.MAX, comm)
    if rank == 0 && (abs(xmin_g) > 1e-8 || abs(xmax_g - fe_Xmax) > 1e-6 || abs(ymin_g) > 1e-8 || abs(ymax_g - fe_Zmax) > 1e-6)
        @warn " problems/MHD/fluxEmergenceSon2025: the mesh spans [$(xmin_g), $(xmax_g)] × [$(ymin_g), $(ymax_g)] but the paper's domain is [0, $(fe_Xmax)] × [0, $(fe_Zmax)]. The perturbation is centered on x = $(0.5*fe_Xmax)."
    end
    sponge_ztop_mhd[] = ymax_g

    ztab, ptab, ρtab, Btab = fe_hydrostatic_table(max(ymax_g, fe_Zmax))

    ch_local = 0.0
    xc       = 0.5*fe_Xmax

    for ip = 1:npoin

        x, z = mesh.x[ip], mesh.y[ip]

        p  = fe_interp(ztab, ptab, z)
        ρ  = fe_interp(ztab, ρtab, z)
        Bx = fe_interp(ztab, Btab, z)
        f  = fe_fsheet(z)

        # Perturbation (paper Eq. 5)
        u = 0.0
        if abs(x - xc) < 0.5*fe_λ
            u = f*fe_A*sin(2.0*π*(x - xc)/fe_λ)
        end
        v  = 0.0
        w  = 0.0
        By = 0.0
        Bz = 0.0
        ψ  = 0.0

        ρE = p/γm1 + 0.5*ρ*(u*u + v*v + w*w) + 0.5*(Bx*Bx + By*By + Bz*Bz) + 0.5*ψ*ψ

        q.qn[ip,1] = ρ
        q.qn[ip,2] = ρ*u
        q.qn[ip,3] = ρ*v
        q.qn[ip,4] = ρE
        q.qn[ip,5] = ρ*w
        q.qn[ip,6] = Bx
        q.qn[ip,7] = By
        q.qn[ip,8] = Bz
        q.qn[ip,9] = ψ
        q.qn[ip,end] = p

        # Reference (magnetostatic) state: the absorbing layer of
        # user_source.jl relaxes towards it, and it is the background of
        # the perturbation output. The perturbation velocity is NOT part of
        # it.
        ρEe = p/γm1 + 0.5*(Bx*Bx)
        q.qe[ip,1] = ρ
        q.qe[ip,2] = 0.0
        q.qe[ip,3] = 0.0
        q.qe[ip,4] = ρEe
        q.qe[ip,5] = 0.0
        q.qe[ip,6] = Bx
        q.qe[ip,7] = 0.0
        q.qe[ip,8] = 0.0
        q.qe[ip,9] = 0.0
        q.qe[ip,end] = p

        # Local maximum wave speed |v| + c_f, with the fast magnetosonic
        # speed maximized over propagation directions: c_f ≤ sqrt(a² + b²),
        # a² = γp/ρ (sound), b² = |B|²/ρ (Alfvén).
        vmag = sqrt(u*u + v*v + w*w)
        cf   = sqrt(γ*p/ρ + (Bx*Bx + By*By + Bz*Bz)/ρ)
        ch_local = max(ch_local, vmag + cf)
    end

    #
    # GLM divergence-cleaning speed: max wave speed of the IC over the whole
    # (global) domain, constant in time.
    #
    c_h_mhd[] = MPI.Allreduce(ch_local, MPI.MAX, comm)

    #
    # Smallest nodal spacing Δh of the mesh, for the GLM damping parameter
    # α_p = Δh c_h/c_p² of user_source.jl (paper Eq. 12). Measured along the
    # element edges: consecutive LGL nodes of the first row/column of every
    # element.
    #
    dh_local = Inf
    for iel = 1:mesh.nelem
        for i = 1:mesh.ngl-1
            ia = mesh.connijk[iel, i,   1]
            ib = mesh.connijk[iel, i+1, 1]
            ja = mesh.connijk[iel, 1,   i]
            jb = mesh.connijk[iel, 1,   i+1]
            dh_local = min(dh_local,
                           sqrt((mesh.x[ib] - mesh.x[ia])^2 + (mesh.y[ib] - mesh.y[ia])^2),
                           sqrt((mesh.x[jb] - mesh.x[ja])^2 + (mesh.y[jb] - mesh.y[ja])^2))
        end
    end
    glm_dh_mhd[] = MPI.Allreduce(dh_local, MPI.MIN, comm)

    #
    # Well-balanced correction table (header item 5): element height and LGL
    # layout read off the first local element's vertical edge, then the
    # residual of the discrete vertical balance tabulated on a virtual column.
    #
    if fe_well_balanced[]
        ngl  = mesh.ngl
        # the local index that runs along y
        y11 = mesh.y[mesh.connijk[1,1,1]]
        jdir_is_j = abs(mesh.y[mesh.connijk[1,1,2]] - y11) > abs(mesh.y[mesh.connijk[1,2,1]] - y11)
        ycol = [jdir_is_j ? mesh.y[mesh.connijk[1,1,j]] : mesh.y[mesh.connijk[1,j,1]] for j = 1:ngl]
        ylo, yhi = extrema(ycol)
        Δy_el = yhi - ylo
        ξ     = sort((2.0 .* (ycol .- ylo) ./ Δy_el) .- 1.0)
        ys, corr = fe_well_balanced_table(Δy_el, ξ, ymax_g, ztab, ptab, ρtab, Btab)
        fe_wb_y[]    = ys
        fe_wb_corr[] = corr
        if rank == 0
            imax = argmax(abs.(corr))
            @info " Well-balanced correction: $(length(ys)) node heights, max |residual| of the vertical balance = $(abs(corr[imax])) at z = $(ys[imax]) (ρg₀ there = $(fe_interp(ztab, ρtab, ys[imax])*g_mhd))"
        end
    else
        fe_wb_y[] = Float64[]; fe_wb_corr[] = Float64[]
    end

    if rank == 0
        @info " Flux sheet plasma beta β_* = $(fe_beta_star[]) (inferred from the paper's Fig. 1(b), see initialize.jl)"
        @info " Initial max B_x/B₀ = $(maximum(Btab)) at z/H₀ = $(ztab[argmax(Btab)]); ρ/ρ₀ at the top = $(ρtab[end]), p/p₀ at the top = $(ptab[end])"
        @info " GLM divergence-cleaning speed c_h = $(c_h_mhd[]); smallest nodal spacing Δh = $(glm_dh_mhd[]); ψ damping rate α_p c_h/Δh = $(glm_alpha_p_mhd[]*c_h_mhd[]/glm_dh_mhd[]) per τ₀"
        @info " Absorbing layer: z_s = $(sponge_zs_mhd[]) to Z_max = $(sponge_ztop_mhd[]), σ_max = $(sponge_sigma_mhd[]) per τ₀"
        @info " Initialize fields for 2D ideal GLM-MHD (flux emergence, Son et al. 2025) ........... DONE"
    end

    return q
end
