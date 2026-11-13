#---------------------------------------------------------------------------------
# Initial state of the 2.5D flux-emergence tests of
#
#   J. E. Leake, M. G. Linton & S. K. Antiochos, ApJ 722:550 (2010), Sec. 2.2,
#   Erratum ApJ 741:125 (2011),
#
# in the code units of user_flux.jl (L0 = 150 km; ρ, p, T in photospheric units;
# B in units of (4π p_ph)^½ = 1030.6 G; p = ρT; g = LK10_G). Mesh: LK10.geo,
# x in [-300, 300], y in [-20, 600] (-45 to 45 Mm, -3 to 90 Mm).
#
# 1. Atmosphere (Eq. 9), hydrostatic, p(0) = ρ(0) = T(0) = 1:
#      T = 1 − g y/(m + 1)                                   y ≤ 0
#      T = 1 + (T_cor − 1)/2 [tanh((y − y_tr)/w_tr) + 1]     y > 0
#    T_cor = 8.6e5/5700 = 150.9, y_tr = 25 (3.75 Mm), w_tr = 5 (0.75 Mm).
#    Below the surface p = T^(m+1) (closed form); above, ln p = −g∫dy/T on a
#    0.002 L0 table. The paper sets m = 1/(γ−1) and calls the polytrope
#    "marginally unstable to convection"; its own Fig. 3(a) has p ∝ T^2.3 in
#    the convection zone, i.e. m ≈ 1.3, which IS slightly unstable. The deck's
#    :lk10_mcz (default 1.3) selects it; 1.5 is the exactly adiabatic polytrope.
#
# 2. Flux tube (Eqs. 10-16), axis along z at (0, y_t), y_t = −12 (−1.8 Mm):
#      B_z = B0 exp(−r²/a²),  B0 = 7800 G,  a = 2 (0.3 Mm),
#      B_θ = Θ(r) r B_z,  Θ1 = c,  Θ2 = c(1 − e^(−r²/c1²)),  Θ3 = c e^(−r²/c2²),
#      c = 1/a, c1 = 1 (0.15 Mm), c2 = 10 (1.5 Mm)
#    (tube 1: constant twist; tube 2: increasing twist, "low tension"; tube 3:
#    decreasing twist). The in-plane field circulates clockwise, (B_x, B_y) =
#    Θ B_z (y − y_t, −x), so that B_x > 0 at the top of the tube: parallel to the
#    inner arcade of the quadrupole and antiparallel to the dipole, as the paper
#    requires (its Fig. 5 arrows). Radial balance (Eq. 16, Heaviside-Lorentz):
#      p1(r) = −B²(r)/2 + ∫_r^∞ B_θ²/r' dr'     (closed form per tube),
#    and the tube is in thermal equilibrium: p = p_bg + p1, ρ = p/T_bg(y).
#
# 3. Coronal field, potential (lk10_field):
#    "quadrupole" (Eq. 17, after Karpen et al. 1996): A_z on the source line
#      y_s = −50 (−7.5 Mm) is q_d (x/L0)² [1 − (x/x_a)²] for |x| ≤ x_a = 140
#      (21 Mm), 0 outside, and A_z is harmonic above it with A_z → 0 at
#      infinity: the Poisson integral, evaluated by Gauss-Legendre quadrature over
#      the source (lk10_quadrupole). q_d = 5.85e7 G cm ("98 G"); 11.7e7 and 3.9e7
#      are the paper's strong and weak quadrupoles. x enters in units of
#      L0 = 150 km: this is the field of the paper's Fig. 3(b) (README).
#    "dipole" (Eq. 18): A_z = d (y − y_d)/r1³, y_d = −100 (−15 Mm), d set by
#      |B(0, 0)| = 104, 52 or 26 G (the paper's Fig. 3(b) and its Sec. 3.4-3.5).
#    "none": field-free corona.
#    B = ∇×(A_z e_z) = (∂A_z/∂y, −∂A_z/∂x).
#
# The reference state qe is the atmosphere plus the coronal field (an exact
# magnetohydrostatic equilibrium, made a discrete one by :lwell_balanced); the
# tube is the departure from it.
#---------------------------------------------------------------------------------

# Case parameters, filled from the deck at initialize (Refs: readable from the REPL)
if !@isdefined(lk10_tube)
    const lk10_tube   = Ref{Int}(2)
    const lk10_field  = Ref{String}("quadrupole")
    const lk10_qd     = Ref{Float64}(5.85e7)     # G cm
    const lk10_Bdip   = Ref{Float64}(52.0)       # G, |B(0,0)| of the dipole
    const lk10_mcz    = Ref{Float64}(1.3)
    const lk10_GL     = Ref{Tuple{Vector{Float64},Vector{Float64}}}((Float64[], Float64[]))
end
if !@isdefined(LK10_YT)
    const LK10_YT  = -12.0     # tube axis height (−1.8 Mm)
    const LK10_A   = 2.0       # tube width a (0.3 Mm)
    const LK10_B0G = 7800.0    # G, axial field at the axis
    const LK10_C1  = 1.0       # c1 (0.15 Mm)
    const LK10_C2  = 10.0      # c2 (1.5 Mm)
    const LK10_XA  = 140.0     # quadrupole half-extent on the source line (21 Mm)
    const LK10_YS  = -50.0     # quadrupole source line (−7.5 Mm)
    const LK10_YD  = -100.0    # dipole source (−15 Mm)
end

#---------------------------------------------------------------------------------
# 1. Atmosphere
#---------------------------------------------------------------------------------
@inline function lk10_temperature(y)
    y <= 0.0 && return 1.0 - LK10_G*y/(lk10_mcz[] + 1.0)
    return 1.0 + 0.5*(LK10_TCOR - 1.0)*(tanh((y - LK10_YTR)/LK10_WTR) + 1.0)
end

# ln p above the surface on a uniform table, trapezoidal rule on 0.002 L0
# (the relative error of p is below 1e-7 at the top of the domain)
function lk10_lnp_table(ymax; dy=0.002)
    n  = max(2, ceil(Int, ymax/dy) + 1)
    ys = collect(range(0.0, ymax, length=n))
    h  = ys[2] - ys[1]
    lnp = zeros(n)
    f0  = 1.0/lk10_temperature(1e-300)          # T(0⁺)
    for k = 2:n
        f1 = 1.0/lk10_temperature(ys[k])
        lnp[k] = lnp[k-1] - LK10_G*0.5*(f0 + f1)*h
        f0 = f1
    end
    return ys, lnp
end

@inline function lk10_pressure(y, ys, lnp)
    if y <= 0.0
        return lk10_temperature(y)^(lk10_mcz[] + 1.0)
    end
    h = ys[2] - ys[1]
    s = y/h
    k = clamp(floor(Int, s), 0, length(ys) - 2)
    w = s - k
    return exp((1.0 - w)*lnp[k+1] + w*lnp[k+2])
end

#---------------------------------------------------------------------------------
# 2. Flux tube: (B_x, B_y, B_z, p1) at (x, y)
#---------------------------------------------------------------------------------
@inline function lk10_twist(r2, tube)
    c = 1.0/LK10_A
    tube == 1 && return c
    tube == 2 && return c*(1.0 - exp(-r2/LK10_C1^2))
    return c*exp(-r2/LK10_C2^2)
end

function lk10_tube_field(x, y, tube)
    B0 = LK10_B0G/LK10_BU
    a2 = LK10_A^2
    c  = 1.0/LK10_A
    Y  = y - LK10_YT
    r2 = x*x + Y*Y
    Bz = B0*exp(-r2/a2)
    Θ  = lk10_twist(r2, tube)
    Bx = Θ*Bz*Y
    By = -Θ*Bz*x
    # I(r) = ∫_r^∞ B_θ²/r' dr' = ½ c² B0² ∫_{r²}^∞ (Θ/c)² e^{−2u/a²} du
    k0 = 2.0/a2
    if tube == 1
        I = 0.5*c*c*B0*B0*exp(-k0*r2)/k0
    elseif tube == 2
        k1 = k0 + 1.0/LK10_C1^2
        k2 = k0 + 2.0/LK10_C1^2
        I  = 0.5*c*c*B0*B0*(exp(-k0*r2)/k0 - 2.0*exp(-k1*r2)/k1 + exp(-k2*r2)/k2)
    else
        k3 = k0 + 2.0/LK10_C2^2
        I  = 0.5*c*c*B0*B0*exp(-k3*r2)/k3
    end
    p1 = -0.5*(Bx*Bx + By*By + Bz*Bz) + I
    return Bx, By, Bz, p1
end

# Vector potential of the tube's in-plane field, A_z(r) = −∫_r^∞ Θ B_z r' dr'
# (so that ∂A_z/∂y = B_x, −∂A_z/∂x = B_y and A_z → 0 far from the tube)
function lk10_tube_potential(x, y, tube)
    B0 = LK10_B0G/LK10_BU
    a2 = LK10_A^2
    c  = 1.0/LK10_A
    r2 = x*x + (y - LK10_YT)^2
    tube == 1 && return -0.5*c*B0*a2*exp(-r2/a2)
    if tube == 2
        k = 1.0/a2 + 1.0/LK10_C1^2
        return -0.5*c*B0*(a2*exp(-r2/a2) - exp(-k*r2)/k)
    end
    k3 = 1.0/a2 + 1.0/LK10_C2^2
    return -0.5*c*B0*exp(-k3*r2)/k3
end

#---------------------------------------------------------------------------------
# 3. Coronal field: (B_x, B_y, A_z) at (x, y)
#---------------------------------------------------------------------------------
# Quadrupole: A_z(x, y) = (1/π) ∫ A_s(ξ) h/((x − ξ)² + h²) dξ, h = y − y_s, over
# the source |ξ| ≤ x_a with A_s = q̂ ξ²(1 − ξ²/x_a²), q̂ = q_d/(B0 L0). The kernel
# is smooth (h ≥ 30 in the domain against a 280-wide source): 600 Gauss-Legendre
# points give the field to round-off.
function lk10_quadrupole(x, y)
    ξs, ws = lk10_GL[]
    q̂  = lk10_qd[]/(LK10_BU*LK10_L0)
    h  = y - LK10_YS
    h2 = h*h
    A = 0.0; dAdh = 0.0; dAdx = 0.0
    @inbounds for k in eachindex(ξs)
        ξ  = ξs[k]
        As = ws[k]*q̂*ξ*ξ*(1.0 - (ξ/LK10_XA)^2)
        s  = x - ξ
        d  = s*s + h2
        A    += As*h/d
        dAdh += As*(s*s - h2)/(d*d)
        dAdx += As*(-2.0*h*s)/(d*d)
    end
    return dAdh/π, -dAdx/π, A/π
end

# Dipole: A_z = d̂ (y − y_d)/r1³ with d̂ = B_dip |y_d|³/(2 B0), so that
# B_x(0, 0) = −B_dip (antiparallel to the top of the tube).
function lk10_dipole(x, y)
    d̂  = (lk10_Bdip[]/LK10_BU)*abs(LK10_YD)^3/2.0
    s  = y - LK10_YD
    r2 = x*x + s*s
    r5 = r2*r2*sqrt(r2)
    return d̂*(r2 - 3.0*s*s)/r5, 3.0*d̂*s*x/r5, d̂*s/(r2*sqrt(r2))
end

function lk10_coronal_field(x, y)
    f = lk10_field[]
    f == "quadrupole" && return lk10_quadrupole(x, y)
    f == "dipole"     && return lk10_dipole(x, y)
    return 0.0, 0.0, 0.0
end

#---------------------------------------------------------------------------------
# Gauss-Legendre nodes on [-x_a, x_a] (Golub-Welsch: eigenvalues of the Jacobi
# matrix), set once
#---------------------------------------------------------------------------------
function lk10_gauss_legendre(n)
    β = [k/sqrt(4.0*k*k - 1.0) for k = 1:n-1]
    J = SymTridiagonal(zeros(n), β)
    ev = eigen(J)
    return ev.values, 2.0 .* ev.vectors[1, :].^2
end

#---------------------------------------------------------------------------------
# Discretely divergence-free in-plane field. The CG operator takes the
# divergence of the nodal field element by element and averages it at shared
# nodes with the lumped-mass weights ω_i ω_j J (strong form + DSS, rhs.jl). If
# the nodal field is itself the lumped-mass average of the element curls of a
# continuous nodal A_z, that discrete divergence is zero to round-off (checked
# on uniform and stretched quadrilaterals), while the nodal interpolant of the
# analytic, divergence-free field has a discrete divergence of 0.2-0.7 B0/L0 for
# the 0.3 Mm tube on 2.5 L0 elements. That error is what the GLM ψ then carries
# away, and in the low-β corona a field error δB is a Lorentz force of B δB/ρ:
# with the interpolated tube the coronal velocities reached 0.3 v0 in 3 s, a
# growing disturbance around the quadrupole null. This is the CG counterpart
# of the constrained-transport initialization of the paper's Lare code, whose
# face-centred B is the discrete curl of A. Shared nodes on partition
# boundaries are summed over the ranks through the global node numbering
# (mesh.ip2gip).
#---------------------------------------------------------------------------------
function lk10_lgl(ngl)
    N = ngl - 1
    x = [-cos(π*k/N) for k = 0:N]
    P = zeros(ngl, ngl)
    for _ = 1:100
        P[:,1] .= 1.0; P[:,2] .= x
        for k = 3:ngl
            P[:,k] .= ((2k - 3) .* x .* P[:,k-1] .- (k - 2) .* P[:,k-2]) ./ (k - 1)
        end
        xo = copy(x)
        x .= xo .- (x .* P[:,ngl] .- P[:,ngl-1]) ./ (ngl .* P[:,ngl])
        maximum(abs.(x .- xo)) < 1e-15 && break
    end
    w = 2.0 ./ (N*(N + 1) .* P[:,ngl].^2)
    D = zeros(ngl, ngl)                 # D[i,k] = l_k'(x_i)
    c = [prod(x[i] - x[k] for k = 1:ngl if k != i) for i = 1:ngl]
    for i = 1:ngl, k = 1:ngl
        i != k && (D[i,k] = c[i]/(c[k]*(x[i] - x[k])))
    end
    for i = 1:ngl
        D[i,i] = -sum(D[i,k] for k = 1:ngl if k != i)
    end
    return x, w, D
end

function lk10_discrete_curl(As::Vector{Vector{Float64}}, mesh, comm)
    ngl = mesh.ngl; nelem = mesh.nelem; npoin = mesh.npoin
    _, ω, D = lk10_lgl(ngl)
    nA  = length(As)
    nprocs = MPI.Comm_size(comm)
    # global node numbering (serial runs may not fill it: then the local one)
    gid = (length(mesh.ip2gip) >= npoin && all(>(0), @view(mesh.ip2gip[1:npoin]))) ?
          Int.(mesh.ip2gip[1:npoin]) :
          (nprocs == 1 ? collect(1:npoin) :
           error(" problems/MHD/fluxEmergenceLeake2010: mesh.ip2gip is not set; cannot assemble the discrete curl across ranks."))
    gn  = MPI.Allreduce(maximum(gid), MPI.MAX, comm)
    acc = zeros(2nA + 1, gn)            # rows: (Bx, By) per potential, then the weight
    X = zeros(ngl, ngl); Y = zeros(ngl, ngl); A = zeros(ngl, ngl)
    for e = 1:nelem
        for j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[e,i,j]
            X[i,j] = mesh.coords[1,ip]; Y[i,j] = mesh.coords[2,ip]
        end
        for j = 1:ngl, i = 1:ngl
            xξ = 0.0; xη = 0.0; yξ = 0.0; yη = 0.0
            for k = 1:ngl
                xξ += D[i,k]*X[k,j]; xη += D[j,k]*X[i,k]
                yξ += D[i,k]*Y[k,j]; yη += D[j,k]*Y[i,k]
            end
            J  = xξ*yη - xη*yξ
            ξx = yη/J; ξy = -xη/J; ηx = -yξ/J; ηy = xξ/J
            w  = ω[i]*ω[j]*abs(J)
            gip = gid[mesh.connijk[e,i,j]]
            for n = 1:nA
                Av = As[n]
                Aξ = 0.0; Aη = 0.0
                for k = 1:ngl
                    Aξ += D[i,k]*Av[mesh.connijk[e,k,j]]
                    Aη += D[j,k]*Av[mesh.connijk[e,i,k]]
                end
                acc[2n-1, gip] += w*(Aξ*ξy + Aη*ηy)        #  B_x =  ∂A/∂y
                acc[2n,   gip] -= w*(Aξ*ξx + Aη*ηx)        #  B_y = −∂A/∂x
            end
            acc[2nA+1, gip] += w
        end
    end
    MPI.Comm_size(comm) > 1 && MPI.Allreduce!(acc, +, comm)
    B = [zeros(2, npoin) for _ = 1:nA]
    for ip = 1:npoin
        gip = gid[ip]
        wsum = acc[2nA+1, gip]
        for n = 1:nA
            B[n][1,ip] = acc[2n-1, gip]/wsum
            B[n][2,ip] = acc[2n,   gip]/wsum
        end
    end
    return B
end

# Density on the background relation p = ρ T₀(ρ) for a given pressure (bisection
# in ln ρ over the T₀(ρ) table; ρT₀(ρ) grows with ρ along the atmosphere)
function lk10_rho_on_background(p)
    lo = lk10_lnrho_lo[]
    hi = lo + (length(lk10_Ttab[]) - 1)*lk10_dlnrho[]
    for _ = 1:100
        mid = 0.5*(lo + hi)
        ρ   = exp(mid)
        ρ*lk10_T0_of_rho(ρ) < p ? (lo = mid) : (hi = mid)
    end
    return exp(0.5*(lo + hi))
end

#---------------------------------------------------------------------------------
# initialize
#---------------------------------------------------------------------------------
function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank == 0 && println(" Initialize fields for 2.5D GLM-MHD flux emergence (Leake, Linton & Antiochos 2010) ........... ")

    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T", "β", "log10rho_cgs", "Bz_G"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    inputs[:backend] != CPU()         && error(" problems/MHD/fluxEmergenceLeake2010: only the CPU backend is supported.")
    inputs[:SOL_VARS_TYPE] != TOTAL() && error(" problems/MHD/fluxEmergenceLeake2010: only SOL_VARS_TYPE = TOTAL() is supported.")

    # case parameters
    lk10_tube[]  = Int(get(inputs, :lk10_tube, 2))
    lk10_field[] = lowercase(String(get(inputs, :lk10_field, "quadrupole")))
    lk10_qd[]    = Float64(get(inputs, :lk10_qd, 5.85e7))
    lk10_Bdip[]  = Float64(get(inputs, :lk10_Bdip, 52.0))
    lk10_mcz[]   = Float64(get(inputs, :lk10_mcz, 1.3))
    lk10_tube[] in (0, 1, 2, 3) || error(" problems/MHD/fluxEmergenceLeake2010: :lk10_tube = $(lk10_tube[]), expected 1, 2, 3 (or 0: no tube).")
    lk10_field[] in ("none", "quadrupole", "dipole") ||
        error(" problems/MHD/fluxEmergenceLeake2010: :lk10_field = \"$(lk10_field[])\", expected \"none\", \"quadrupole\" or \"dipole\".")
    lk10_relax[]        = Bool(get(inputs, :lk10_relax, true))
    lk10_relax_ymin[]   = Float64(get(inputs, :lk10_relax_ymin, -Inf))
    lk10_glm_alpha[]    = Float64(get(inputs, :lk10_glm_alpha, 1.0))
    lk10_damp_xd[]      = Float64(get(inputs, :lk10_damp_xd, 250.0))
    lk10_damp_yd[]      = Float64(get(inputs, :lk10_damp_yd, 500.0))
    lk10_damp_rate[]    = Float64(get(inputs, :lk10_damp_rate, 1.0))
    if lk10_field[] == "quadrupole"
        ξ, w = lk10_gauss_legendre(600)
        lk10_GL[] = (LK10_XA .* ξ, LK10_XA .* w)
    end

    npoin  = mesh.npoin
    xmin_g = MPI.Allreduce(minimum(@view(mesh.coords[1,1:npoin])), MPI.MIN, comm)
    xmax_g = MPI.Allreduce(maximum(@view(mesh.coords[1,1:npoin])), MPI.MAX, comm)
    ymin_g = MPI.Allreduce(minimum(@view(mesh.coords[2,1:npoin])), MPI.MIN, comm)
    ymax_g = MPI.Allreduce(maximum(@view(mesh.coords[2,1:npoin])), MPI.MAX, comm)
    lk10_xmax[] = max(abs(xmin_g), abs(xmax_g))
    lk10_ymax[] = ymax_g

    # atmosphere table above the surface
    ys_tab, lnp_tab = lk10_lnp_table(max(ymax_g, 1.0) + 1.0)

    # T₀(ρ) table for the energy relaxation: the background temperature against
    # ln ρ, uniform in ln ρ between the densities of the top and of the bottom
    let nfine = 400001
        yf  = collect(range(min(ymin_g, 0.0) - 1.0, ymax_g + 1.0, length=nfine))
        lnρ = [log(lk10_pressure(y, ys_tab, lnp_tab)/lk10_temperature(y)) for y in yf]   # decreasing
        Tf  = [lk10_temperature(y) for y in yf]
        lo, hi = lnρ[end], lnρ[1]
        ntab = 20001
        dl   = (hi - lo)/(ntab - 1)
        tab  = zeros(ntab)
        rl = reverse(lnρ); rT = reverse(Tf)                 # increasing in ln ρ
        for k = 1:ntab
            l = lo + (k - 1)*dl
            j = clamp(searchsortedlast(rl, l), 1, nfine - 1)
            w = (l - rl[j])/(rl[j+1] - rl[j])
            tab[k] = (1.0 - w)*rT[j] + w*rT[j+1]
        end
        lk10_lnrho_lo[] = lo; lk10_dlnrho[] = dl; lk10_Ttab[] = tab
    end

    γm1   = γ_mhd - 1.0
    tube  = lk10_tube[]
    tube_eq = lowercase(String(get(inputs, :lk10_tube_eq, "thermal")))
    tube_eq in ("thermal", "relaxed") ||
        error(" problems/MHD/fluxEmergenceLeake2010: :lk10_tube_eq = \"$(tube_eq)\", expected \"thermal\" or \"relaxed\".")
    relaxed_tube = tube_eq == "relaxed"

    # In-plane fields as discrete curls of the nodal vector potentials (see
    # lk10_discrete_curl): the coronal field alone (reference) and tube + corona
    # (initial state). B_z is the nodal value of the analytic axial field.
    Acor = zeros(npoin); Atot = zeros(npoin)
    for ip = 1:npoin
        x, y = mesh.coords[1,ip], mesh.coords[2,ip]
        Acor[ip] = lk10_coronal_field(x, y)[3]
        Atot[ip] = Acor[ip] + (tube == 0 ? 0.0 : lk10_tube_potential(x, y, tube))
    end
    Bcor, Btot = lk10_discrete_curl([Acor, Atot], mesh, comm)

    ch_local = 0.0; cmax_x = 0.0; cmax_y = 0.0
    pmin_local = Inf
    for ip = 1:npoin
        x, y = mesh.coords[1,ip], mesh.coords[2,ip]
        T   = lk10_temperature(y)
        pbg = lk10_pressure(y, ys_tab, lnp_tab)
        ρbg = pbg/T
        Bxc, Byc = Bcor[1,ip], Bcor[2,ip]

        # reference state: atmosphere + coronal (potential) field
        q.qe[ip,1] = ρbg
        q.qe[ip,2] = 0.0; q.qe[ip,3] = 0.0; q.qe[ip,5] = 0.0
        q.qe[ip,4] = pbg/γm1 + 0.5*(Bxc*Bxc + Byc*Byc)
        q.qe[ip,6] = Bxc; q.qe[ip,7] = Byc; q.qe[ip,8] = 0.0
        q.qe[ip,9] = 0.0
        q.qe[ip,end] = pbg

        # initial state: + the buoyant tube (p1 from the analytic tube field,
        # Eq. 16), in thermal equilibrium with its surroundings (the paper), or,
        # with :lk10_tube_eq => "relaxed", on the background relation ρT₀(ρ) = p
        # that the energy relaxation drives it to within ~1 s
        _, _, Bzt, p1 = tube == 0 ? (0.0, 0.0, 0.0, 0.0) : lk10_tube_field(x, y, tube)
        p  = pbg + p1
        ρ  = (relaxed_tube && abs(p1) > 1e-12*pbg) ? lk10_rho_on_background(p) : p/T
        Bx = Btot[1,ip]; By = Btot[2,ip]; Bz = Bzt
        q.qn[ip,1] = ρ
        q.qn[ip,2] = 0.0; q.qn[ip,3] = 0.0; q.qn[ip,5] = 0.0
        q.qn[ip,4] = p/γm1 + 0.5*(Bx*Bx + By*By + Bz*Bz)
        q.qn[ip,6] = Bx; q.qn[ip,7] = By; q.qn[ip,8] = Bz
        q.qn[ip,9] = 0.0
        q.qn[ip,end] = p

        pmin_local = min(pmin_local, p)
        cf = sqrt((γ_mhd*p + Bx*Bx + By*By + Bz*Bz)/ρ)      # fast-speed bound at rest
        if cf > ch_local
            ch_local = cf; cmax_x = x; cmax_y = y
        end
    end

    # GLM cleaning speed: the fastest wave of the initial state (or the deck's
    # :lk10_ch), held constant
    cmax = MPI.Allreduce(ch_local, MPI.MAX, comm)
    chin = Float64(get(inputs, :lk10_ch, 0.0))
    c_h_mhd[] = chin > 0.0 ? chin : cmax
    pmin = MPI.Allreduce(pmin_local, MPI.MIN, comm)
    pmin <= 0.0 && error(" problems/MHD/fluxEmergenceLeake2010: non-positive initial pressure ($(pmin)).")

    # smallest nodal spacing (Dedner damping rate α_p c_h/Δh)
    dh_local = Inf
    for iel = 1:mesh.nelem
        for i = 1:mesh.ngl-1
            ia = mesh.connijk[iel,i,1];  ib = mesh.connijk[iel,i+1,1]
            ja = mesh.connijk[iel,1,i];  jb = mesh.connijk[iel,1,i+1]
            dh_local = min(dh_local,
                           hypot(mesh.coords[1,ib] - mesh.coords[1,ia], mesh.coords[2,ib] - mesh.coords[2,ia]),
                           hypot(mesh.coords[1,jb] - mesh.coords[1,ja], mesh.coords[2,jb] - mesh.coords[2,ja]))
        end
    end
    lk10_glm_dh[] = MPI.Allreduce(dh_local, MPI.MIN, comm)

    if rank == 0
        Bt  = LK10_B0G/LK10_BU
        _, _, _, p1c = tube == 0 ? (0.0, 0.0, 0.0, 0.0) : lk10_tube_field(0.0, LK10_YT, tube)
        pbt = lk10_pressure(LK10_YT, ys_tab, lnp_tab)
        println(" Units: L0 = 150 km, v0 = $(round(LK10_V0/1e5, digits=4)) km/s, t0 = $(round(LK10_T0, digits=4)) s, ",
                "ρ_ph = $(round(LK10_RHO0, sigdigits=5)) g/cm³, p_ph = $(round(LK10_PPH, sigdigits=5)) dyn/cm², ",
                "B0 = $(round(LK10_BU, digits=2)) G, g = $(round(LK10_G, digits=5))")
        println(" Atmosphere: m_cz = $(lk10_mcz[]) (p ∝ T^$(lk10_mcz[] + 1) below y = 0), T_cor = $(round(LK10_TCOR, digits=3)), ",
                "ρ(-20) = $(round(lk10_pressure(-20.0, ys_tab, lnp_tab)/lk10_temperature(-20.0), sigdigits=5)), ",
                "ρ(25) = $(round(lk10_pressure(25.0, ys_tab, lnp_tab)/lk10_temperature(25.0), sigdigits=5)), ",
                "ρ($(round(ymax_g))) = $(round(lk10_pressure(ymax_g, ys_tab, lnp_tab)/lk10_temperature(ymax_g), sigdigits=5))")
        if tube > 0
            ρbt = pbt/lk10_temperature(LK10_YT)
            ρt  = relaxed_tube ? lk10_rho_on_background(pbt + p1c) : (pbt + p1c)/lk10_temperature(LK10_YT)
            println(" Tube $(tube) ($(tube_eq)): B0 = $(round(Bt, digits=4)) (7800 G), a = $(LK10_A), y_t = $(LK10_YT); at the axis ",
                    "p_bg = $(round(pbt, sigdigits=5)), p1 = $(round(p1c, sigdigits=5)), β = $(round((pbt + p1c)/(0.5*Bt^2), sigdigits=4)), ",
                    "density deficit $(round((1 - ρt/ρbt)*100, digits=1)) %")
        else
            println(" No tube (:lk10_tube => 0): the background alone, which must stay at rest.")
        end
        if lk10_field[] == "quadrupole"
            # null on x = 0 and the surface field
            yy = range(0.0, 60.0, length=6001)
            Bx0 = [lk10_quadrupole(0.0, y)[1] for y in yy]
            inull = findfirst(k -> sign(Bx0[k]) != sign(Bx0[k+1]), 1:length(yy)-1)
            xs = range(0.0, 200.0, length=2001)
            Bs = [hypot(lk10_quadrupole(x, 0.0)[1:2]...) for x in xs]
            println(" Quadrupole: q_d = $(lk10_qd[]) G cm; |B|(x_a/2, 0) = $(round(Bs[argmin(abs.(xs .- 70.0))]*LK10_BU, digits=2)) G, ",
                    "max |B| on y = 0: $(round(maximum(Bs)*LK10_BU, digits=2)) G at x = $(round(xs[argmax(Bs)]*0.15, digits=2)) Mm, ",
                    "null at y = $(inull === nothing ? "none" : round(yy[inull]*0.15, digits=3)) Mm")
        elseif lk10_field[] == "dipole"
            println(" Dipole: |B(0,0)| = $(lk10_Bdip[]) G, d = $(round(lk10_Bdip[]*(15e8)^3/2, sigdigits=4)) G cm³, y_d = -15 Mm")
        else
            println(" Field-free corona.")
        end
        println(" Fastest wave of the initial state: c_f = $(round(cmax, digits=3)) at (x, y) = ($(round(cmax_x, digits=1)), $(round(cmax_y, digits=1)))",
                " [$(round(cmax*LK10_V0/1e5, digits=1)) km/s]")
        println(" GLM: c_h = $(round(c_h_mhd[], digits=3)), smallest nodal spacing Δh = $(round(lk10_glm_dh[], digits=4)), ",
                "ψ damping rate α_p c_h/Δh = $(round(lk10_glm_alpha[]*c_h_mhd[]/lk10_glm_dh[], digits=2)) per t0; ",
                "c_h Δt/Δh = $(round(c_h_mhd[]*inputs[:Δt]/lk10_glm_dh[], digits=3))")
        println(" Energy relaxation: $(lk10_relax[] ? "on" : "OFF"), τ = $(round(lk10_tau0[], digits=4)) ρ^(-$(lk10_tau_exp[])) t0; ",
                "above y = $(lk10_relax_ymin[]); damping layers |x| > $(lk10_damp_xd[]), y > $(lk10_damp_yd[]), rate $(lk10_damp_rate[])/t0 at the boundary")
        println(" Dissipation: ν = $(round(LK10_NU, sigdigits=4)) (kinematic), anomalous η₀ = $(round(LK10_ETA0, sigdigits=4)), ",
                "j_crit = $(round(LK10_JCRIT, sigdigits=4)) (code units)")
        println(" Initialize fields for 2.5D GLM-MHD flux emergence (Leake et al. 2010) ........... DONE")
        flush(stdout)
    end

    return q
end
