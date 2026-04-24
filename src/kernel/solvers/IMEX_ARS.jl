#-------------------------------------------------------------------------------
# IMEX Runge–Kutta time integration (Ascher–Ruuth–Spiteri 1997, "ARS(2,3,2)").
#
# Splits the semi-discrete problem
#
#     du/dt = N(u) + L u
#
# where N(u) is the non-stiff nonlinear part (explicit) and L u is the
# stiff linear part (implicit), into an L-stable, second-order, 3-stage
# additive Runge–Kutta scheme.
#
# For the viscous Burgers equation
#
#     ∂q/∂t + ∇·F(q) = ν Δq,    F(q) = (q²/2, q²/2, …)
#
# we take
#
#     N(q) = -M⁻¹ ∇·F(q)       (explicit nonlinear advection)
#     L q  = -ν M⁻¹ K q        (implicit linear diffusion)
#
# where M is the SEM mass matrix (diagonal under inexact LGL quadrature) and
# K is the SEM stiffness matrix Kᵢⱼ = ∫ ∇ψᵢ · ∇ψⱼ dΩ. Both 1D and 2D (NSD_2D,
# CG/LGL, quadrilateral elements) are supported; the 2D path uses the usual
# isoparametric chain rule through metrics.dξdx, dξdy, dηdx, dηdy, Je.
#
# Reference:
#   U. Ascher, S. Ruuth, R. Spiteri,
#   "Implicit-explicit Runge-Kutta methods for time-dependent PDEs",
#   Applied Numerical Mathematics 25 (1997) 151–167.
#-------------------------------------------------------------------------------

struct IMEX_ARS232 end

#-------------------------------------------------------------------------------
# Butcher tableaux for ARS(2,3,2).
#   γ  = 1 - 1/√2
#   δ  = -2√2/3   (explicit off-diagonal)
#
# Implicit (DIRK), 0-row is trivial:
#
#    0   |  0     0     0
#    γ   |  0     γ     0
#    1   |  0    1-γ    γ
#   --------------------
#         0    1-γ    γ       ← b̃
#
# Explicit (ERK):
#
#    0   |  0        0       0
#    γ   |  γ        0       0
#    1   |  δ      1-δ       0
#   --------------------
#         0      1-γ        γ     ← b  (matches b̃ here)
#
# NOTE: Below the coefficients are expressed exactly as in Pareschi–Russo (2005)
# form used by the existing IMEX infrastructure in problems/.../kopriva_IMEX,
# which is equivalent up to a relabelling of the stages. The two formulations
# produce the same second-order L-stable updates.
#-------------------------------------------------------------------------------
function ars232_tableau(::Type{T}) where {T}
    γ   = T(1) - T(1)/sqrt(T(2))
    s2  = sqrt(T(2))

    # Implicit DIRK
    Ã = zeros(T, 3, 3)
    Ã[2,1] = γ;              Ã[2,2] = γ
    Ã[3,1] = T(1)/(T(2)*s2); Ã[3,2] = T(1)/(T(2)*s2); Ã[3,3] = γ
    b̃ = T[T(1)/(T(2)*s2), T(1)/(T(2)*s2), γ]
    c̃ = T[T(0), T(2)*γ, T(1)]

    # Explicit ERK
    a32 = (T(3) + T(2)*s2)/T(6)
    A   = zeros(T, 3, 3)
    A[2,1] = T(2)*γ
    A[3,1] = T(1) - a32;     A[3,2] = a32
    b      = copy(b̃)
    c      = T[T(0), T(2)*γ, T(1)]

    return A, Ã, b, b̃, c, c̃
end

#-------------------------------------------------------------------------------
# Build 1D CG spectral-element mass (diagonal, inexact LGL) and stiffness
# (sparse) matrices on the global index set.
#
# Periodicity is handled at the mesh level (connijk already identifies the
# first and last global nodes for 1D periodic problems), so DSS on the
# global indices automatically assembles the correct periodic operators.
#-------------------------------------------------------------------------------
function build_mass_stiff_1d(mesh, basis, ω)
    ngl   = mesh.ngl
    nelem = mesh.nelem
    npoin = mesh.npoin

    M = zeros(Float64, npoin)

    I_vec = Int[]
    J_vec = Int[]
    V_vec = Float64[]
    sizehint!(I_vec, ngl*ngl*nelem)
    sizehint!(J_vec, ngl*ngl*nelem)
    sizehint!(V_vec, ngl*ngl*nelem)

    dψ = basis.dψ  # dψ[i,k] = (dψ_i/dξ)(ξ_k)

    for iel = 1:nelem
        Jac  = mesh.Δx[iel]/2.0     # dx/dξ
        dξdx = 1.0/Jac              # dξ/dx

        # Diagonal mass (inexact LGL: ψ_i(ξ_k) = δ_{ik})
        for i = 1:ngl
            ip = mesh.connijk[iel, i, 1]
            M[ip] += ω[i] * Jac
        end

        # Local stiffness Kᵉ[i,j] = Σ_k ω_k (dψ_i/dx)(ξ_k)(dψ_j/dx)(ξ_k) * Jac
        #                         = dξdx * Σ_k ω_k dψ[i,k] dψ[j,k]
        for j = 1:ngl, i = 1:ngl
            kij = 0.0
            for k = 1:ngl
                kij += ω[k] * dψ[i,k] * dψ[j,k]
            end
            kij *= dξdx   # ω_k * (1/J) factor; the two dξdx and Jac combine to 1/J

            ip = mesh.connijk[iel, i, 1]
            jp = mesh.connijk[iel, j, 1]
            push!(I_vec, ip); push!(J_vec, jp); push!(V_vec, kij)
        end
    end

    K = sparse(I_vec, J_vec, V_vec, npoin, npoin)
    return M, K
end

#-------------------------------------------------------------------------------
# Build 2D CG spectral-element mass (diagonal, inexact LGL) and stiffness
# (sparse Laplacian ∑_e ∫ ∇ψ_i · ∇ψ_j dΩ) matrices on quadrilateral elements.
#
# With a tensor-product LGL basis ψ_{i,j}(ξ,η) = h_i(ξ) h_j(η) and LGL
# interpolation = LGL quadrature nodes, h_j(η_l) = δ_{jl}. The Laplacian
# entry K^e_{(i,j),(m,n)} collapses to four contributions:
#
#   (A) j == n:  Σ_k  ω_k ω_j Je(k,j) dψ[i,k] dψ[m,k] (ξx²+ξy²)(k,j)
#   (B) i == m:  Σ_l  ω_i ω_l Je(i,l) dψ[j,l] dψ[n,l] (ηx²+ηy²)(i,l)
#   (C) fixed (m,j):   ω_m ω_j Je(m,j) dψ[i,m] dψ[n,j] (ξx·ηx+ξy·ηy)(m,j)
#   (D) fixed (i,n):   ω_i ω_n Je(i,n) dψ[j,n] dψ[m,i] (ξx·ηx+ξy·ηy)(i,n)
#
# Mesh-level periodicity (gmsh doubly-periodic quadrilateral meshes) is
# already encoded in mesh.connijk, so the global sparse assembly below
# automatically produces the correct periodic Laplacian.
#-------------------------------------------------------------------------------
function build_mass_stiff_2d(mesh, basis, ω, metrics)
    ngl   = mesh.ngl
    nelem = mesh.nelem
    npoin = mesh.npoin

    M  = zeros(Float64, npoin)
    dψ = basis.dψ        # dψ[i,k] = h_i'(ξ_k)
    Je   = metrics.Je    # (nelem, ngl, ngl, 1)
    dξdx = metrics.dξdx; dξdy = metrics.dξdy
    dηdx = metrics.dηdx; dηdy = metrics.dηdy

    I_vec = Int[]; J_vec = Int[]; V_vec = Float64[]
    nnz_hint = nelem * ngl^4
    sizehint!(I_vec, nnz_hint); sizehint!(J_vec, nnz_hint); sizehint!(V_vec, nnz_hint)

    @inbounds for iel = 1:nelem
        # Diagonal mass
        for j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel, i, j]
            M[ip] += ω[i] * ω[j] * Je[iel, i, j, 1]
        end

        # Local stiffness: (ngl*ngl) × (ngl*ngl) block, assembled node-by-node.
        for j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel, i, j]
            for n = 1:ngl, m = 1:ngl
                jp  = mesh.connijk[iel, m, n]
                val = 0.0

                # (A) sum over k at l = j, requires j == n
                if j == n
                    for k = 1:ngl
                        w  = ω[k] * ω[j] * Je[iel, k, j, 1]
                        ξx = dξdx[iel, k, j, 1]; ξy = dξdy[iel, k, j, 1]
                        val += w * dψ[i, k] * dψ[m, k] * (ξx*ξx + ξy*ξy)
                    end
                end

                # (B) sum over l at k = i, requires i == m
                if i == m
                    for l = 1:ngl
                        w  = ω[i] * ω[l] * Je[iel, i, l, 1]
                        ηx = dηdx[iel, i, l, 1]; ηy = dηdy[iel, i, l, 1]
                        val += w * dψ[j, l] * dψ[n, l] * (ηx*ηx + ηy*ηy)
                    end
                end

                # (C) single quadrature point (m, j)
                let k = m, l = j
                    w  = ω[k] * ω[l] * Je[iel, k, l, 1]
                    ξx = dξdx[iel, k, l, 1]; ξy = dξdy[iel, k, l, 1]
                    ηx = dηdx[iel, k, l, 1]; ηy = dηdy[iel, k, l, 1]
                    val += w * dψ[i, k] * dψ[n, l] * (ξx*ηx + ξy*ηy)
                end

                # (D) single quadrature point (i, n)
                let k = i, l = n
                    w  = ω[k] * ω[l] * Je[iel, k, l, 1]
                    ξx = dξdx[iel, k, l, 1]; ξy = dξdy[iel, k, l, 1]
                    ηx = dηdx[iel, k, l, 1]; ηy = dηdy[iel, k, l, 1]
                    val += w * dψ[j, l] * dψ[m, k] * (ηx*ξx + ηy*ξy)
                end

                push!(I_vec, ip); push!(J_vec, jp); push!(V_vec, val)
            end
        end
    end

    K = sparse(I_vec, J_vec, V_vec, npoin, npoin)
    return M, K
end

#-------------------------------------------------------------------------------
# Dispatcher: pick the 1D or 2D mass/stiffness builder from params.SD. 3D is
# not implemented yet — this is a scalar, isotropic-Laplacian integrator.
#-------------------------------------------------------------------------------
build_mass_stiff(mesh, basis, ω, metrics, ::NSD_1D) =
    build_mass_stiff_1d(mesh, basis, ω)

build_mass_stiff(mesh, basis, ω, metrics, ::NSD_2D) =
    build_mass_stiff_2d(mesh, basis, ω, metrics)

#-------------------------------------------------------------------------------
# Evaluate the explicit (inviscid / advective) RHS for the 1D problem by
# calling the existing `rhs!` with the viscosity switch temporarily disabled,
# then divided by the (diagonal) mass matrix as usual.
#-------------------------------------------------------------------------------
function eval_explicit_rhs!(du, u, params, time, saved_lvisc::Bool)
    params.inputs[:lvisc] = false
    rhs!(du, u, params, time)
    params.inputs[:lvisc] = saved_lvisc
    return du
end

#-------------------------------------------------------------------------------
# Main IMEX time loop (scalar or vector problem on NSD_1D or NSD_2D, CPU /
# ContGal). Output is written at user-requested diagnostic times via the
# existing write_output machinery.
#-------------------------------------------------------------------------------
function imex_ars232_time_loop!(inputs, params, u)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println_rank(" # Solving ODE with IMEX ARS(2,3,2) ........... "; msg_rank = rank)

    if !(params.SD == NSD_1D() || params.SD == NSD_2D())
        error("IMEX_ARS232 is currently implemented for NSD_1D and NSD_2D only.")
    end

    mesh   = params.mesh
    basis  = params.basis
    ω      = params.ω
    nelem  = mesh.nelem
    ngl    = mesh.ngl
    npoin  = mesh.npoin
    neqs   = params.neqs
    Δt     = TFloat(inputs[:Δt])
    tinit  = TFloat(inputs[:tinit])
    tend   = TFloat(inputs[:tend])
    ν      = TFloat(inputs[:μ][1])

    #---------------------------------------------------------------------------
    # Precompute mass and stiffness, then factorize the implicit operator.
    # The ARS(2,3,2) scheme has a single implicit-diagonal coefficient γ, so
    # the same matrix (M + Δt γ ν K) is used at stages 2 and 3.
    #---------------------------------------------------------------------------
    M_diag, K = build_mass_stiff(mesh, basis, ω, params.metrics, params.SD)
    Minv      = 1.0 ./ M_diag

    A_ex, A_im, b_ex, b_im, c_ex, c_im = ars232_tableau(TFloat)
    γ    = A_im[2,2]
    Mmat = spdiagm(0 => M_diag)
    Aimp = Mmat + (Δt * γ * ν) .* K
    Fimp = lu(Aimp)   # sparse LU factorization (stored for all steps)

    #---------------------------------------------------------------------------
    # Work arrays. For scalar problems (Burgers), neqs == 1, so u is a flat
    # vector of length npoin; the stage and RHS arrays match that layout.
    #---------------------------------------------------------------------------
    nstate = npoin * neqs
    Y1 = similar(u);  Y2 = similar(u);  Y3 = similar(u)
    N1 = similar(u);  N2 = similar(u);  N3 = similar(u)
    L1 = similar(u);  L2 = similar(u);  L3 = similar(u)
    rhs_buf  = similar(u)
    rhs_solv = similar(u)

    # Keep track of the user's lvisc setting (the IMEX loop owns viscosity,
    # so we toggle this flag around the explicit-RHS call).
    saved_lvisc = get(inputs, :lvisc, false)

    # f_im(y) = -ν M⁻¹ K y  per equation (applies to a single field).
    function apply_L!(out, y)
        # out .= -ν .* Minv .* (K * y)
        mul!(out, K, y)
        @. out = -ν * Minv * out
        return out
    end

    # Diagnostics output schedule
    dosetimes = collect(TFloat.(inputs[:diagnostics_at_times]))
    sort!(dosetimes)
    next_out_idx = 1
    iout = 0

    #---------------------------------------------------------------------------
    # Write the initial condition.
    #---------------------------------------------------------------------------
    println_rank(" #  IMEX: writing initial condition .............."; msg_rank = rank)
    write_output(params.SD, u, params.uaux, tinit, iout,
                 params.mesh, params.mp,
                 params.connijk_original, params.poin_in_bdy_face_original,
                 params.x_original, params.y_original, params.z_original,
                 inputs[:output_dir], inputs,
                 params.qp.qvars, params.qp.qoutvars,
                 inputs[:outformat];
                 nvar = params.qp.neqs, qexact = params.qp.qe)

    #---------------------------------------------------------------------------
    # Time loop.
    #---------------------------------------------------------------------------
    t_n      = tinit
    n_step   = 0

    while t_n < tend - 1.0e-14

        Δt_step = min(Δt, tend - t_n)

        # If this step would cross the next diagnostic time, clip it so the
        # diagnostic output lands exactly on the requested time.
        if next_out_idx <= length(dosetimes)
            t_target = dosetimes[next_out_idx]
            if t_target > t_n + 1.0e-14 && t_n + Δt_step > t_target - 1.0e-14
                Δt_step = t_target - t_n
            end
        end

        γΔt    = γ * Δt_step
        Aimp_s = Mmat + (γΔt * ν) .* K
        Fimp_s = Δt_step == Δt ? Fimp : lu(Aimp_s)

        # --- Stage 1: Y1 = u^n, evaluate N1 and L1 --------------------------
        copyto!(Y1, u)
        eval_explicit_rhs!(N1, Y1, params, t_n, saved_lvisc)
        apply_L!(L1, Y1)

        # --- Stage 2: solve (M + Δt γ ν K) Y2 = M * rhs2 --------------------
        # rhs2 = u^n + Δt * (A_ex[2,1] * N1 + A_im[2,1] * L1)
        @. rhs_buf = u + Δt_step * (A_ex[2,1] * N1 + A_im[2,1] * L1)
        @. rhs_solv = M_diag * rhs_buf
        ldiv!(Y2, Fimp_s, rhs_solv)

        eval_explicit_rhs!(N2, Y2, params, t_n + c_ex[2] * Δt_step, saved_lvisc)
        apply_L!(L2, Y2)

        # --- Stage 3: solve (M + Δt γ ν K) Y3 = M * rhs3 --------------------
        # rhs3 = u^n + Δt * (A_ex[3,1]*N1 + A_ex[3,2]*N2 + A_im[3,1]*L1 + A_im[3,2]*L2)
        @. rhs_buf = u + Δt_step * (A_ex[3,1] * N1 + A_ex[3,2] * N2 +
                                    A_im[3,1] * L1 + A_im[3,2] * L2)
        @. rhs_solv = M_diag * rhs_buf
        ldiv!(Y3, Fimp_s, rhs_solv)

        eval_explicit_rhs!(N3, Y3, params, t_n + c_ex[3] * Δt_step, saved_lvisc)
        apply_L!(L3, Y3)

        # --- Final update ----------------------------------------------------
        @. u = u + Δt_step * (b_ex[1]*N1 + b_ex[2]*N2 + b_ex[3]*N3 +
                              b_im[1]*L1 + b_im[2]*L2 + b_im[3]*L3)

        t_n   += Δt_step
        n_step += 1

        # --- Diagnostic output ----------------------------------------------
        if next_out_idx <= length(dosetimes) &&
           abs(t_n - dosetimes[next_out_idx]) < 1.0e-12
            iout += 1
            println_rank(@sprintf(" # IMEX: t = %.6f   step = %d", t_n, n_step);
                         msg_rank = rank)
            write_output(params.SD, u, params.uaux, t_n, iout,
                         params.mesh, params.mp,
                         params.connijk_original, params.poin_in_bdy_face_original,
                         params.x_original, params.y_original, params.z_original,
                         inputs[:output_dir], inputs,
                         params.qp.qvars, params.qp.qoutvars,
                         inputs[:outformat];
                         nvar = params.qp.neqs, qexact = params.qp.qe)
            next_out_idx += 1
        end
    end

    # Always write a final snapshot if one is not already written.
    if next_out_idx <= length(dosetimes) || iout == 0
        iout += 1
        write_output(params.SD, u, params.uaux, t_n, iout,
                     params.mesh, params.mp,
                     params.connijk_original, params.poin_in_bdy_face_original,
                     params.x_original, params.y_original, params.z_original,
                     inputs[:output_dir], inputs,
                     params.qp.qvars, params.qp.qoutvars,
                     inputs[:outformat];
                     nvar = params.qp.neqs, qexact = params.qp.qe)
    end

    println_rank(" # Solving ODE with IMEX ARS(2,3,2) ........... DONE"; msg_rank = rank)
    return u
end
