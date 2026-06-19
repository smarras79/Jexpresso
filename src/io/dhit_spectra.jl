import FFTW

# ============================================================================
# Decaying Homogeneous Isotropic Turbulence (DHIT) energy spectra
#
# Computes 1D and 3D kinetic-energy spectra of the DGSEM velocity field the way
# Flad & Gassner (J. Comput. Phys. 350, 2017, "On the use of kinetic energy
# preserving DG-schemes for large eddy simulation") report them for the DHIT
# test case:
#
#   * The piecewise-polynomial DG solution (LGL collocation nodes) is
#     super-sampled onto a *global equidistant* Cartesian grid of
#         N_uni = nel_per_dir * (N+1)
#     points per direction — exactly the "DOF per direction" of the paper
#     (e.g. 6 cells x N=7  ->  48^3 DOF; 18 cells x N=7  ->  144^3 DOF).
#     The interpolation is an element-local tensor-product Lagrange
#     interpolation from the LGL nodes to the equidistant points.
#   * A 3D FFT is taken of the three velocity components on that uniform grid.
#   * The 3D KE spectrum is the spherical-shell sum
#         E(k) = sum_{ k-1/2 <= |kappa| < k+1/2 } 1/2 (|u^(kappa)|^2 + ...),
#     reported for integer shells up to the Nyquist wavenumber N_uni/2.
#   * The 1D spectrum is the energy summed over the two perpendicular
#     wavenumber planes, averaged over the three homogeneous directions.
#
# Both spectra integrate to the resolved kinetic energy  sum_k E(k) = <|u|^2>/2.
#
# This routine assumes a *uniform, axis-aligned, periodic* hexahedral mesh
# (the DHIT cube), which is what the test case uses. It is MPI-aware: every
# rank interpolates its own elements onto their equidistant sub-blocks, the
# sub-blocks are gathered to rank 0, and rank 0 assembles the global grid,
# performs the FFTs, and writes the output files.
# ============================================================================

# 1D Lagrange interpolation matrix M[m, i] = L_i(ξt_m) from nodal points ξ
# (length ngl) to target points ξt (length nt).
function _lagrange_interp_matrix(ξ::Vector{Float64}, ξt::Vector{Float64})
    ngl = length(ξ)
    nt  = length(ξt)
    M   = zeros(Float64, nt, ngl)
    @inbounds for m in 1:nt
        for i in 1:ngl
            val = 1.0
            for l in 1:ngl
                l == i && continue
                val *= (ξt[m] - ξ[l]) / (ξ[i] - ξ[l])
            end
            M[m, i] = val
        end
    end
    return M
end

# Reference LGL node positions (in [-1,1], ordered by the connijk index) for
# each axis, derived directly from the coordinates of a reference element so we
# do not depend on the basis carrying its quadrature points.
function _ref_nodes_from_element(connijk, ngl::Int, x, y, z, iel::Int)
    ξx = zeros(Float64, ngl); ξy = zeros(Float64, ngl); ξz = zeros(Float64, ngl)
    for i in 1:ngl; ξx[i] = x[connijk[iel, i, 1, 1]]; end
    for j in 1:ngl; ξy[j] = y[connijk[iel, 1, j, 1]]; end
    for k in 1:ngl; ξz[k] = z[connijk[iel, 1, 1, k]]; end
    norm!(v) = (lo = minimum(v); hi = maximum(v); @. v = 2.0*(v - lo)/(hi - lo) - 1.0)
    norm!(ξx); norm!(ξy); norm!(ξz)
    return ξx, ξy, ξz
end

"""
    compute_dhit_spectra(u, params, t, iout)

Compute and write the 1D and 3D DHIT kinetic-energy spectra at output index
`iout` / time `t`. No-op unless `params.inputs[:dhit_spectra] == true`.
"""
function compute_dhit_spectra(u, params, t, iout)
    get(params.inputs, :dhit_spectra, false) || return nothing

    comm  = get_mpi_comm()
    rank  = MPI.Comm_rank(comm)
    mesh  = params.mesh
    ngl   = mesh.ngl
    neqs  = params.neqs
    npoin = mesh.npoin

    # Recover primitive velocities at every node from the conserved state.
    uaux = params.uaux
    u2uaux!(@view(uaux[:, :]), u, neqs, npoin)

    # --- Global domain bounds and element count per direction ---
    x = Array(mesh.x); y = Array(mesh.y); z = Array(mesh.z)
    xmin = MPI.Allreduce(minimum(x), MPI.MIN, comm)
    xmax = MPI.Allreduce(maximum(x), MPI.MAX, comm)
    ymin = MPI.Allreduce(minimum(y), MPI.MIN, comm)
    ymax = MPI.Allreduce(maximum(y), MPI.MAX, comm)
    zmin = MPI.Allreduce(minimum(z), MPI.MIN, comm)
    zmax = MPI.Allreduce(maximum(z), MPI.MAX, comm)
    Lx = xmax - xmin; Ly = ymax - ymin; Lz = zmax - zmin

    gnelem   = mesh.gnelem
    nel_dir  = round(Int, cbrt(gnelem))
    if nel_dir^3 != gnelem
        rank == 0 && @warn "DHIT spectra: global element count $gnelem is not a perfect cube; spectra skipped."
        return nothing
    end
    hex = Lx / nel_dir; hey = Ly / nel_dir; hez = Lz / nel_dir
    N_uni = nel_dir * ngl                      # equidistant points per direction (= DOF/dir)

    # --- Element-local interpolation: LGL nodes -> equidistant sub-grid ---
    # Targets inside one element (ref coords): ξt_m = -1 + 2(m-1)/ngl, m=1..ngl,
    # so the ngl sub-points of consecutive elements tile a periodic uniform grid.
    ξt = [-1.0 + 2.0*(m-1)/ngl for m in 1:ngl]
    connijk = Array(mesh.connijk)
    local_buf = Float64[]
    if mesh.nelem > 0
        ξx, ξy, ξz = _ref_nodes_from_element(connijk, ngl, x, y, z, 1)
        Mx = _lagrange_interp_matrix(ξx, ξt)
        My = _lagrange_interp_matrix(ξy, ξt)
        Mz = _lagrange_interp_matrix(ξz, ξt)

        uel = zeros(Float64, ngl, ngl, ngl, 3)
        t1  = zeros(Float64, ngl, ngl, ngl, 3)   # after x-contraction
        t2  = zeros(Float64, ngl, ngl, ngl, 3)   # after y-contraction
        blk = zeros(Float64, ngl, ngl, ngl, 3)   # equidistant sub-block
        sizehint!(local_buf, mesh.nelem * (3 + 3*ngl^3))

        for iel in 1:mesh.nelem
            # structured element index (0-based) from its lower corner
            xlo = x[connijk[iel, 1, 1, 1]]; ylo = y[connijk[iel, 1, 1, 1]]; zlo = z[connijk[iel, 1, 1, 1]]
            @inbounds for i in 1:ngl, j in 1:ngl, k in 1:ngl
                ip = connijk[iel, i, j, k]
                xlo = min(xlo, x[ip]); ylo = min(ylo, y[ip]); zlo = min(zlo, z[ip])
            end
            ex = round(Int, (xlo - xmin)/hex)
            ey = round(Int, (ylo - ymin)/hey)
            ez = round(Int, (zlo - zmin)/hez)

            @inbounds for k in 1:ngl, j in 1:ngl, i in 1:ngl
                ip = connijk[iel, i, j, k]
                ρ  = uaux[ip, 1]
                uel[i, j, k, 1] = uaux[ip, 2]/ρ
                uel[i, j, k, 2] = uaux[ip, 3]/ρ
                uel[i, j, k, 3] = uaux[ip, 4]/ρ
            end

            # tensor-product interpolation: contract x, then y, then z
            fill!(t1, 0.0); fill!(t2, 0.0); fill!(blk, 0.0)
            @inbounds for c in 1:3, k in 1:ngl, j in 1:ngl, m in 1:ngl
                s = 0.0
                for i in 1:ngl; s += Mx[m, i]*uel[i, j, k, c]; end
                t1[m, j, k, c] = s
            end
            @inbounds for c in 1:3, k in 1:ngl, n in 1:ngl, m in 1:ngl
                s = 0.0
                for j in 1:ngl; s += My[n, j]*t1[m, j, k, c]; end
                t2[m, n, k, c] = s
            end
            @inbounds for c in 1:3, p in 1:ngl, n in 1:ngl, m in 1:ngl
                s = 0.0
                for k in 1:ngl; s += Mz[p, k]*t2[m, n, k, c]; end
                blk[m, n, p, c] = s
            end

            push!(local_buf, Float64(ex), Float64(ey), Float64(ez))
            @inbounds for c in 1:3, p in 1:ngl, n in 1:ngl, m in 1:ngl
                push!(local_buf, blk[m, n, p, c])
            end
        end
    end

    # --- Gather equidistant sub-blocks to rank 0 ---
    all_bufs = MPI.gather(local_buf, comm)
    rank != 0 && return nothing

    n   = N_uni
    Ug  = zeros(Float64, n, n, n)
    Vg  = zeros(Float64, n, n, n)
    Wg  = zeros(Float64, n, n, n)
    rec = 3 + 3*ngl^3
    for buf in all_bufs
        nrec = length(buf) ÷ rec
        for e in 1:nrec
            base = (e-1)*rec
            ex = round(Int, buf[base+1]); ey = round(Int, buf[base+2]); ez = round(Int, buf[base+3])
            off = base + 3
            @inbounds for c in 1:3, p in 1:ngl, n2 in 1:ngl, m in 1:ngl
                off += 1
                gi = ex*ngl + m; gj = ey*ngl + n2; gk = ez*ngl + p
                if     c == 1; Ug[gi, gj, gk] = buf[off]
                elseif c == 2; Vg[gi, gj, gk] = buf[off]
                else           Wg[gi, gj, gk] = buf[off]
                end
            end
        end
    end

    # --- FFTs and spectra ---
    Uh = FFTW.fft(Ug); Vh = FFTW.fft(Vg); Wh = FFTW.fft(Wg)
    invN3 = 1.0 / (Float64(n)^3)
    dk    = 2π / Lx                              # spectral resolution (=1 for L=2π)

    kidx(i) = (i-1) <= n÷2 ? (i-1) : (i-1-n)     # signed integer wavenumber index
    kmax3d  = floor(Int, sqrt(3.0)*(n÷2)) + 1
    E3d  = zeros(Float64, kmax3d + 1)            # index 0..kmax3d
    kmax1d = n÷2
    E1x = zeros(Float64, kmax1d + 1)
    E1y = zeros(Float64, kmax1d + 1)
    E1z = zeros(Float64, kmax1d + 1)

    @inbounds for kz in 1:n
        k3 = kidx(kz)
        for ky in 1:n
            k2 = kidx(ky)
            for kx in 1:n
                k1 = kidx(kx)
                au = Uh[kx, ky, kz]*invN3
                av = Vh[kx, ky, kz]*invN3
                aw = Wh[kx, ky, kz]*invN3
                e  = 0.5*(abs2(au) + abs2(av) + abs2(aw))
                e == 0.0 && continue
                # 3D spherical-shell binning
                m = round(Int, sqrt(Float64(k1*k1 + k2*k2 + k3*k3)))
                E3d[m + 1] += e
                # 1D spectra: energy in each direction's |wavenumber| plane
                E1x[abs(k1) + 1] += e
                E1y[abs(k2) + 1] += e
                E1z[abs(k3) + 1] += e
            end
        end
    end

    TKE = sum(E3d)

    # --- Dissipation-rate spectrum D(k) = 2 ν k² E(k) and band integrals ---
    # ν is the (fixed) molecular kinematic viscosity used for both DNS and LES,
    # so the energy/dissipation integrals are directly comparable, exactly as in
    # Flad & Gassner (2017), Fig. 9. ρ_ref = 1 for this unit-density DHIT setup.
    inputs = params.inputs
    μmol = get(inputs, :mu_molecular, length(get(inputs, :μ, [])) >= 2 ? inputs[:μ][2] : 0.0)
    ν    = Float64(get(inputs, :dhit_nu, μmol))                 # kinematic ν = μ/ρ_ref, ρ_ref=1
    kcut = Float64(get(inputs, :dhit_kcut, 16.0))               # integration cutoff wavenumber
    mcut = max(1, round(Int, kcut/dk))

    Dk = zeros(Float64, kmax3d + 1)                             # 2 ν k² E(k)
    @inbounds for m in 1:kmax3d
        kp = m*dk
        Dk[m + 1] = 2.0*ν*kp*kp*E3d[m + 1]
    end
    # Band integrals up to kcut (shells already carry the per-shell energy, so a
    # plain sum is the Fourier-space integral ∫_1^{kcut} … dk):
    Ekin_band = 0.0; eps_band = 0.0
    @inbounds for m in 1:min(mcut, kmax3d)
        Ekin_band += E3d[m + 1]
        eps_band  += Dk[m + 1]
    end
    eps_total = sum(Dk)

    # --- Write output ---
    outdir = inputs[:output_dir]
    f3 = joinpath(outdir, @sprintf("dhit_spectrum_3D_%06d.dat", iout))
    open(f3, "w") do io
        @printf(io, "# DHIT 3D KE spectrum   time=%.6e  N_uni=%d  Nyquist=%d  dk=%.6e  nu=%.6e  TKE=%.8e\n",
                t, n, n÷2, dk, ν, TKE)
        println(io, "# k            E(k)           D(k)=2*nu*k^2*E(k)")
        for m in 1:kmax3d
            @printf(io, "%.6e  %.8e  %.8e\n", m*dk, E3d[m + 1], Dk[m + 1])
        end
    end

    # --- Time series of band integrals (Fig. 9): appended across output times ---
    fint = joinpath(outdir, "dhit_integrals.dat")
    write_header = !isfile(fint) || iout == 0
    open(fint, write_header ? "w" : "a") do io
        if write_header
            @printf(io, "# DHIT energy/dissipation integrals   nu=%.6e  kcut=%.3f  dk=%.6e\n", ν, kcut, dk)
            println(io, "# time          Ekin(k<=kcut)   eps(k<=kcut)    Ekin_total      eps_total")
        end
        @printf(io, "%.6e  %.8e  %.8e  %.8e  %.8e\n", t, Ekin_band, eps_band, TKE, eps_total)
    end

    f1 = joinpath(outdir, @sprintf("dhit_spectrum_1D_%06d.dat", iout))
    open(f1, "w") do io
        @printf(io, "# DHIT 1D KE spectra (per-direction and average)   time=%.6e  N_uni=%d  dk=%.6e\n",
                t, n, dk)
        println(io, "# k            E1_avg         E1_x           E1_y           E1_z")
        for m in 1:kmax1d
            avg = (E1x[m + 1] + E1y[m + 1] + E1z[m + 1]) / 3.0
            @printf(io, "%.6e  %.8e  %.8e  %.8e  %.8e\n",
                    m*dk, avg, E1x[m + 1], E1y[m + 1], E1z[m + 1])
        end
    end

    println(" # DHIT spectra written: ", basename(f3), ", ", basename(f1),
            "  (N_uni=", n, ", TKE=", @sprintf("%.6e", TKE), ")")
    return nothing
end
