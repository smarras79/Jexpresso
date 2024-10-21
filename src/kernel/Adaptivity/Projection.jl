include(joinpath( "..", "bases", "basis_structs.jl"))

# Define a mutable struct to hold the Legendre polynomial and derivatives
mutable struct LegendrePoly
    p0::Float64   # Nth order Legendre polynomial
    p0_1::Float64 # First derivative of p0
    p0_2::Float64 # Second derivative of p0
    p1::Float64   # (N-1)th order Legendre polynomial
    p1_1::Float64 # First derivative of p1
    p1_2::Float64 # Second derivative of p1
    p2::Float64   # (N-2)th order Legendre polynomial
    p2_1::Float64 # First derivative of p2
    p2_2::Float64 # Second derivative of p2
    p00::Float64  # (N+1)th order Legendre polynomial
    p00_1::Float64 # First derivative of p00
    p00_2::Float64 # Second derivative of p00
end

# Function to construct the Legendre polynomial and its derivatives
function legendre_poly_loc(n::Int, x::Float64)::LegendrePoly
    # Initialize all polynomials and their derivatives
    p2, p2_1, p2_2 = 0.0, 0.0, 0.0
    p1, p1_1, p1_2 = 0.0, 0.0, 0.0
    p0, p0_1, p0_2 = 1.0, 0.0, 0.0
    p00, p00_1, p00_2 = 1.0, 0.0, 0.0

    # Construct the Nth order Legendre polynomial using recurrence relations
    for j in 1:n
        # Save previous values
        p2, p2_1, p2_2 = p1, p1_1, p1_2
        p1, p1_1, p1_2 = p0, p0_1, p0_2

        # Recurrence coefficients
        a = (2.0 * j - 1.0) / j
        b = (j - 1.0) / j
        # Update p0 and its derivatives
        p0 = a * x * p1 - b * p2
        p0_1 = a * (p1 + x * p1_1) - b * p2_1
        p0_2 = a * (2.0 * p1_1 + x * p1_2) - b * p2_2

        # For p00 (N+1 order polynomial)
        a = (2.0 * j + 1.0) / (j + 1.0)
        b = j / (j + 1.0)
        p00 = a * x * p0 - b * p1
        p00_1 = a * (p0 + x * p0_1) - b * p1_1
        p00_2 = a * (2.0 * p0_1 + x * p0_2) - b * p1_2
    end

    # Return the Legendre polynomial and its derivatives as a struct
    return LegendrePoly(p0, p0_1, p0_2, p1, p1_1, p1_2, p2, p2_1, p2_2, p00, p00_1, p00_2)
end


# Lagrange Polynomial Construction
function lagrange_poly(L, xq, xgl, nq, ngl)
    L .= 1.0

    for i in 1:ngl
        for j in 1:ngl
            if i != j
                for k in 1:nq
                    L[i, k] *= (xq[k] - xgl[j]) / (xgl[i] - xgl[j])
                end
            end
        end
    end
end

# Main Function (translated from Fortran subroutine create_2d_projection_matrices_numa2d)
function scatter_gather_projection!(plane, nglx, ngly, nglz)
    ngl1, ngl2 = (plane == 1) ? (nglx, ngly) :
                 (plane == 2) ? (nglx, nglz) :
                 (plane == 3) ? (ngly, nglz) : (0, 0)

    nq1 = ngl1 + 1
    nq2 = ngl2 + 1
    nngl = ngl1 * ngl2

    Psg = zeros(Float64, nngl, nngl, 8)
    M = zeros(Float64, nngl, nngl)
    xlgl1 = basis_structs_ξ_ω!(LGL(), ngl1-1, CPU())
    xlgl2 = basis_structs_ξ_ω!(LGL(), ngl2-1, CPU())
    qlgl1 = basis_structs_ξ_ω!(LGL(), nq1-1, CPU())
    qlgl2 = basis_structs_ξ_ω!(LGL(), nq2-1, CPU())

    xgl1, wgl1 = xlgl1.ξ, xlgl1.ω
    xgl2, wgl2 = xlgl2.ξ, xlgl2.ω
    xq1, wq1 = qlgl1.ξ, qlgl1.ω
    xq2, wq2 = qlgl2.ξ, qlgl2.ω

    xq3, xq4 = copy(xq1), copy(xq2)
    L1o, L2o = zeros(Float64, ngl1, nq1), zeros(Float64, ngl2, nq2)
    L1, L2, L3, L4 = zeros(Float64, ngl1, nq1), zeros(Float64, ngl2, nq2), zeros(Float64, ngl1, nq1), zeros(Float64, ngl2, nq2)

    lagrange_poly(L1o, xq1, xgl1, nq1, ngl1)
    lagrange_poly(L2o, xq2, xgl2, nq2, ngl2)

    # Scaling factors
    s = 0.5
    o1 = -0.5
    o2 = 0.5

    # Offset quadrature points
    xq1 .= o1 .+ s .* xq1
    xq3 .= o2 .+ s .* xq3
    xq2 .= o1 .+ s .* xq2
    xq4 .= o2 .+ s .* xq4

    lagrange_poly(L1, xq1, xgl1, nq1, ngl1)
    lagrange_poly(L2, xq2, xgl2, nq2, ngl2)
    lagrange_poly(L3, xq3, xgl1, nq1, ngl1)
    lagrange_poly(L4, xq4, xgl2, nq2, ngl2)

    # Initialize matrices
    M .= 0.0
    Psg .= 0.0

    # Compute projection matrices
    for l in 1:ngl1, k in 1:ngl2
        mm = (l - 1) * ngl2 + k
        for j in 1:ngl1, i in 1:ngl2
            nn = (j - 1) * ngl2 + i
            for q in 1:nq1, p in 1:nq2
                M[nn, mm] += L2o[i, p] * L1o[j, q] * L2o[k, p] * L1o[l, q] * wq2[p] * wq1[q]
                Psg[nn, mm, 1] += L2[i, p] * L1[j, q] * L2o[k, p] * L1o[l, q] * wq2[p] * wq1[q]
                Psg[nn, mm, 2] += L4[i, p] * L1[j, q] * L2o[k, p] * L1o[l, q] * wq2[p] * wq1[q]
                Psg[nn, mm, 3] += L2[i, p] * L3[j, q] * L2o[k, p] * L1o[l, q] * wq2[p] * wq1[q]
                Psg[nn, mm, 4] += L4[i, p] * L3[j, q] * L2o[k, p] * L1o[l, q] * wq2[p] * wq1[q]
                Psg[nn, mm, 5] += L2o[i, p] * L1o[j, q] * L2[k, p] * L1[l, q] * wq2[p] * wq1[q] * s * s
                Psg[nn, mm, 6] += L2o[i, p] * L1o[j, q] * L4[k, p] * L1[l, q] * wq2[p] * wq1[q] * s * s
                Psg[nn, mm, 7] += L2o[i, p] * L1o[j, q] * L2[k, p] * L3[l, q] * wq2[p] * wq1[q] * s * s
                Psg[nn, mm, 8] += L2o[i, p] * L1o[j, q] * L4[k, p] * L3[l, q] * wq2[p] * wq1[q] * s * s
            end
        end
    end

    # Matrix inversion using LU decomposition
    M_inv = inv(M)

    # Final projection matrix computation
    for i in 1:8
        Psg[:, :, i] = Psg[:, :, i] * M_inv
    end

    # Fix for 2D
    if ngl2 == 1
        Psg[:, :, 1] = 0.5 * (Psg[:, :, 1] + Psg[:, :, 2])
        Psg[:, :, 2] = 0.5 * (Psg[:, :, 3] + Psg[:, :, 4])
        Psg[:, :, 3] = Psg[:, :, 5] + Psg[:, :, 6]
        Psg[:, :, 4] = Psg[:, :, 7] + Psg[:, :, 8]
    elseif ngl1 == 1
        Psg[:, :, 1] = 0.5 * (Psg[:, :, 1] + Psg[:, :, 3])
        Psg[:, :, 2] = 0.5 * (Psg[:, :, 2] + Psg[:, :, 4])
        Psg[:, :, 3] = Psg[:, :, 5] + Psg[:, :, 7]
        Psg[:, :, 4] = Psg[:, :, 6] + Psg[:, :, 8]
    end

    return Psg
end


using Test

function test_create_2d_projection_matrices_numa2d()
    # Test parameters
    nglx = 2
    ngly = 2
    nglz = 2
    plane = 1  # Testing for plane 1 (X-Y plane)
    
    # Initialize Psg to be an empty array
    Psg = zeros(Float64, 0, 0, 0)  # Will be resized inside the function
    
    # Call the function with test parameters
    Psg = scatter_gather_projection!(plane, nglx, ngly, nglz)
    
    # Test if the output matrix has the correct dimensions
    nngl = nglx * ngly
    @test size(Psg) == (nngl, nngl, 8)

    # Test if the matrix contains non-zero values
    @test count(!iszero, Psg) > 0  # Ensure that there are non-zero values in Psg

    # Check specific properties of the Psg matrix if known (e.g., symmetry, positive values)
    # For simplicity, we can check that none of the values in Psg are NaN or Inf
    @test all(isfinite, Psg)  # Ensure no NaN or Inf values in Psg
    @info Psg

    println("All tests for create_2d_projection_matrices_numa2d passed successfully!")
end

function build_interpolation(ra, rb, wa)
    Nra = length(ra)
    Nrb = length(rb)

    # Initialize the interpolation matrix
    I_a2b = zeros(Float64, Nrb, Nra)

    # Loop over all points in rb
    for k in 1:Nrb
        d = 0.0
        # Loop over all points in ra
        for j in 1:Nra
            # Check if ra and rb values are approximately equal
            if isapprox(rb[k], ra[j], atol=sqrt(eps(Float64)))
                # Set the entire row to 0
                I_a2b[k, :] .= 0.0
                # Set the matching element to 1
                I_a2b[k, j] = 1.0
                d = 1.0
                break
            end
            # Compute the value for I_a2b
            I_a2b[k, j] = wa[j] / (rb[k] - ra[j])
            d += wa[j] / (rb[k] - ra[j])
        end

        # Normalize the row by dividing by d
        if d != 1.0  # Avoid redundant division for the case where d was set to 1.0
            for j in 1:Nra
                I_a2b[k, j] /= d
            end
        end
    end

    return I_a2b
end


function barycentric_weights(ra)
    Nra = length(ra)
    wa = ones(Float64, Nra)  # Initialize weights to 1

    # Loop over all points in ra to compute barycentric weights
    for k in 1:Nra
        for j in 1:Nra
            if j != k
                wa[k] *= (ra[k] - ra[j])
            end
        end
        wa[k] = 1.0 / wa[k]  # Final weight is the reciprocal
    end

    return wa
end


function build_projection_1d(ξa)
    Np = length(ξa)

    # Create the barycentric weights (we assume a similar function exists in Julia)
    ωa = barycentric_weights(ξa)

    # Create top and bottom grids
    ξa_t = (ξa .+ 1) ./ 2
    ξa_b = (ξa .- 1) ./ 2

    # Initialize interpolation matrices
    interp = zeros(Float64, Np, Np, 2)

    # Build interpolation matrices
    interp[:, :, 2] = build_interpolation(ξa, ξa_b, ωa)
    interp[:, :, 1] = build_interpolation(ξa, ξa_t, ωa)

    # Build the mass matrix and its inverse
    lgb =  basis_structs_ξ_ω!(LG(),Np-1, CPU())  # Similar to Fortran's legendre_gauss function
    ξb, ωb_gl = lgb.ξ, lgb.ω
    # Get barycentric weights for ξb
    ωb = barycentric_weights(ξb)

    # Build interpolation matrices I_a2b and I_b2a
    I_a2b = build_interpolation(ξa, ξb, ωa)
    I_b2a = build_interpolation(ξb, ξa, ωb)

    # Initialize mass matrix M and its inverse MI
    M = diagm(ωb_gl)   # Mass matrix
    MI = diagm(1.0 ./ ωb_gl)  # Inverse of the mass matrix

    # Update M and MI using interpolation matrices
    M = transpose(I_a2b) * M * I_a2b
    MI = I_b2a * MI * transpose(I_b2a)

    # Create the projection matrices
    project = zeros(Float64, Np, Np, 2)
    project[:, :, 1] = MI * (transpose(interp[:, :, 1]) * M) / 2
    project[:, :, 2] = MI * (transpose(interp[:, :, 2]) * M) / 2

    return interp, project
end



function p8est_transfer_q!(q_dst, q_src, lvl_src, lvl_dst, mesh_dst, mesh_src, n2o_ele_map, SD::NSD_3D)

    # Dimensions of the source and destination grids (global variables nglx, ngly, nglz are assumed defined)
    nx = mesh_dst.ngl
    ny = mesh_dst.ngl
    nz = mesh_dst.ngl
    np = nx * ny * nz  # Number of points
    nf = size(q_dst, 2)  # Number of fields
    child_order = [[1,2,2],[2,2,2],[1,2,1],[2,2,1],[1,1,2],[2,1,2],[1,1,1],[2,1,1]]

    num_elem = length(lvl_dst)  # Number of elements in destination

    # Temporary storage for interpolation and projection
    qx = zeros(Float64, np, nf)
    qxy = zeros(Float64, np, nf)

    lgl = basis_structs_ξ_ω!(LGL(), mesh_dst.nop, CPU())

    # Ensure the interpolation matrices are computed
    # if !@isdefined(interp_x)
    #     global interp_x, project_x, interp_y, project_y, interp_z, project_z
    interp_x, project_x = build_projection_1d(lgl.ξ)
    interp_y, project_y = build_projection_1d(lgl.ξ)
    interp_z, project_z = build_projection_1d(lgl.ξ)
    # end

    # Initialize source and destination element indices
    k_dst, o_dst = 1, 0
    k_src, o_src = 1, 0
    ip_dst = 0
    ip_src = 0
    @info n2o_ele_map
    # Main transfer loop
    if typeof(n2o_ele_map) == Vector
        while k_dst <= num_elem
            k_src = n2o_ele_map[k_dst]
            
            # Source and Destination are at the same level
            if lvl_src[k_src] == lvl_dst[k_dst]
                for k=1:nz, j=1:ny, i=1:nx
                    for f in 1:nf
                        ip_dst = mesh_dst.connijk[k_dst,i,j,k]
                        ip_src = mesh_src.connijk[k_src,i,j,k]
                        q_dst[ip_dst, f] = q_src[ip_src, f]
                    end
                end
                k_dst += 1

            # Destination is refined (hanging node)
            elseif lvl_src[k_src] < lvl_dst[k_dst]
                # Interpolate for each child
                for order in child_order
                    cx = order[1]
                    cy = order[2]
                    cz = order[3]
                    qx .= 0.0
                    qxy .= 0.0

                    # Zero out q_dst for this element
                    for k=1:nz, j=1:ny, i=1:nx
                        for f in 1:nf
                            ip_dst = mesh_dst.connijk[k_dst,i,j,k]
                            ip_src = mesh_src.connijk[k_src,i,j,k]
                            q_dst[ip_dst, f] = 0.0
                        end
                    end

                    # Interpolation in x: (I ⊗ I ⊗ A) * v
                    for iz in 0:nz-1, iy in 0:ny-1, ix in 0:nx-1
                        for jx in 0:nx-1
                            for f in 1:nf
                                ip_src = mesh_src.connijk[k_src,ix+1,iy+1,iz+1]
                                qx[(iz * ny + iy) * nx + jx + 1, f] += interp_x[jx + 1, ix + 1, cx] * q_src[ip_src, f]
                            end
                        end
                    end

                    # Interpolation in y: (I ⊗ A ⊗ I) * v
                    for iy in 0:ny-1, jy in 0:ny-1
                        for iz in 0:nx-1, ix in 0:nx-1
                            for f in 1:nf
                                # ip_dst = mesh_dst.connijk[k_dst,ix+1,jy+1]
                                qxy[(iz * ny + jy) * nx + ix + 1, f] += interp_y[jy + 1, iy + 1, cy] * qx[(iz * ny + iy) * nx + ix + 1, f]
                            end
                        end
                    end

                    # # Interpolation in z: (A ⊗ I ⊗ I) * v
                    for iz in 0:nz-1, jz in 0:nz-1
                        for iy in 0:ny-1, ix in 0:nx-1
                            for f in 1:nf
                                ip_dst = mesh_dst.connijk[k_dst,ix+1,iy+1,jz+1]
                                q_dst[ip_dst, f] += interp_z[jz + 1, iz + 1, cz] * qxy[(iz * ny + iy) * nx + ix + 1, f]
                            end
                        end
                    end
                    k_dst += 1
                    k_src = n2o_ele_map[k_dst]
                end
                # k_src += 1
            end
        end
    else
        while k_dst <= num_elem
            k_src_vec = n2o_ele_map[k_dst]
            k_src = k_src_vec[1]
            # Source and Destination are at the same level
            if lvl_src[k_src] == lvl_dst[k_dst]
                for k=1:nz, j=1:ny, i=1:nx
                    for f in 1:nf
                        ip_dst = mesh_dst.connijk[k_dst,i,j,k]
                        ip_src = mesh_src.connijk[k_src,i,j,k]
                        q_dst[ip_dst, f] = q_src[ip_src, f]
                    end
                end
                k_dst += 1

            # Destination is refined (hanging node)
            elseif lvl_src[k_src] < lvl_dst[k_dst]
                # Interpolate for each child
                for order in child_order
                    cx = order[1]
                    cy = order[2]
                    cz = order[3]
                    qx .= 0.0
                    qxy .= 0.0

                    # Zero out q_dst for this element
                    for k=1:nz, j=1:ny, i=1:nx
                        for f in 1:nf
                            ip_dst = mesh_dst.connijk[k_dst,i,j,k]
                            ip_src = mesh_src.connijk[k_src,i,j,k]
                            q_dst[ip_dst, f] = 0.0
                        end
                    end

                    # Interpolation in x: (I ⊗ I ⊗ A) * v
                    for iz in 0:nz-1, iy in 0:ny-1, ix in 0:nx-1
                        for jx in 0:nx-1
                            for f in 1:nf
                                ip_src = mesh_src.connijk[k_src,ix+1,iy+1,iz+1]
                                # @info k_src
                                qx[(iz * ny + iy) * nx + jx + 1, f] += interp_x[jx + 1, ix + 1, cx] * q_src[ip_src, f]
                            end
                        end
                    end

                    # Interpolation in y: (I ⊗ A ⊗ I) * v
                    for iz in 0:nz-1, iy in 0:ny-1, jy in 0:ny-1
                        for ix in 0:nx-1
                            for f in 1:nf
                                # ip_dst = mesh_dst.connijk[k_dst,ix+1,jy+1,iz+1]
                                qxy[(iz * ny + jy) * nx + ix + 1, f] += interp_y[jy + 1, iy + 1, cy] * qx[(iz * ny + iy) * nx + ix + 1, f]
                            end
                        end
                    end

                    # # Interpolation in z: (A ⊗ I ⊗ I) * v
                    for iz in 0:nz-1, jz in 0:nz-1
                        for iy in 0:ny-1, ix in 0:nx-1
                            for f in 1:nf
                                ip_dst = mesh_dst.connijk[k_dst,ix+1,iy+1,jz+1]
                                q_dst[ip_dst, f] += interp_z[jz + 1, iz + 1, cz] * qxy[(iz * ny + iy) * nx + ix + 1, f]
                            end
                        end
                    end
                    k_dst += 1
                    k_src = n2o_ele_map[k_dst][1]
                end
                # k_src += 1
            else
                # Source is refined (hanging node)
                # Zero out q_dst for this element
                for k=1:nz, j=1:ny, i=1:nx
                    for f in 1:nf
                        ip_dst = mesh_dst.connijk[k_dst,i,j,k]
                        q_dst[ip_dst, f] = 0.0
                    end
                end

                # Interpolation with projection in x, y, and z
                for order in 1:8
                    k_src_i = k_src_vec[order]
                    cx = child_order[order][1]
                    cy = child_order[order][2]
                    cz = child_order[order][3]
                    qx .= 0.0
                    qxy .= 0.0

                    # Projection in x: (I ⊗ I ⊗ A) * v
                    for iz in 0:nz-1, iy in 0:ny-1, ix in 0:nx-1
                        for jx in 0:nx-1
                            for f in 1:nf
                                ip_src = mesh_src.connijk[k_src_i,ix+1,iy+1,iz+1]
                                qx[(iz * ny + iy) * nx + jx + 1, f] += project_x[jx + 1, ix + 1, cx] * q_src[ip_src, f]
                            end
                        end
                    end

                    # Projection in y: (I ⊗ A ⊗ I) * v
                    for iz in 0:nz-1, iy in 0:ny-1, jy in 0:ny-1
                        for ix in 0:nx-1
                            for f in 1:nf
                                # ip_dst = mesh_dst.connijk[k_dst,ix+1,jy+1]
                                qxy[(iz * ny + jy) * nx + ix + 1, f] += project_y[jy + 1, iy + 1, cy] * qx[(iz * ny + iy) * nx + ix + 1, f]
                            end
                        end
                    end

                    #     # Projection in z: (A ⊗ I ⊗ I) * v
                    for iz in 0:nz-1, jz in 0:nz-1
                        for iy in 0:ny-1, ix in 0:nx-1
                            for f in 1:nf
                                ip_dst = mesh_dst.connijk[k_dst,ix+1,iy+1,jz+1]
                                q_dst[ip_dst, f] += project_z[jz + 1, iz + 1, cz] * qxy[(iz * ny + iy) * nx + ix + 1, f]
                            end
                        end
                    end

                end
                k_dst += 1
            end
        end
    end
end



function p8est_transfer_q!(q_dst, q_src, lvl_src, lvl_dst, mesh_dst, mesh_src, n2o_ele_map, SD::NSD_2D)

    # Dimensions of the source and destination grids (global variables nglx, ngly, nglz are assumed defined)
    nx = mesh_dst.ngl
    ny = mesh_dst.ngl
    np = nx * ny  # Number of points
    nf = size(q_dst, 2)  # Number of fields
    child_order = [[2,1],[2,2],[1,1],[1,2]]

    num_elem = length(lvl_dst)  # Number of elements in destination

    # Temporary storage for interpolation and projection
    qx = zeros(Float64, np, nf)
    # qxy = zeros(Float64, np, nf)

    lgl = basis_structs_ξ_ω!(LGL(), mesh_dst.nop, CPU())

    # Ensure the interpolation matrices are computed
    # if !@isdefined(interp_x)
        # global interp_x, project_x, interp_y, project_y
    interp_x, project_x = build_projection_1d(lgl.ξ)
    interp_y, project_y = build_projection_1d(lgl.ξ)
    # end

    # Initialize source and destination element indices
    k_dst, o_dst = 1, 0
    k_src, o_src = 1, 0
    ip_dst = 0
    ip_src = 0
    # Main transfer loop
    if typeof(n2o_ele_map) == Vector
        while k_dst <= num_elem
            k_src = n2o_ele_map[k_dst]
            
            # Source and Destination are at the same level
            if lvl_src[k_src] == lvl_dst[k_dst]
                for j=1:mesh_dst.ngl, i=1:mesh_dst.ngl
                    for f in 1:nf
                        ip_dst = mesh_dst.connijk[k_dst,i,j]
                        ip_src = mesh_src.connijk[k_src,i,j]
                        q_dst[ip_dst, f] = q_src[ip_src, f]
                    end
                end
                k_dst += 1

            # Destination is refined (hanging node)
            elseif lvl_src[k_src] < lvl_dst[k_dst]
                # Interpolate for each child
                for order in child_order
                    cy = order[2]
                    cx = order[1]
                    qx .= 0.0
                    # qxy .= 0.0

                    # Zero out q_dst for this element
                    for j=1:mesh_dst.ngl, i=1:mesh_dst.ngl
                        for f in 1:nf
                            ip_dst = mesh_dst.connijk[k_dst,i,j]
                            ip_src = mesh_src.connijk[k_src,i,j]
                            q_dst[ip_dst, f] = 0.0
                        end
                    end

                    # Interpolation in x: (I ⊗ I ⊗ A) * v
                    for iy in 0:ny-1, ix in 0:nx-1
                        for jx in 0:nx-1
                            for f in 1:nf
                                ip_src = mesh_src.connijk[k_src,ix+1,iy+1]
                                qx[iy * nx + jx + 1, f] += interp_x[jx + 1, ix + 1, cx] * q_src[ip_src, f]
                            end
                        end
                    end

                    # Interpolation in y: (I ⊗ A ⊗ I) * v
                    for iy in 0:ny-1, jy in 0:ny-1
                        for ix in 0:nx-1
                            for f in 1:nf
                                ip_dst = mesh_dst.connijk[k_dst,ix+1,jy+1]
                                q_dst[ip_dst, f] += interp_y[jy + 1, iy + 1, cy] * qx[iy * nx + ix + 1, f]
                            end
                        end
                    end

                    # # Interpolation in z: (A ⊗ I ⊗ I) * v
                    # for iz in 0:nz-1, jz in 0:nz-1
                    #     for iy in 0:ny-1, ix in 0:nx-1
                    #         for f in 1:nf
                    #             q_dst[f, n, o_dst + (jz * ny + iy) * nx + ix + 1] += interp_z[jz + 1, iz + 1, cz] * qxy[f, (iz * ny + iy) * nx + ix + 1]
                    #         end
                    #     end
                    # end
                    k_dst += 1
                    k_src = n2o_ele_map[k_dst]
                end
                # k_src += 1
            end
        end
    else
        while k_dst <= num_elem
            k_src_vec = n2o_ele_map[k_dst]
            k_src = k_src_vec[1]
            # Source and Destination are at the same level
            if lvl_src[k_src] == lvl_dst[k_dst]
                for j=1:mesh_dst.ngl, i=1:mesh_dst.ngl
                    for f in 1:nf
                        ip_dst = mesh_dst.connijk[k_dst,i,j]
                        ip_src = mesh_src.connijk[k_src,i,j]
                        q_dst[ip_dst, f] = q_src[ip_src, f]
                    end
                end
                k_dst += 1

            # Destination is refined (hanging node)
            elseif lvl_src[k_src] < lvl_dst[k_dst]
                # Interpolate for each child
                for order in child_order
                    cy = order[2]
                    cx = order[1]
                    qx .= 0.0
                    # qxy .= 0.0

                    # Zero out q_dst for this element
                    for j=1:mesh_dst.ngl, i=1:mesh_dst.ngl
                        for f in 1:nf
                            ip_dst = mesh_dst.connijk[k_dst,i,j]
                            ip_src = mesh_src.connijk[k_src,i,j]
                            q_dst[ip_dst, f] = 0.0
                        end
                    end

                    # Interpolation in x: (I ⊗ I ⊗ A) * v
                    for iy in 0:ny-1, ix in 0:nx-1
                        for jx in 0:nx-1
                            for f in 1:nf
                                ip_src = mesh_src.connijk[k_src,ix+1,iy+1]
                                # @info k_src
                                qx[iy * nx + jx + 1, f] += interp_x[jx + 1, ix + 1, cx] * q_src[ip_src, f]
                            end
                        end
                    end

                    # Interpolation in y: (I ⊗ A ⊗ I) * v
                    for iy in 0:ny-1, jy in 0:ny-1
                        for ix in 0:nx-1
                            for f in 1:nf
                                ip_dst = mesh_dst.connijk[k_dst,ix+1,jy+1]
                                q_dst[ip_dst, f] += interp_y[jy + 1, iy + 1, cy] * qx[iy * nx + ix + 1, f]
                            end
                        end
                    end

                    # # Interpolation in z: (A ⊗ I ⊗ I) * v
                    # for iz in 0:nz-1, jz in 0:nz-1
                    #     for iy in 0:ny-1, ix in 0:nx-1
                    #         for f in 1:nf
                    #             q_dst[f, n, o_dst + (jz * ny + iy) * nx + ix + 1] += interp_z[jz + 1, iz + 1, cz] * qxy[f, (iz * ny + iy) * nx + ix + 1]
                    #         end
                    #     end
                    # end
                    k_dst += 1
                    k_src = n2o_ele_map[k_dst][1]
                end
                # k_src += 1
            else
                # Source is refined (hanging node)
                # Zero out q_dst for this element
                for j=1:mesh_dst.ngl, i=1:mesh_dst.ngl
                    for f in 1:nf
                        ip_dst = mesh_dst.connijk[k_dst,i,j]
                        q_dst[ip_dst, f] = 0.0
                    end
                end

                # Interpolation with projection in x, y, and z
                for order in 1:4
                    k_src_i = k_src_vec[order]
                    cy = child_order[order][2]
                    cx = child_order[order][1]
                    qx .= 0.0
                    # qxy .= 0.0

                    # Projection in x: (I ⊗ I ⊗ A) * v
                    for iy in 0:ny-1, ix in 0:nx-1
                        for jx in 0:nx-1
                            for f in 1:nf
                                ip_src = mesh_src.connijk[k_src_i,ix+1,iy+1]
                                qx[iy * nx + jx + 1, f] += project_x[jx + 1, ix + 1, cx] * q_src[ip_src, f]
                            end
                        end
                    end

                    # Projection in y: (I ⊗ A ⊗ I) * v
                    for iy in 0:ny-1, jy in 0:ny-1
                        for ix in 0:nx-1
                            for f in 1:nf
                                ip_dst = mesh_dst.connijk[k_dst,ix+1,jy+1]
                                q_dst[ip_dst, f] += project_y[jy + 1, iy + 1, cy] * qx[iy * nx + ix + 1, f]
                            end
                        end
                    end

                    #     # Projection in z: (A ⊗ I ⊗ I) * v
                    #     for iz in 0:nz-1, jz in 0:nz-1
                    #         for iy in 0:ny-1, ix in 0:nx-1
                    #             for f in 1:nf
                    #                 q_dst[f, n, o_dst + (jz * ny + iy) * nx + ix + 1] += project_z[jz + 1, iz + 1, cz] * qxy[f, (iz * ny + iy) * nx + ix + 1]
                    #             end
                    #         end
                    #     end

                end
                k_dst += 1
            end
        end
    end
end




const get_d_to_face_to_parent_face = Gridap.Adaptivity.get_d_to_face_to_parent_face


function mod_mesh_adaptive!(partitioned_model_coarse, ref_coarse_flags, omesh, mesh::St_mesh, inputs::Dict, nparts, distribute)

    # determine backend
    backend = CPU()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    
    #
    # Read GMSH grid from file
    #      
    parts  = distribute(LinearIndices((nparts,)))
    mesh.parts = distribute(LinearIndices((nparts,)))
    mesh.nparts = nparts
    ladpative = 1
    if ladpative == 0
        gmodel = GmshDiscreteModel(inputs[:gmsh_filename], renumber=true)
        partitioned_model = GmshDiscreteModel(parts, inputs[:gmsh_filename], renumber=true)
    elseif ladpative == 1
        # ref_coarse_flags=map(parts,partition(get_cell_gids(partitioned_model_coarse.dmodel))) do rank,indices
        #     flags=zeros(Cint,length(indices))
        #     flags.=nothing_flag
        #     if rank == 3
        #         flags .= refine_flag
        #     end
        #     flags
        # end
        partitioned_model,glue_adapt=Gridap.Adaptivity.adapt(partitioned_model_coarse,ref_coarse_flags);
    end



    cell_gids = local_views(partition(get_cell_gids(partitioned_model))).item_ref[]
    dmodel = local_views(partitioned_model.dmodel.models).item_ref[]
    model  = DiscreteModelPortion(dmodel, own_to_local(cell_gids))
    topology      = get_grid_topology(model)
    gtopology      = get_grid_topology(model)
    mesh.nsd      = num_cell_dims(model)
    
    POIN_flg = 0
    EDGE_flg = 1
    FACE_flg = 2
    ELEM_flg = 3
    
    if mesh.nsd == 3
        mesh.SD = NSD_3D()

        mesh.NNODES_EL  = 8
        mesh.NEDGES_EL  = 12
        mesh.NFACES_EL  = 6
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 4
    elseif mesh.nsd == 2

        mesh.SD = NSD_2D()
        ELEM_flg = FACE_flg
        
        mesh.NNODES_EL  = 4
        mesh.NEDGES_EL  = 4
        mesh.NFACES_EL  = 1
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 4
    elseif mesh.nsd == 1
        mesh.SD = NSD_1D()

        mesh.NNODES_EL  = 2
        mesh.NEDGES_EL  = 1
        mesh.NFACES_EL  = 0
        mesh.EDGE_NODES = 2
        mesh.FACE_NODES = 0
    else
        error( " WRONG NSD: This is not theoretical physics: we only handle 1, 2, or 3 dimensions!")
    end
    

    vtk_directory = "./refine/" 
    writevtk(partitioned_model, vtk_directory)



    p2pp = Geometry.get_face_to_parent_face(model,0)

    pgids = local_views(partition(get_face_gids(partitioned_model,POIN_flg))).item_ref[]
    edgids = local_views(partition(get_face_gids(partitioned_model,EDGE_flg))).item_ref[]
    fgids = local_views(partition(get_face_gids(partitioned_model,FACE_flg))).item_ref[]
    elgids = local_views(partition(get_face_gids(partitioned_model,ELEM_flg))).item_ref[]
    point2ppoint = local_to_global(pgids)[p2pp]
    edge2pedge = local_to_global(edgids)
    face2pface = local_to_global(fgids)
    elm2pelm = local_to_global(elgids)


    mesh.gnpoin_linear = num_faces(partitioned_model.dmodel,POIN_flg)    
    mesh.gnpoin        = mesh.gnpoin_linear         #This will be updated for the high order grid
    mesh.gnedges       = num_faces(partitioned_model,EDGE_flg)
    mesh.gnfaces       = num_faces(partitioned_model,FACE_flg)   
    mesh.gnelem        = num_faces(partitioned_model,ELEM_flg)


    mesh.npoin_linear = num_faces(model,POIN_flg)    
    mesh.npoin        = mesh.npoin_linear         #This will be updated for the high order grid
    mesh.nedges       = num_faces(model,EDGE_flg)
    mesh.nfaces       = num_faces(model,FACE_flg)   
    mesh.nelem        = num_faces(model,ELEM_flg)


    
    if (ladpative == 1)
        mesh.nelem_bdy    = length(MyGeometry.get_boundary_cells(model,mesh.nsd))
        mesh.nfaces_bdy   = length(MyGeometry.get_boundary_faces(model,mesh.nsd,FACE_flg))
        mesh.nedges_bdy   = length(MyGeometry.get_boundary_faces(model,mesh.nsd,EDGE_flg))
    else 
        mesh.nelem_bdy    = count(get_isboundary_face(topology,mesh.nsd))
        mesh.nfaces_bdy   = count(get_isboundary_face(topology,mesh.nsd-1))
        mesh.nedges_bdy   = count(get_isboundary_face(topology,mesh.nsd-2))
    end
    mesh.nelem_int    = mesh.nelem - mesh.nelem_bdy
    mesh.nfaces_int   = mesh.nfaces - mesh.nfaces_bdy
    mesh.nedges_int   = mesh.nedges - mesh.nedges_bdy

    #get_isboundary_face(topology,mesh.nsd-1)
    println(" # GMSH LINEAR GRID PROPERTIES")
    println(" # N. Global points         : ", mesh.gnpoin_linear)
    println(" # N. Global elements       : ", mesh.gnelem)
    println(" # N. Global edges          : ", mesh.gnedges)
    println(" # N. Global faces          : ", mesh.gnfaces)    
    println(" # N. points         : ", mesh.npoin_linear)
    println(" # N. elements       : ", mesh.nelem)
    println(" # N. edges          : ", mesh.nedges)
    println(" # N. faces          : ", mesh.nfaces)    
    println(" # N. internal elem  : ", mesh.nelem_int)
    println(" # N. internal edges : ", mesh.nedges_int) 
    println(" # N. internal faces : ", mesh.nfaces_int)    
    println(" # N. boundary elem  : ", mesh.nelem_bdy)
    println(" # N. boundary edges : ", mesh.nedges_bdy)
    println(" # N. boundary faces : ", mesh.nfaces_bdy)
    println(" # GMSH LINEAR GRID PROPERTIES ...................... END")
    
    ngl                     = mesh.nop + 1
    tot_linear_poin         = mesh.npoin_linear
    
    tot_edges_internal_nodes = mesh.nedges*(ngl-2)
    tot_faces_internal_nodes = mesh.nfaces*(ngl-2)*(ngl-2)
    tot_vol_internal_nodes   = mesh.nelem*(ngl-2)^(mesh.nsd)
    
    el_edges_internal_nodes = mesh.NEDGES_EL*(ngl-2)
    el_faces_internal_nodes = mesh.NFACES_EL*(ngl-2)*(ngl-2)
    el_vol_internal_nodes   = (ngl-2)^(mesh.nsd)
    
    #Update number of grid points from linear count to total high-order points
    mesh.npoin = tot_linear_poin + tot_edges_internal_nodes + tot_faces_internal_nodes + (mesh.nsd - 2)*tot_vol_internal_nodes
    
    if (mesh.nop > 1)
        println(" # GMSH HIGH-ORDER GRID PROPERTIES")
        println(" # N. edges internal points   : ", tot_edges_internal_nodes)
        println(" # N. faces internal points   : ", tot_faces_internal_nodes)
        println(" # N. volumes internal points : ", tot_vol_internal_nodes)
        println(" # N. total high order points : ", mesh.npoin)
        println(" # GMSH HIGH-ORDER GRID PROPERTIES ...................... END")
    end
    
    #
    # Resize as needed
    #
    mesh.x = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
    mesh.y = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))
    mesh.z = KernelAbstractions.zeros(backend, TFloat, Int64(mesh.npoin))

    mesh.ip2gip = KernelAbstractions.zeros(backend, TInt, Int64(mesh.npoin))
    mesh.gip2owner = KernelAbstractions.ones(backend, TInt, Int64(mesh.npoin))*local_views(parts).item_ref[]
    
    
    mesh.conn_edge_el = KernelAbstractions.zeros(backend, TInt, 2, Int64(mesh.NEDGES_EL), Int64(mesh.nelem))    
    mesh.conn_face_el = KernelAbstractions.zeros(backend, TInt,  4, Int64(mesh.NFACES_EL), Int64(mesh.nelem))  
    mesh.bdy_edge_in_elem = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy))  
    mesh.poin_in_edge = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges), Int64(mesh.ngl))
    mesh.poin_in_bdy_edge = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy), Int64(mesh.ngl))
    
    mesh.poin_in_face     = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nfaces), Int64(mesh.ngl), Int64(mesh.ngl))
    mesh.edge_type        = Array{Union{Nothing, String}}(nothing, Int64(mesh.nedges))
    mesh.bdy_edge_type    = Array{Union{Nothing, String}}(nothing, Int64(mesh.nedges_bdy))
    mesh.bdy_edge_type_id = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nedges_bdy))  
    
    if mesh.nsd > 2
        mesh.poin_in_bdy_face = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nfaces_bdy), Int64(mesh.ngl), Int64(mesh.ngl))
        mesh.face_type = Array{Union{Nothing, String}}(nothing, Int64(mesh.nfaces))
        mesh.bdy_face_type = Array{Union{Nothing, String}}(nothing, Int64(mesh.nfaces_bdy))
        mesh.bdy_face_in_elem = KernelAbstractions.zeros(backend, TInt,  Int64(mesh.nfaces_bdy))
    end
    mesh.npoin_el         = mesh.NNODES_EL + el_edges_internal_nodes + el_faces_internal_nodes + (mesh.nsd - 2)*el_vol_internal_nodes
    mesh.conn = KernelAbstractions.zeros(backend,TInt, Int64(mesh.nelem), Int64(mesh.npoin_el))
    
    #
    # Connectivity matrices
    #
    mesh.cell_node_ids     = get_cell_node_ids(get_grid(model))
    mesh.conn_unique_faces = get_face_nodes(model, FACE_flg) #faces --> 4 nodes
    mesh.conn_unique_edges = get_face_nodes(model, EDGE_flg) #edges --> 2 nodes

    mesh.cell_edge_ids     = get_faces(topology, mesh.nsd, 1) #edge map from local to global numbering i.e. iedge_g = cell_edge_ids[1:NELEM][1:NEDGES_EL]
    mesh.cell_face_ids     = get_faces(topology, mesh.nsd, mesh.nsd-1) #face map from local to global numbering i.e. iface_g = cell_face_ids[1:NELEM][1:NFACE_EL]
    mesh.face_edge_ids     = get_faces(topology,mesh.nsd-1, 1)
    mesh.edge_g_color::Array{Int64, 1} = zeros(Int64, mesh.nedges)
    # @info rank,node_global

    #
    # element refinement level
    #
    mesh.ad_lvl = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem))

    glue = local_views(glue_adapt).item_ref[]
    r_c_flags = local_views(ref_coarse_flags).item_ref[]
    for i = 1:omesh.nelem
        elem_idx = glue.o2n_faces_map[i]
        for j in elem_idx
            mesh.ad_lvl[j] = omesh.ad_lvl[i]
            if r_c_flags[i] == 1
                mesh.ad_lvl[j] += 1
            elseif r_c_flags[i] == 2
                mesh.ad_lvl[j] -= 1
            end
        end
    end
    @info rank, mesh.ad_lvl


    if (mesh.nsd == 1)
        nothing
    elseif (mesh.nsd == 2)
    
        mesh.connijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl),1)
    
        for iel = 1:mesh.nelem
            mesh.conn[iel, 1] = mesh.cell_node_ids[iel][1]
            mesh.conn[iel, 2] = mesh.cell_node_ids[iel][2]
            mesh.conn[iel, 3] = mesh.cell_node_ids[iel][4]
            mesh.conn[iel, 4] = mesh.cell_node_ids[iel][3]

            #
            # 3-----4
            # |     |
            # |     |
            # 1-----2
            #
            mesh.connijk[iel, 1,      1] = mesh.cell_node_ids[iel][2]
            mesh.connijk[iel, 1,    ngl] = mesh.cell_node_ids[iel][1]
            mesh.connijk[iel, ngl,  ngl] = mesh.cell_node_ids[iel][3]
            mesh.connijk[iel, ngl,    1] = mesh.cell_node_ids[iel][4]
            
            # @printf(" [1,1] [ngl, 1] [1, ngl] [ngl, ngl] %d %d %d %d\n", mesh.connijk[iel, 1, 1], mesh.connijk[iel, ngl, 1] , mesh.connijk[iel, 1,ngl], mesh.connijk[iel, ngl, ngl] )
        end
        #
        # Fill in elements dictionary needed by NodeOrdering.jl
        #
        elements = Dict(
            kk => mesh.conn[kk, 1:4]
            for kk = 1:mesh.nelem)
        element_types = Dict(
            kk => :Quad4
            for kk = 1:mesh.nelem)
        
        #
        # Rewrite coordinates in RCM order:
        #
        # filename = "./COORDS_LO_" + rank + ".dat" 
        open("./COORDS_LO_$rank.dat", "w") do f
            for ip = 1:mesh.npoin_linear
                
                mesh.x[ip] = get_node_coordinates(get_grid(model))[ip][1]
                mesh.y[ip] = get_node_coordinates(get_grid(model))[ip][2]
                
                mesh.ip2gip[ip] = point2ppoint[ip]
                # mesh.gip2owner[ip] = 1
                @printf(f, " %.6f %.6f 0.000000 %d %d\n", mesh.x[ip],  mesh.y[ip], ip, point2ppoint[ip])
            end
        end #f

    elseif (mesh.nsd == 3)
        
        mesh.connijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl), Int64(mesh.ngl))
        mesh.conn_edgesijk = KernelAbstractions.zeros(backend, TInt, Int64(mesh.nelem), Int64(mesh.NEDGES_EL))

        for iel = 1:mesh.nelem
            #CGNS numbering: OK ref: HEXA...
            mesh.conn[iel, 1] = mesh.cell_node_ids[iel][1]#9
            mesh.conn[iel, 2] = mesh.cell_node_ids[iel][5]#11
            mesh.conn[iel, 3] = mesh.cell_node_ids[iel][6]#6
            mesh.conn[iel, 4] = mesh.cell_node_ids[iel][2]#1
            mesh.conn[iel, 5] = mesh.cell_node_ids[iel][3]#10
            mesh.conn[iel, 6] = mesh.cell_node_ids[iel][7]#12
            mesh.conn[iel, 7] = mesh.cell_node_ids[iel][8]#5
            mesh.conn[iel, 8] = mesh.cell_node_ids[iel][4]#4

            #OK
            mesh.connijk[iel, 1, 1, 1]       = mesh.cell_node_ids[iel][2]
            mesh.connijk[iel, ngl, 1, 1]     = mesh.cell_node_ids[iel][1]
            mesh.connijk[iel, ngl, ngl, 1]   = mesh.cell_node_ids[iel][5]
            mesh.connijk[iel, 1, ngl, 1]     = mesh.cell_node_ids[iel][6]
            mesh.connijk[iel, 1, 1, ngl]     = mesh.cell_node_ids[iel][4]
            mesh.connijk[iel, ngl, 1, ngl]   = mesh.cell_node_ids[iel][3]
            mesh.connijk[iel, ngl, ngl, ngl] = mesh.cell_node_ids[iel][7]
            mesh.connijk[iel, 1, ngl, ngl]   = mesh.cell_node_ids[iel][8]
            
        end
        
        #
        # Fill in elements dictionary needed by NodeOrdering.jl
        #
        elements = Dict(
            kk => mesh.conn[kk, 1:8]
            for kk = 1:mesh.nelem)
        element_types = Dict(
            kk => :Hexa8
            for kk = 1:mesh.nelem)
        
        #
        #Use NodeNumbering.jl
        #
        #adjacency = create_adjacency_graph(elements, element_types)
        #degrees = node_degrees(adjacency)
        #neworder = RCM(adjacency, degrees, tot_linear_poin, tot_linear_poin)
        #finalorder = renumbering(neworder)
        #RCM_adjacency = create_RCM_adjacency(adjacency, finalorder)
        #newmatrix = adjacency_visualization(RCM_adjacency)
        #display(UnicodePlots.heatmap(newmatrix))
        
        
        #
        # Rewrite coordinates in RCM order:
        #
        open("./COORDS_LO_$rank.dat", "w") do f
            #open("./COORDS_LO.dat", "w") do f
            for ip = 1:mesh.npoin_linear
                mesh.x[ip] = get_node_coordinates(get_grid(model))[ip][1]
                mesh.y[ip] = get_node_coordinates(get_grid(model))[ip][2]
                mesh.z[ip] = get_node_coordinates(get_grid(model))[ip][3]
                mesh.ip2gip[ip] = point2ppoint[ip]
                # mesh.gip2owner[ip] = 1
                @printf(f, " %.6f %.6f %.6f %d %d\n", mesh.x[ip],  mesh.y[ip], mesh.z[ip], ip, point2ppoint[ip])
            end
        end #f
    end


    #
    # Add high-order points to edges, faces, and elements (volumes)
    #
    # initialize LGL struct and buyild Gauss-Lobatto-xxx points
    lgl = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.nop, backend)

    println(" # POPULATE GRID with SPECTRAL NODES ............................ ")
    #
    # Edges
    #
    populate_conn_edge_el!(mesh, mesh.SD)
    add_high_order_nodes_edges!(mesh, lgl, mesh.SD, backend, edge2pedge)

    #
    # Faces
    #
    populate_conn_face_el!(mesh, mesh.SD)
    add_high_order_nodes_faces!(mesh, lgl, mesh.SD, face2pface)

    #
    # Volume
    #
    # NOTICE: in 2D we consider only edges. faces are the elements.
    #         
    add_high_order_nodes_volumes!(mesh, lgl, mesh.SD, elm2pelm)

    mesh.gnpoin = MPI.Allreduce(maximum(mesh.ip2gip), MPI.MAX, comm)
    # @info mesh.gnpoin
    mesh.gip2owner = find_gip_owner(mesh.ip2gip)
    
    for ip = mesh.npoin_linear+1:mesh.npoin
        mesh.x[ip] = mesh.x_ho[ip]
        mesh.y[ip] = mesh.y_ho[ip]
        mesh.z[ip] = 0.0
        if (mesh.nsd > 2)
            mesh.z[ip] = mesh.z_ho[ip]
        end
    end
    
    mesh.xmax = maximum(mesh.x)
    mesh.xmin = minimum(mesh.x)
    mesh.ymax = maximum(mesh.y)
    mesh.ymin = minimum(mesh.y)
    if (mesh.nsd > 2)
        mesh.zmax = maximum(mesh.z)
        mesh.zmin = minimum(mesh.z)
    end

    for ip = 1: mesh.npoin
        if mesh.gip2owner[ip] != rank+1
            # @info mesh.x[ip], mesh.y[ip], mesh.gip2owner[ip], rank+1
        end
    end

    #----------------------------------------------------------------------
    # Extract boundary edges and faces nodes:
    #----------------------------------------------------------------------
    #
    # Bdy edges
    #
    if mesh.nsd == 2
        # isboundary_edge = compute_isboundary_face(topology, EDGE_flg)
        isboundary_edge = fill(false, mesh.nedges)  
        
        # @info isboundary_edge
        #
        # Get labels contained in the current GMSH grid:
        #
        n_semi_inf = 0
        labels = get_face_labeling(model)
        # @info rank, labels.tag_to_name
        for ilabel in labels.tag_to_name
            edges_to_tag  = get_face_tag_index(labels,ilabel,EDGE_flg)
            idx_edges_inflow = findall( x -> x == 1, edges_to_tag)
            #    
            # Tag the boundary edge with its type as defined in the user-provided GMSH file:
            #
            for idx in idx_edges_inflow
                mesh.edge_type[idx] = ilabel
                isboundary_edge[idx] = true
            end
            # @info mesh.edge_type
        end
        # @info rank, isboundary_edge
        iedge_bdy = 1
        for iedge = 1:mesh.nedges #total nedges
            if isboundary_edge[iedge] == true
                # if rank == 1
                #     @info mesh.x[mesh.poin_in_edge[iedge, 1]], mesh.y[mesh.poin_in_edge[iedge, 1]]
                # end
                for igl = 1:mesh.ngl
                    mesh.poin_in_bdy_edge[iedge_bdy, igl] = mesh.poin_in_edge[iedge, igl]
                    mesh.bdy_edge_type[iedge_bdy] = mesh.edge_type[iedge]
                end
                if (mesh.bdy_edge_type[iedge_bdy] == "Laguerre")
                    n_semi_inf += 1
                end
                iedge_bdy += 1
            end
        end
        for iel = 1:mesh.nelem
            for iedge_bdy = 1:mesh.nedges_bdy
                if issubset(mesh.poin_in_bdy_edge[iedge_bdy, :], mesh.connijk[iel, :, :])
                    mesh.bdy_edge_in_elem[iedge_bdy] = iel
                end
            end
        end
        # build mesh data structs for Laguerre semi-infinite elements
        if ("Laguerre" in mesh.bdy_edge_type)
            gr = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta],backend) 
            factorx = inputs[:xfac_laguerre]#0.1
            factory = inputs[:yfac_laguerre]#0.025
            mesh.connijk_lag = KernelAbstractions.zeros(backend, TInt, Int64(n_semi_inf), Int64(mesh.ngl), Int64(mesh.ngr),1)
            bdy_normals = zeros(n_semi_inf, 2)
            bdy_tangents = zeros(n_semi_inf, 2)
            e_iter = 1
            iter = mesh.npoin + 1
            x_new = KernelAbstractions.zeros(backend, TFloat, mesh.npoin + n_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)+mesh.ngr-1)
            y_new = KernelAbstractions.zeros(backend, TFloat, mesh.npoin + n_semi_inf*(mesh.ngl-1)*(mesh.ngr-1)+mesh.ngr-1)
            x_new[1:mesh.npoin] .= mesh.x[:]
            y_new[1:mesh.npoin] .= mesh.y[:]
            for iedge = 1:size(mesh.bdy_edge_type,1)
                if (mesh.bdy_edge_type[iedge] == "Laguerre") 
                    iel = mesh.bdy_edge_in_elem[iedge]
                    #find tangent and normal vectors to the boundary
                    ip = mesh.poin_in_bdy_edge[iedge,1]
                    ip1 = mesh.poin_in_bdy_edge[iedge,2]
                    #tangent vector 
                    x = mesh.x[ip]
                    x1 = mesh.x[ip1]
                    y = mesh.y[ip]
                    y1 = mesh.y[ip1]
                    tan = [x-x1, y-y1]
                    # deduce normal vector components
                    if (tan[2] > 1e-7)
                        x2 = 1.0
                        y2 = -x2*tan[1]/tan[2]
                    else
                        y2 = 1.0
                        x2 = -y2*tan[2]/tan[1]
                    end
                    nor = [x2,y2]
                    # generate unit versions of tangent and normal vectors
                    modu = sqrt(tan[1]^2+tan[2]^2)
                    tan = tan * (1/modu)
                    modu = sqrt(nor[1]^2+nor[2]^2)
                    nor = nor * (1/modu)
                    #make sure normal is outward facing
                    l = 1
                    m = 1
                    l1 = 1
                    m1 = 1
                    for ii=1:mesh.ngl
                        for jj=1:mesh.ngl
                            if (mesh.connijk[iel,ii,jj] == ip)
                                l=ii
                                m=jj
                            end
                            if (mesh.connijk[iel,ii,jj] == ip1)
                                l1 = ii
                                m1 = jj
                            end
                        end
                    end
                    if (l == l1)
                        ip2 = mesh.connijk[iel,3,m]
                    else
                        ip2 = mesh.connijk[iel,l,3]
                    end
                    v = [mesh.x[ip2]-x, mesh.y[ip2]-y]
                    if (dot(v,nor) > 0.0)
                        nor .= -nor
                    end
                    bdy_normals[e_iter,:] .= nor
                    bdy_tangents[e_iter,:] .= tan
                    for i=1:mesh.ngl
                        ip = mesh.poin_in_bdy_edge[iedge,i]
                        mesh.connijk_lag[e_iter,i,1] = ip
                        for j=2:mesh.ngr
			                if (inputs[:xscale]==1.0)
                                x_temp = mesh.x[ip] + nor[1]*gr.ξ[j]*factorx
                            else
                                x_temp = mesh.x[ip] + nor[1]*gr.ξ[j]*factorx/(inputs[:xscale] * 0.5)
			                end
                            if (inputs[:yscale] == 1.0)
			                    y_temp = mesh.y[ip] + nor[2]*gr.ξ[j]*factory
			                else 
                                y_temp = mesh.y[ip] + nor[2]*gr.ξ[j]*factory/(inputs[:yscale] * 0.5)
                            end
			                matched = 0
                            if (i == mesh.ngl || i == 1)
                                iter_end = 0
                                while (matched == 0 && iter_end == 0)
                                    for e_check = 1:n_semi_inf
                                        for i_check =1:mesh.ngl
                                            for j_check =1:mesh.ngr
                                                ip_check = mesh.connijk_lag[e_check,i_check,j_check]
                                                if (ip_check != 0.0 && e_check != e_iter)
                                                    if (AlmostEqual(x_temp,x_new[ip_check]) && AlmostEqual(y_temp,y_new[ip_check]))
                                                        mesh.connijk_lag[e_iter,i,j] = ip_check
                                                        matched = 1
                                                    end
                                                end
                                            end
                                        end 
                                    end
                                    iter_end = 1
                                end    
                            else
                                x_new[iter] = x_temp#mesh.x[ip] + nor[1]*gr.ξ[j]*factorx
                                y_new[iter] = y_temp#mesh.y[ip] + nor[2]*gr.ξ[j]*factory
                                mesh.connijk_lag[e_iter,i,j] = iter
                                iter += 1
                                matched = 1
                            end
                            if (matched == 0)
                                x_new[iter] = x_temp#mesh.x[ip] + nor[1]*gr.ξ[j]*factorx
                                y_new[iter] = y_temp#mesh.y[ip] + nor[2]*gr.ξ[j]*factory
                                mesh.connijk_lag[e_iter,i,j] = iter
                                iter += 1 
                            end
                            
                            #@info nor[1],nor[2],x_new[iter],y_new[iter], mesh.x[ip],mesh.y[ip]
                        end
                    end
                    e_iter += 1
                end
            end
            #@info mesh.npoin, iter - 1, mesh.ngr, n_semi_inf, e_iter - 1
            mesh.npoin_original = mesh.npoin
            mesh.npoin = iter -1
            mesh.x = x_new
            mesh.y = y_new
            mesh.z = KernelAbstractions.zeros(backend, TFloat, mesh.npoin)
            mesh.nelem_semi_inf = n_semi_inf
        end

    elseif mesh.nsd > 2
        # isboundary_face = compute_isboundary_face(topology, FACE_flg)
        isboundary_face = fill(false, mesh.nfaces)  
        #
        # Get labels contained in the current GMSH grid:
        #
        labels = get_face_labeling(model)
        for ilabel in labels.tag_to_name
            faces_to_tag  = get_face_tag_index(labels,ilabel,FACE_flg)
            idx_faces_inflow = findall( x -> x == 1, faces_to_tag)
            #    
            # Tag the boundary edge with its type as defined in the user-provided GMSH file:
            #
            for idx in idx_faces_inflow
                mesh.face_type[idx] = ilabel
                isboundary_face[idx] = true
            end
        end
        get_bdy_poin_in_face_on_edges!(mesh, @view(isboundary_face[:]), mesh.SD)
        iface_bdy = 1
        for iface in findall(x -> x == true, isboundary_face) #total nedges
            # if isboundary_face[iface] == true
            for igl = 1:mesh.ngl
                for jgl = 1:mesh.ngl
                    mesh.poin_in_bdy_face[iface_bdy, igl,jgl] = mesh.poin_in_face[iface, igl,jgl]
                    mesh.bdy_face_type[iface_bdy] = mesh.face_type[iface]
                    #@info "face point number", mesh.poin_in_face[iface,igl,jgl],iface,igl,jgl
                end
            end
            iface_bdy += 1
            # end
        end
        for iel = 1:mesh.nelem
            for iface_bdy = 1:mesh.nfaces_bdy
                if issubset(mesh.poin_in_bdy_face[iface_bdy, :,:], mesh.connijk[iel, :, :, :])
                    mesh.bdy_face_in_elem[iface_bdy] = iel
                end
            end
        end
        #=for iface =1:mesh.nfaces_bdy
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    ip = mesh.poin_in_bdy_face[iface,i,j]
                    @info "bdy points coords", mesh.x[ip],mesh.y[ip],mesh.z[ip]
                end
            end
        end=#
    end

    #----------------------------------------------------------------------
    # END Extract boundary edges and faces nodes
    #----------------------------------------------------------------------
    #
    #
    # Free memory of obsolete arrays
    #
    mesh.x_ho = zeros(1)
    mesh.y_ho = zeros(1)
    mesh.z_ho = zeros(1)
    #resize!(mesh.x_ho, 1)
    #resize!(mesh.y_ho, 1)
    #resize!(mesh.z_ho, 1)
    GC.gc()
    #

    write_vtk_grid_only(mesh.SD, mesh, "VTK_grid_refined", "./", parts, nparts)


    #show(stdout, "text/plain", mesh.conn')
    println(" # POPULATE GRID with SPECTRAL NODES ............................ DONE")
    # @info glue.n2o_cell_to_child_id
    @info rank, glue.o2n_faces_map, glue.n2o_cell_to_child_id
    return partitioned_model, glue.n2o_faces_map[mesh.nsd+1]
end

function projection_solutions(q_src, ref_coarse_flags, partitioned_model, omesh, nmesh, inputs, nparts, distribute)


    # Read gmsh grid using the GridapGmsh reader
    # @info q_src[:,4]
    partitioned_model_refined, n2o_ele_map = mod_mesh_adaptive!(partitioned_model, ref_coarse_flags, omesh, nmesh, inputs, nparts, distribute)
    q_dst = KernelAbstractions.zeros(CPU(),  TFloat, (nmesh.npoin, size(q_src, 2)))
    p8est_transfer_q!(q_dst, q_src, omesh.ad_lvl, nmesh.ad_lvl, nmesh, omesh, n2o_ele_map, nmesh.SD)
    return q_dst, partitioned_model_refined
end



function test_projection_solutions(omesh, qp, partitioned_model, inputs, nparts, distribute)


    q_src = KernelAbstractions.zeros(CPU(),  TFloat, (omesh.npoin, qp.neqs+1))
    q_src .= qp.qn .- qp.qe
    
    nmesh = St_mesh{TInt,TFloat, CPU()}(nsd=TInt(inputs[:nsd]),
    nop=TInt(inputs[:nop]),
    ngr=TInt(inputs[:nop_laguerre]+1),
    SD=NSD_1D())

    ref_coarse_flags=map(omesh.parts,partition(get_cell_gids(partitioned_model.dmodel))) do rank,indices
        flags=zeros(Cint,length(indices))
        flags.=nothing_flag
        for el_i = 1:omesh.nelem
            for i = 1:omesh.ngl
                for j = 1:omesh.ngl, k = 1:omesh.ngl
                    ip = omesh.connijk[el_i,i,j,k]
                    if abs(q_src[ip,5]) > 0.0
                        flags[el_i] = refine_flag
                    end
                end
            end
        end
        # @info "ref_coarse_flags", rank, flags
        flags
    end

    q_dst, partitioned_model_refined = projection_solutions(q_src, ref_coarse_flags, partitioned_model, omesh, nmesh, inputs, nparts, distribute)
    outvarsref = ("rho_ref", "uρ_ref", "vρ_ref", "wρ_ref", "theta_ref", "p_ref")
    write_vtk_ref(nmesh.SD, nmesh, q_dst, "adapt_initial_state_first", inputs[:output_dir]; nvar=length(qp.qn[1,:]), outvarsref=outvarsref)

    # qp_refined1 = initialize(nmesh.SD, inputs[:equations], nmesh, inputs, inputs[:output_dir], TFloat)

    


    ref_coarse_flags2=map(nmesh.parts,partition(get_cell_gids(partitioned_model_refined.dmodel))) do rank,indices
        flags=zeros(Cint,length(indices))
        flags.=nothing_flag
        for el_i = 1:nmesh.nelem
            for i = 1:nmesh.ngl
                for j = 1:omesh.ngl, k = 1:omesh.ngl
                    ip = nmesh.connijk[el_i,i,j,k]
                    if abs(q_dst[ip,5]) > 0.2
                        flags[el_i] = refine_flag
                    end
                    if (abs(q_dst[ip,5]) < 0.05 && abs(q_dst[ip,5] > 0.01))
                        flags[el_i] = coarsen_flag
                    end
                end
            end
        end
        @info "ref_coarse_flags", rank, flags
        flags
    end

    nmesh2 = St_mesh{TInt,TFloat, CPU()}(nsd=TInt(inputs[:nsd]),
    nop=TInt(inputs[:nop]),
    ngr=TInt(inputs[:nop_laguerre]+1),
    SD=NSD_1D())
    q_dst2, partitioned_model_refined2 = projection_solutions(q_dst, ref_coarse_flags2, partitioned_model_refined, nmesh, nmesh2, inputs, nparts, distribute)
    # @info n2o_ele_map2
    write_vtk_ref(nmesh2.SD, nmesh2, q_dst2, "adapt_initial_state_second", inputs[:output_dir]; nvar=length(qp.qn[1,:]), outvarsref=outvarsref)
    
    comm = MPI.COMM_WORLD

    MPI.Barrier(comm)
    @mystop("my stop at mesh.jl L135")
end