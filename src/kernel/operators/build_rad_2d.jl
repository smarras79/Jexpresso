using SparseArrays
using PETSc
function build_radiative_transfer_problem(mesh, inputs, neqs, ngl, dψ, ψ, ω, Je, dξdx, dξdy, dηdx, dηdy, nx, ny, elem_to_edge, 
        extra_mesh, QT::Inexact, SD::NSD_2D, AD::ContGal)
    comm = MPI.COMM_WORLD
    npoin = mesh.npoin
    nelem = mesh.nelem
    if (inputs[:adaptive_extra_meshes])
        extra_meshes_coords = [Vector{Float64}(undef, size(extra_mesh[e].extra_coords,1)) for i in 1:nelem]
        extra_meshes_connijk = [Array{Int}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop+1) for i in 1:nelem]
        extra_meshes_extra_Je = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop+1) for i in 1:nelem]
        extra_meshes_extra_nops = [Array{Float64}(undef, extra_mesh[e].extra_nelem) for i in 1:nelem]
        extra_meshes_extra_npoins = zeros(Int, nelem)
        extra_meshes_extra_nelems = zeros(Int, nelem)
        extra_meshes_ref_level = [Array{Int}(undef, extra_mesh[e].extra_nelem) for i in 1:nelem]
        npoin_ang_total = 0
        for e=1:nelem
            extra_meshes_coords[e] = extra_mesh[e].extra_coords
            extra_meshes_connijk[e] = extra_mesh[e].extra_connijk
            extra_meshes_extra_Je[e] = extra_mesh[e].extra_metrics.Je
            extra_meshes_extra_npoins[e] = extra_mesh[e].extra_npoin
            extra_meshes_extra_nelems[e] = extra_mesh[e].extra_nelem
            extra_meshes_extra_nops[e] = extra_mesh[e].extra_nop
            extra_meshes_ref_level[e] = extra_mesh[e].ref_level
            npoin_ang_total += mesh.ngl*mesh.ngl*extra_mesh[e].extra_npoin
        end
        connijk_spa = [Array{Int}(undef, ngl, ngl, maximum(connijk_ang[iel])) for iel = 1:nelem]
        nc_mat = adaptive_spatial_angular_numbering_2D_1D!(connijk_spa,nelem, ngl, connijk, extra_meshes_connijk_ang, extra_meshes_extra_nops, extra_meshes_extra_nelems,
                                                  extra_meshes_coords, mesh.x, mesh.y,extra_meshes_ref_level)
        @time LHS = sparse_lhs_assembly_2Dby1D_adaptive(extra_meshes_ref_level, ω, Je, mesh.connijk, extra_mesh[1].ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                    extra_meshes_extra_Je,
                                    extra_meshes_coords, extra_meshes_extra_nops, npoin_ang_total, nelem, ngl, extra_meshes_extra_nelem,
                                   dξdx, dξdy, dηdx, dηdy, extra_mesh.extra_npoin, inputs[:rad_HG_g], nc_mat)
        @time M = sparse_mass_assembly_2Dby1D(extra_meshes_ref_level, ω, Je, mesh.connijk, extra_mesh[1].ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_extra_connijk,
                                    extra_meshes_extra_Je,
                                    extra_meshes_coords, extra_meshes_extra_nops, npoin_ang_total, nelem, ngl, extra_meshes_extra_nelem,
                                   extra_meshes_extra_npoin, nc_mat)
        total_ip = size(LHS,1)
        pointwise_interaction = LHS * ones(Float64,total_ip)
        criterion = compute_adaptivity_criterion(pointwise_interaction, nelem, ngl, mesh.connijk, extra_meshes_connijk, extra_meshes_extra_nop, extra_meshes_extra_nelem, extra_meshes_coords)
        adapt_angular_grid_2Dby1D(criterion,thresholds,LHS,M,extra_meshes_ref_level,nelem,ngl,nelem_ang,nop_ang) 
        if !(maximum(ref_level) == 0)
            connijk_spa = [Array{Int}(undef, ngl, ngl, maximum(connijk_ang[iel])) for iel = 1:nelem]
            nc_mat = adaptive_spatial_angular_numbering_2D_1D!(connijk_spa,nelem, ngl, connijk, extra_meshes_connijk_ang, extra_meshes_extra_nops, extra_meshes_extra_nelems,
                                                  extra_meshes_coords, mesh.x, mesh.y,extra_meshes_ref_level)
            @time LHS = sparse_lhs_assembly_2Dby1D_adaptive(extra_meshes_ref_level, ω, Je, mesh.connijk, extra_mesh[1].ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                        extra_meshes_extra_Je,
                                        extra_meshes_coords, extra_meshes_extra_nops, npoin_ang_total, nelem, ngl, extra_meshes_extra_nelem,
                                        dξdx, dξdy, dηdx, dηdy, extra_mesh.extra_npoin, inputs[:rad_HG_g], nc_mat)
            @time M = sparse_mass_assembly_2Dby1D(extra_meshes_ref_level, ω, Je, mesh.connijk, extra_mesh[1].ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_extra_connijk,
                                        extra_meshes_extra_Je,
                                        extra_meshes_coords, extra_meshes_extra_nops, npoin_ang_total, nelem, ngl, extra_meshes_extra_nelem,
                                        extra_meshes_extra_npoin, nc_mat)
        end
    else
        npoin_ang_total = npoin*extra_mesh.extra_npoin
        @time LHS = sparse_lhs_assembly_2Dby1D(ω, Je, mesh.connijk, extra_mesh.ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh.ψ, extra_mesh.extra_connijk, 
                                    extra_mesh.extra_metrics.Je, 
                                    extra_mesh.extra_coords, extra_mesh.extra_nop, npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem,
                                   dξdx, dξdy, dηdx, dηdy, extra_mesh.extra_npoin, inputs[:rad_HG_g])
        @info "assembled LHS"
        @time M = sparse_mass_assembly_2Dby1D(ω, Je, mesh.connijk, extra_mesh.ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh.ψ, extra_mesh.extra_connijk,
                                    extra_mesh.extra_metrics.Je,
                                    extra_mesh.extra_coords, extra_mesh.extra_nop, npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem,
                                   extra_mesh.extra_npoin)
        @info "assembled Mass matrix"
        @info nnz(M), nnz(LHS), npoin_ang_total^2, nnz(M)/npoin_ang_total^2, nnz(LHS)/npoin_ang_total^2
        
        # inexact integration makes M diagonal, build the sparse inverse to save space
        ip2gip_extra, gip2owner_extra, gnpoin = setup_global_numbering_extra_dim(mesh.ip2gip, mesh.gip2owner, npoin, extra_mesh.extra_npoin, npoin_ang_total)
        Md = diag(M)
        pM = setup_assembler(SD, Md, ip2gip_extra, gip2owner_extra)
        if  pM != nothing 
            assemble_mpi!(Md,pM)
            M = Diagonal(Md)
            M = sparse(M)
        
        #    assemble_mpi!(LHS,pM)
        end

        I_vec = Vector{Int}()
        J_vec = Vector{Int}()
        V_vec = Vector{Float64}()
        max_entries = npoin_ang_total^2
        sizehint!(I_vec, Int64(round(max_entries*0.0001)))
        sizehint!(J_vec, Int64(round(max_entries*0.0001)))
        sizehint!(V_vec, Int64(round(max_entries*0.0001)))

        for ip=1:npoin_ang_total
            val = 1/M[ip,ip]
            push!(I_vec, ip)
            push!(J_vec, ip)
            push!(V_vec, val)
        end
        M_inv = sparse(I_vec, J_vec, V_vec)
        @info size(M_inv), size(LHS)    

        #@time M_inv = M \ Matrix(I, size(M)) #M\Diagonal(ones(npoin_ang_total))
        #M_inv = sparse(M_inv)
        M = nothing
        LHS_sp = nothing
        LHS_el_ang = nothing
        M_rad_el = nothing
        LHS_ang_spat = nothing
        Mass_ang_spat = nothing
        E_rad_el = nothing
        P_rad_el = nothing
        S_rad_el = nothing
        M_sp = nothing
        GC.gc()
    
        A = M_inv * LHS
        
        M_inv = nothing
        LHS = nothing
        GC.gc()
        RHS = zeros(TFloat, npoin_ang_total)
        ref = zeros(TFloat, npoin_ang_total)
    end

    @time for iel=1:nelem
        for i=1:ngl
            for j =1:ngl
                ip = mesh.connijk[iel,i,j]
                x = mesh.x[ip]
                y = mesh.y[ip]
                if (inputs[:adaptive_extra_meshes])

                else

                    for e_ext = 1:extra_mesh.extra_nelem
                        for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                            ip_ext = extra_mesh.extra_connijk[e_ext,iθ]
                            θ = extra_mesh.extra_coords[1,ip_ext]
                            ip_g = (ip-1) * extra_mesh.extra_npoin + ip_ext
                            κip = 10*exp(-((x-3/2)/3)^2)*exp(-y/2)
                            σip = 0.1*κip
                            gip = exp(-((1. / 3) * (x - (3 / 3.)))^2)#exp(-((x-3/3)/3)^2)
                            dgip = -(2. / 3^2) * (x - (3 / 3.)) * gip
                            hip = exp(-4. * (2 - y) / 2)#exp(-4*(2-y)/2)
                            dhip = (4. / 2) * hip
                            sip = exp(-((96 / (2. * π)) * (θ - (7. * π / 5.)))^2)#exp(-((96/(2*π))*(θ-7*π/5))^2)
                            uip = gip*hip*sip
                            ref[ip_g] = uip
                            propip = (cos(θ)*dgip*hip+sin(θ)*gip*dhip)*sip#cos(θ)*hip*sip*gip*((-2*x)/9 + 2/9) + sin(θ)*gip*sip*2*hip
                            if (ip in mesh.poin_in_bdy_edge)
                                applied = false
                                iedge = elem_to_edge[iel,i,j,1]
                                edge_i = elem_to_edge[iel,i,j,2]
                                
                                if (cos(θ) * nx[iedge,edge_i] + sin(θ)*ny[iedge,edge_i] >= 0.0)
                                    RHS[ip_g] = user_rad_bc(x,y,θ)#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                    A[ip_g,:] .= 0.0
                                    A[ip_g,ip_g] = 1.0
                                    applied = true
                                end
                                #=if(y == 2.0)
                                    if (sin(θ)* (-1) >= 0.0)
                                        RHS[ip_g] = user_rad_bc(x,y,θ)#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                        A[ip_g,:] .= 0.0
                                        A[ip_g,ip_g] = 1.0
                                        applied = true
                                    end
                                end
                                if(y == 0.0)
                                    if (sin(θ)*(1) >= 0.0)
                                        RHS[ip_g] = user_rad_bc(x,y,θ)#uip
                                        A[ip_g,:] .= 0.0
                                        A[ip_g,ip_g] = 1.0
                                        applied = true
                                    end
                                end
                                if(x == 3.0)
                                    if (cos(θ)*(-1) >= 0.0)
                                        RHS[ip_g] = user_rad_bc(x,y,θ)#uip
                                        A[ip_g,:] .= 0.0
                                        A[ip_g,ip_g] = 1.0
                                        applied = true
                                    end
                                end
                                if(x == 0.0)
                                    if (cos(θ)*(1) >= 0.0 )
                                        RHS[ip_g] = user_rad_bc(x,y,θ)#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                        A[ip_g,:] .= 0.0
                                        A[ip_g,ip_g] = 1.0
                                        applied = true
                                    end
                                end=#
                                if (applied == false)
                                    RHS[ip_g] = user_rhs(x,y,θ)#(-gip*hip*(user_f!(x,y,θ))*σip + κip*uip +  propip)
                                end
                                #if (y == mesh.ymax) || (x == mesh.xmin)
                                #    RHS[ip_g] = exp(-((192/(2*π))*(θ-7*π/4))^2)#uip
                                #    A[ip_g,:] .= 0.0
                                #    A[ip_g,ip_g] = 1.0
                                #else
                                #    RHS[ip_g] = uip #uip
                                #    A[ip_g,:] .= 0.0
                                #    A[ip_g,ip_g] = 1.0
                                #end
                            else
                                RHS[ip_g] = user_rhs(x,y,θ)#(-gip*hip*(user_f!(x,y,θ))*σip + κip*uip +  propip) 
                            end
                        end
                    end
                end
            end
        end
    end
    
    #x = real.(eigvals(Array(A)))
    #y = imag.(eigvals(Array(A)))
    #display(Makie.scatter(x, y, label="e-values"))

    #A_inv = inv(A)
    @info "built RHS"
    #@info RHS
    @info "solving system"
    As = sparse(A)
    A = nothing
    GC.gc()
    @info typeof(As)
    #@time solution = As \ RHS
    #=PETSc.initialize()
    A_petsc = PETSc.MatSeqAIJ(As)#PETSc.MatSeqAIJ(comm, ones(Float64, gnpoin, gnpoin))  
    b = PETSc.VecSeq(comm, zeros(Float64, gnpoin))
    x = PETSc.VecSeq(comm, zeros(Float64, gnpoin))
=#
    @time solution = solve_parallel_lsqr(ip2gip_extra, gip2owner_extra, As, RHS, gnpoin, npoin_ang_total, pM)
    #=for ip = 1:npoin_ang_total
        gip = ip2gip_extra[ip]
        for jp in nzrange(As,ip)
            gjp = ip2gip_extra[jp]
            A_petsc[gip,gjp] = As[ip,jp]
        end
        b[gip] = RHS[ip]
    end=#



  #  PETSc.assemble(A_petsc)
    #PETSc.assemble(b)
    #=ksp = PETSc.KSP(A_petsc; 
                    ksp_type="gmres",     # Solver type
                    pc_type="ilu",     # Preconditioner
                    ksp_gmres_restart= 100,
                    pc_factor_levels = 1,
                    pc_factor_shift_type = "nonzero",
                    pc_factor_shift_amount = 1e-12,
                    ksp_rtol=1e-8,        # Relative tolerance
                    ksp_max_it=npoin_ang_total,      # Max iterations
                    ksp_monitor=true)=#
    #this one actually works
   #= @time ksp = PETSc.KSP(A_petsc;
                    ksp_type="lsqr",     # Solver type
                    pc_type="none",     # Preconditioner
                    ksp_rtol=1e-14,        # Relative tolerance
                    ksp_max_it=npoin_ang_total,      # Max iterations
                    ksp_monitor=false)
    # Extract local solution
    =#
    #x = ksp\RHS#PETSc.LocalVector(x)
    @info maximum(solution), minimum(solution) 
    @info "done radiation solved"
    @info "dof", npoin_ang_total
    #=@info "absolute errors, inf, L1, L2", maximum(abs.(solution - ref)), sum(abs.(solution-ref))/npoin_ang_total, sqrt(sum((solution-ref).^2))/npoin_ang_total
    @info "relative errors, inf, L1 ,L2", maximum(abs.(solution - ref))/maximum(abs.(ref)), sum(abs.(solution-ref))/sum(abs.(ref)), sqrt(sum((solution-ref).^2))/sqrt(sum((ref).^2))
    @info "norms", maximum(abs.(solution - ref)), maximum(abs.(ref)), sum(abs.(solution-ref)), sum(abs.(ref)), sqrt(sum((solution-ref).^2)), sqrt(sum((ref).^2))
    =#A = nothing
    RHS = nothing
    GC.gc()
    @info "integrating solution and reference in angle"
    int_sol = zeros(TFloat, npoin,1)
    int_ref = zeros(TFloat, npoin,1)
    L2_err = 0.0
    L2_ref = 0.0
    for iel=1:nelem
        for i=1:ngl
            for j =1:ngl
                ip = mesh.connijk[iel,i,j]
                x = mesh.x[ip]
                y = mesh.y[ip]
                div = 1
                if (i==1 && j == 1) || (i==ngl && j==ngl) || (i==1 && j == ngl) || (i==ngl && j == 1)
                    if !(ip in  mesh.poin_in_bdy_edge)
                        div = 4
                    else
                        if !( (x == mesh.xmin && y == mesh.ymin) || (x == mesh.xmin && y == mesh.ymax) || (x == mesh.xmax && y == mesh.ymin) || (x == mesh.xmax && y == mesh.ymax))
                            div = 2
                        end
                    end
                
                elseif (i==1 || j==1 || i==ngl || j==ngl)
                    if !(ip in  mesh.poin_in_bdy_edge)
                        div = 2
                    end 
                end
                for e_ext = 1:extra_mesh.extra_nelem
                    for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                        if (iθ == 1 || iθ == extra_mesh.extra_nop[e_ext]+1)
                            div1 = div *2
                        else
                            div1 = div
                        end
                        ip_ext = extra_mesh.extra_connijk[e_ext,iθ]
                        θ = extra_mesh.extra_coords[1,ip_ext]
                        ip_g = (ip-1) * extra_mesh.extra_npoin + ip_ext
                        int_sol[ip] += solution[ip_g]*extra_mesh.extra_metrics.Je[e_ext,iθ]*extra_mesh.ωθ[iθ]/div1
                        int_ref[ip] += (ref[ip_g]-solution[ip_g])*extra_mesh.extra_metrics.Je[e_ext,iθ]*extra_mesh.ωθ[iθ]/div1
                        L2_ref += (ref[ip_g])^2*extra_mesh.extra_metrics.Je[e_ext,iθ]*extra_mesh.ωθ[iθ]*ω[i]*ω[j]*Je[iel,i,j]/div1
                        L2_err += (ref[ip_g]-solution[ip_g])^2*extra_mesh.extra_metrics.Je[e_ext,iθ]*extra_mesh.ωθ[iθ]*ω[i]*ω[j]*Je[iel,i,j]/div1
                    end
                end
            end
        end
    end
    #@info "new L2 norms", sqrt(L2_ref), sqrt(L2_err), sqrt(L2_err/L2_ref)
    #plot_triangulation(NSD_2D(), mesh, int_ref[:], "ref",  inputs[:output_dir], inputs; iout=1, nvar=1)
    plot_triangulation(NSD_2D(), mesh, int_sol[:], "sol",  inputs[:output_dir], inputs; iout=2, nvar=1)
end

function sparse_lhs_assembly_2Dby1D(ω, Je, connijk, ωθ, x, y, ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang,
                                   dξdx, dξdy, dηdx, dηdy, npoin_ang, rad_HG_g)

    max_entries = npoin_ang_total^2
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()

    sizehint!(I_vec, Int64(round(max_entries*0.0009)))
    sizehint!(J_vec, Int64(round(max_entries*0.0009)))
    sizehint!(V_vec, Int64(round(max_entries*0.0009)))
    HG, error = quadgk(v -> (1-rad_HG_g^2)/((1+rad_HG_g^2-2*rad_HG_g*cos(v))^(3/2)), 0, 2*π, rtol=1e-13, atol = 1e-13)    
    for iel=1:nelem
        for j=1:ngl
            for i=1:ngl
                ip = connijk[iel,i,j]
                ωJac = ω[i]*ω[j]*Je[iel,i,j]
                dξdx_ij = dξdx[iel,i,j]
                dξdy_ij = dξdy[iel,i,j]
                dηdx_ij = dηdx[iel,i,j]
                dηdy_ij = dηdy[iel,i,j]
                κ = user_extinction(x[ip],y[ip])
                σ = user_scattering_coef(x[ip],y[ip])
                for e_ext = 1:nelem_ang
                    for iθ = 1:nop_ang[e_ext]+1
                        ip_ext = connijk_ang[e_ext,iθ]
                        ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]
                        
                        for n=1:ngl
                            for m=1:ngl
                                jp = connijk[iel,m,n]
                                for jθ = 1:nop_ang[e_ext]+1
                                    jp_ext = connijk_ang[e_ext,jθ]
                                    extinction = κ*ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]*ψ_ang[iθ,jθ]
                                    i_angle = connijk_ang[e_ext,iθ]

                                    θ = coords_ang[1,i_angle]
                                    propagation = (ωJac*ωJac_rad)*ψ_ang[iθ,jθ]*(ψ[n,j]*dψ[m,i]*dξdx_ij*cos(θ) + ψ[m,i]*dψ[n,j]*dηdy_ij*sin(θ))
                                    intϕ = 0.0 
                                    for e_ext_scatter = 1:nelem_ang
                                        for kθ = 1:nop_ang[e_ext]+1
                                            div = 1
                                            if (kθ == nop_ang[e_ext]+1 || kθ == 1)
                                                div =2
                                            end
                                            ipθ = connijk_ang[e_ext_scatter,kθ]
                                            θ1 = coords_ang[1,ipθ]
                                            Φ = user_scattering_functions(θ,θ1,HG)
                                            ωJac_rad_scatter = ωθ[kθ]*Je_ang[e_ext_scatter,kθ]
                                            intϕ +=   ωJac_rad_scatter*Φ/div
                                        end
                                    end
                                    scattering = ψ_ang[iθ,jθ]*intϕ * ψ[i,m] * ψ[j,n] * ωJac*ωJac_rad*σ
                                    val = extinction + propagation - scattering
                                    idx_ip = (ip-1)*(npoin_ang) + ip_ext
                                    idx_jp = (jp-1)*(npoin_ang) + jp_ext
                                    if abs(val) > eps(Float64)  # Skip near-zero entries
                                        push!(I_vec, idx_ip)
                                        push!(J_vec, idx_jp)
                                        push!(V_vec, val)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return sparse(I_vec, J_vec, V_vec)
end

function sparse_lhs_assembly_2Dby1D_adaptive(ref_level, ω, Je, connijk, ωθ, x, y, ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang,
                                   dξdx, dξdy, dηdx, dηdy, npoin_ang, rad_HG_g, connijk_spa, nc_mat)

    adapted = !(maximum(ref_level) ==0)
    max_entries = npoin_ang_total^2
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()
    nc_mat = rowvals(nc_mat')
    sizehint!(I_vec, Int64(round(max_entries*0.0009)))
    sizehint!(J_vec, Int64(round(max_entries*0.0009)))
    sizehint!(V_vec, Int64(round(max_entries*0.0009)))
    HG, error = quadgk(v -> (1-rad_HG_g^2)/((1+rad_HG_g^2-2*rad_HG_g*cos(v))^(3/2)), 0, 2*π, rtol=1e-13, atol = 1e-13)
    for iel=1:nelem
        for j=1:ngl
            for i=1:ngl
                ip = connijk[iel,i,j]
                ωJac = ω[i]*ω[j]*Je[iel,i,j]
                dξdx_ij = dξdx[iel,i,j]
                dξdy_ij = dξdy[iel,i,j]
                dηdx_ij = dηdx[iel,i,j]
                dηdy_ij = dηdy[iel,i,j]
                κ = user_extinction(x[ip],y[ip])
                σ = user_scattering_coef(x[ip],y[ip])
                for e_ext = 1:nelem_ang[iel]
                    for iθ = 1:nop_ang[iel][e_ext]+1
                        ip_ext = connijk_ang[iel][e_ext,iθ]
                        ωJac_rad = ωθ[iθ]*Je_ang[iel][e_ext,iθ]

                        for n=1:ngl
                            for m=1:ngl
                                jp = connijk[iel,m,n]
                                for jθ = 1:nop_ang[iel][e_ext]+1
                                    jp_ext = connijk_ang[iel][e_ext,jθ]
                                    extinction = κ*ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]*ψ_ang[iθ,jθ]
                                    i_angle = connijk_ang[iel][e_ext,iθ]

                                    θ = coords_ang[iel][i_angle]
                                    propagation = (ωJac*ωJac_rad)*ψ_ang[iθ,jθ]*(ψ[n,j]*dψ[m,i]*dξdx_ij*cos(θ) + ψ[m,i]*dψ[n,j]*dηdy_ij*sin(θ))
                                    intϕ = 0.0
                                    for e_ext_scatter = 1:nelem_ang[iel]
                                        for kθ = 1:nop_ang[iel][e_ext]+1
                                            div = 1
                                            if (kθ == nop_ang[iel][e_ext]+1 || kθ == 1)
                                                div =2
                                            end
                                            ipθ = connijk_ang[iel][e_ext_scatter,kθ]
                                            θ1 = coords_ang[iel][ipθ]
                                            Φ = user_scattering_functions(θ,θ1,HG)
                                            ωJac_rad_scatter = ωθ[kθ]*Je_ang[iel][e_ext_scatter,kθ]
                                            intϕ +=   ωJac_rad_scatter*Φ/div
                                        end
                                    end
                                    scattering = ψ_ang[iθ,jθ]*intϕ * ψ[i,m] * ψ[j,n] * ωJac*ωJac_rad*σ
                                    val = extinction + propagation - scattering
                                    idx_ip = connijk_spa[iel,i,j,ip_ext]
                                    idx_jp = connijk_spa[iel,m,n,jp_ext]
                                    if abs(val) > eps(Float64)  # Skip near-zero entries
                                        push!(I_vec, idx_ip)
                                        push!(J_vec, idx_jp)
                                        push!(V_vec, val)
                                        if (adapted)
                                            for nc_i in nzrange(nc_mat',idx_ip)
                                                push!(I_vec, idx_ip)
                                                push!(J_vec, nc_rows[nc_i])
                                                push!(V_vec, val*nc_mat[nc_i,nc_rows[nc_i]])
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return sparse(I_vec, J_vec, V_vec)
end

function sparse_mass_assembly_2Dby1D(ω, Je, connijk, ωθ, x, y, ψ, dψ, ψ_ang, 
        connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang, npoin_ang)
  
    max_entries = npoin_ang_total^2
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()

    sizehint!(I_vec, Int64(round(max_entries*0.0001)))
    sizehint!(J_vec, Int64(round(max_entries*0.0001)))
    sizehint!(V_vec, Int64(round(max_entries*0.0001)))

    for iel=1:nelem
        for j=1:ngl
            for i=1:ngl
                ip = connijk[iel,i,j]
                ωJac = ω[i]*ω[j]*Je[iel,i,j]
                for e_ext = 1:nelem_ang
                    for iθ = 1:nop_ang[e_ext]+1
                        ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]
                        ip_ext = connijk_ang[e_ext,iθ]
                        for n=1:ngl
                            for m=1:ngl
                                jp = connijk[iel,m,n]
                                for jθ = 1:nop_ang[e_ext]+1
                                    jp_ext = connijk_ang[e_ext,jθ]

                                    val = ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]*ψ_ang[iθ,jθ]
                                    idx_ip = (ip-1)*(npoin_ang) + ip_ext
                                    idx_jp = (jp-1)*(npoin_ang) + jp_ext
                                    if abs(val) > eps(Float64)  # Skip near-zero entries
                                        push!(I_vec, idx_ip)
                                        push!(J_vec, idx_jp)
                                        push!(V_vec, val)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return sparse(I_vec, J_vec, V_vec)
end

function sparse_mass_assembly_2Dby1D_adaptive(ref_level, ω, Je, connijk, ωθ, x, y, ψ, dψ, ψ_ang,
        connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang, npoin_ang, connijk_spa,nc_mat)

    nc_rows = rowvals(nc_mat')
    max_entries = npoin_ang_total^2
    adapted = !(maximum(ref_level) ==0)
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()

    sizehint!(I_vec, Int64(round(max_entries*0.0001)))
    sizehint!(J_vec, Int64(round(max_entries*0.0001)))
    sizehint!(V_vec, Int64(round(max_entries*0.0001)))

    for iel=1:nelem
        for j=1:ngl
            for i=1:ngl
                ip = connijk[iel,i,j]
                ωJac = ω[i]*ω[j]*Je[iel,i,j]
                for e_ext = 1:nelem_ang[iel]
                    for iθ = 1:nop_ang[iel][e_ext]+1
                        ωJac_rad = ωθ[iθ]*Je_ang[iel][e_ext,iθ]
                        ip_ext = connijk_ang[iel][e_ext,iθ]
                        for n=1:ngl
                            for m=1:ngl
                                jp = connijk[iel,m,n]
                                for jθ = 1:nop_ang[iel][e_ext]+1
                                    jp_ext = connijk_ang[iel][e_ext,jθ]

                                    val = ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]*ψ_ang[iθ,jθ]
                                    idx_ip = connijk_spa[iel,i,j,ip_ext]
                                    idx_jp = connijk_spa[iel,m,n,jp_ext]
                                    if abs(val) > eps(Float64)  # Skip near-zero entries
                                        push!(I_vec, idx_ip)
                                        push!(J_vec, idx_jp)
                                        push!(V_vec, val)
                                        if (adapted) 
                                            for nc_i in nzrange(nc_mat',idx_ip)
                                                push!(I_vec, idx_ip)
                                                push!(J_vec, nc_rows[nc_i])
                                                push!(V_vec, val*nc_mat[nc_i,nc_rows[nc_i]])
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return sparse(I_vec, J_vec, V_vec)
end

function compute_adaptivity_criterion(pointwise_interaction, nelem, ngl, connijk, connijk_ang, nop_ang, nelem_ang, coords_ang)
    criterion = [Vector{Float64},(undef, nelem_ang[e]) for e=1:nelem]
    for iel=1:nelem 
        for j=1:ngl
            for i=1:ngl 
                ip = connijk[iel,i,j]
                for e_ext = 1:nelem_ang[iel]
                    for iθ = 1:nop_ang[iel][e_ext]+1
                        ip_ext = connijk_ang[iel][e_ext,iθ]
                        idx_ip = (ip-1)*(npoin_ang) + ip_ext
                        θ =  coords_ang[iel][ip_ext]
                        e_ext1 = 0
                        e_ext2 = 0
                        iθ1 = 0
                        iθ2 = 0
                        if(iθ > 1 && iθ < nop_ang[iel][e_ext]+1)
                            e_ext1 = e_ext
                            e_ext2 = e_ext
                            iθ1 = iθ-1
                            iθ2 = iθ+1
                        elseif (iθ == 1)
                            if(e_ext == 1)
                                e_ext1 = nelem_ang[iel]
                                e_ext2 = e_ext
                                iθ1 = nop_ang[iel][e_ext1]+1
                                iθ2 = iθ+1
                            else
                                e_ext1 = e_ext-1
                                e_ext2 = e_ext
                                iθ1 = nop_ang[iel][e_ext1]+1
                                iθ2 = iθ+1
                            end
                        elseif (iθ = nop_ang[iel][e_ext]+1)
                            if(e_ext == nelem_ang[iel])
                                e_ext1 = e_ext
                                e_ext2 = 1
                                iθ1 = iθ-1
                                iθ2 = 1
                            else
                                e_ext1 = e_ext
                                e_ext2 = e_ext+1
                                iθ1 = iθ-1
                                iθ2 = 1
                            end
                        end
                        ip_ext1 = connijk_ang[iel][e_ext1,iθ1]
                        ip_ext2 = connijk_ang[iel][e_ext2,iθ2]
                        idx_ip1 = (ip-1)*(npoin_ang) + ip_ext1
                        idx_ip2 = (ip-1)*(npoin_ang) + ip_ext2
                        θ1 = coords_ang[iel][ip_ext1]
                        θ2 = coords_ang[iel][ip_ext2]
                        Δθ = θ2-θ1
                        if (Δθ < 0)
                            Δθ += 2*π
                        end
                        criterion[e][e_ext] += (pointwise_interaction[idx_ip2]-pointwise_interaction[idx_ip1])/Δθ
                    end
                    criterion[e][e_ext] = criterion[e][e_ext]/(nop_ang[iel][e_ext]+1)
                end

            end
        end
    end
    return criterion
end

function adapt_angular_grid_2Dby1D(criterion,thresholds,ref_level,nelem,ngl,nelem_ang,nop_ang)
    lgl = basis_structs_ξ_ω!(LGL(), nop, backend)
    neighbors = zeros(Int,nelem,8,2)
    ang_adapted = zeros(Int, nelem)
    #loop through all spatial elements
    for iel = 1:nelem
        #loop through angular elements
        original_e_ext = nelem_ang[iel] #save original number of angular elements
        ang_adapt = false
        for e_ext = 1:original_e_ext
            #determine if angular element is to be adapted
            level= min(ref_level[iel,e_ext]+1, size(thresholds,1))

            if criterion[iel][e_ext] > thresholds[level] && level < size(thresholds,1)
                ref_level[iel,e_ext] += 1
                ang_adapted[iel] = 1
                # Make new angular elements and reconstruct extra_mesh arrays
                ang_adapt = true
                nelem_ang[iel] += 1
                npoin_ang[iel] += nop_ang[iel][e_ext]-1
                θmax = coords_ang[iel][1,connijk_ang[iel][e_ext,nop_ang[iel][e_ext]]]
                θmin = coords_ang[iel][1,connijk_ang[iel][e_ext,1]]
                θ12 = (θmax + θmin)/2
                connijk_ang_new = zeros(Int, nelem_ang[iel])
                coords_new = zeros(Float64, npoin_ang[iel]) 
                metrics = allocate_metrics(NSD_1D(), nelem_ang[iel], 0, nop_ang[iel][e_ext], TFloat, backend)
                nop_ang_new = zeros(Int,nelem_ang[iel])
                nop_ang_new .= nop_ang[iel][e_ext]
                iter = 1
                criterion_new = zeros(Int,nelem_ang[iel])
                ref_level_new = zeros(Int,nelem_ang[iel])
                #populate the elements coming before
                for e_ext1=1:e_ext-1
                    for i=1:nop_ang[iel][e_ext1]
                        connijk_ang_new[e_ext1,i] = iter
                        coords_new[iter] = coords_ang[iel][1,connijk_ang[iel][e_ext1,i]]
                        metrics.dxdξ[e_ext1, i, 1]  = metrics_ang[iel].dxdξ[e_ext1, i, 1]
                        metrics.Je[e_ext1, i, 1]  = metrics_ang[iel].Je[e_ext1, i, 1]
                        metrics.dξdx[e_ext1, i, 1]  = metrics_ang[iel].dξdx[e_ext1, i, 1]
                        criterion_new[e_ext] = criterion[iel][e_ext1]
                        ref_level_new[e_ext] = ref_level[iel][e_ext1]
                        if (iter != 1) || (e_ext==1) iter +=1 end
                    end
                end
                #populate for the new elements
                lgl = basis_structs_ξ_ω!(LGL(), nop, backend)
                for i=1:nop_ang[iel][e_ext]
                    ξ = lgl.ξ[i]
                    coords_new[iter] = θmin*(1.0-ξ)*0.5+θ12*(1.0 + ξ)*0.5
                    connijk_ang_new[e_ext,i]    = iter
                    metrics.dxdξ[e_ext, i, 1]   = (θ12-θmin)/2
                    metrics.Je[e_ext, i, 1]     = metrics.dxdξ[e_ext, i, 1]
                    metrics.dξdx[e_ext, i, 1]  = 1.0/metrics.Je[e_ext, i, 1]
                    criterion_new[e_ext] = criterion[iel][e_ext]
                    ref_level_new[e_ext] = ref_level[iel][e_ext]
                    if (iter != 1) iter +=1 end
                end
                for i=1:nop_ang[iel][e_ext+1]
                    ξ = lgl.ξ[i]
                    coords_new[iter] = θ12*(1.0-ξ)*0.5+θmax*(1.0 + ξ)*0.5
                    connijk_ang_new[e_ext+1,i] = iter
                    metrics.dxdξ[e_ext+1, i, 1]   = (θmax-θ12)/2
                    metrics.Je[e_ext+1, i, 1]     = metrics.dxdξ[e_ext+1, i, 1]
                    metrics.dξdx[e_ext+1, i, 1]  = 1.0/metrics.Je[e_ext+1, i, 1]
                    criterion_new[e_ext] = criterion[iel][e_ext]
                    ref_level_new[e_ext] = ref_level[iel][e_ext]
                    if (iter != 1) iter +=1 end
                end

                for e_ext1=e_ext+2:nelem_ang[iel]
                    for i=1:nop_ang[iel][e_ext1]
                        connijk_ang_new[e_ext1,i] = iter
                        coords_new[iter] = coords_ang[iel][1,connijk_ang[iel][e_ext1-1,i]]
                        metrics.dxdξ[e_ext1, i, 1]  = metrics_ang[iel].dxdξ[e_ext1-1, i, 1]
                        metrics.Je[e_ext1, i, 1]  = metrics_ang[iel].Je[e_ext1-1, i, 1]
                        metrics.dξdx[e_ext1, i, 1]  = metrics_ang[iel].dξdx[e_ext1-1, i, 1]
                        criterion_new[e_ext] = criterion[iel][e_ext1-1]
                        ref_level_new[e_ext] = ref_level[iel][e_ext1-1]
                        if (iter != 1) iter +=1 end
                    end 
                end

                for iel = 1:nelem
                    for i = 1:nop+1
                        for k = 1:nop+1
                            metrics.dxdξ[iel, k, 1]  = Δθe[iel]/2
                            metrics.Je[iel, k, 1]   = metrics.dxdξ[iel, k, 1]
                            metrics.dξdx[iel, k, 1] = 1.0/metrics.Je[iel, k, 1]
                        end
                    end
                end
                connijk_ang[iel] = connijk_ang_new
                coords_ang[iel] = coords_new
                metrics_ang[iel] = metrics
                nop_ang[iel] = nop_ang_new
                criterion[iel] = criterion_new
                ref_level[iel] = ref_level_new
            end 
        end
        if (ang_adapt)
            # Find and save spatial neighbords for non-conforming assembly
            #First find element corners
            xmin = 10^10
            xmax = -10^10
            ymin = 10^10
            ymax = -10^10
        
            for i=1:ngl
                for j=1:ngl
                    ip = connijk[iel,i,j]
                    if (x[ip] < xmin) xmin = x[ip] end
                    if (y[ip] < ymin) ymin = y[ip] end
                    if (x[ip] > xmax) xmax = x[ip] end
                    if (y[ip] > ymax) ymax = y[ip] end
                end
            end
            match_bdy = 0
            if (xmin == xmin_grid) match_bdy += 1 end
            if (xmax == xmax_grid) match_bdy += 1 end
            if (ymin == ymin_grid) match_bdy += 1 end
            if (ymax == ymax_grid) match_bdy += 1 end

            iter = 1
            found_neighbors = 0
        
            while (iter =< nelem && found_neighbors <8)
                #find corners for comparison
                xmin_i = 10^10
                xmax_i = -10^10
                ymin_i = 10^10
                ymax_i = -10^10
                for i=1:ngl
                    for j=1:ngl
                        ip = connijk[iter,i,j]
                        if (x[ip] < xmin_i) xmin_i = x[ip] end
                        if (y[ip] < ymin_i) ymin_i = y[ip] end
                        if (x[ip] > xmax_i) xmax_i = x[ip] end
                        if (y[ip] > ymax_i) ymax_i = y[ip] end
                    end
                end
                if (ymax == ymin_i || ymin == ymax_i || xmin == xmax_i || xmax == xmin_i)
                    found_neighbors += 1
                    neighbors[iel,found_neighbors,1] = iter
                    #check for conformity here
                    if !(adapted_ang[iel] == 0 && adapted_ang[iter] == 0) #if no angular refinement has taken place no need to check conformity
                        if (nelem_ang[iel] != nelem_ang[iter] || coords_ang[iel] != coords_ang[iter]) #non conforming
                            neighbors[iel,found_neighbors,2] = 1 # Save information that these neighbors are non conforming
                        end
                    end
                end
                if (match_bdy == 1 && found_neighbors == 5)
                    found_neighbors = 8
                elseif (match_bdy == 2 && found_neighbors == 3)
                        found_neighbors = 8
                end
                iter += 1
            end
        end
    end

end

function adaptive_spatial_angular_numbering_2D_1D!(connijk_spa,nelem, ngl, connijk, connijk_ang, nop_ang, nelem_ang, ang_coords, x, y,ref_level)
    points = []
    x_points = []
    y_points = []
    θ_points = []
    iter = 1
    interp_sources = zeros(Float64,nop_ang[1][1])
    interp_targets = zeros(Float64,nop_ang[1][1])
    ω = zeros(Float64,nop_ang[1][1])
    L = zeros(Float64,nop_ang[1][1],nop_ang[1][1])
    for iel = 1:nelem
        for i=1:ngl
            for j=1:ngl
                ip = connijk[iel,i,j]
                x_p = x[ip]
                y_p = y[ip]
                for e_ext=1:nelem_ang[iel]
                    for iθ = 1:nop_ang[iel]
                        ip_ang = connijk_ang[iel][e_ext,iθ]
                        θ_p = ang_coords[iel][1,ip_ang]
                        if (x_p in x_points && y_p in y_points && θ_p in θ_points)
                            found = false
                            iter1 = 1
                            while (found == false && iter1 < iter)
                                if (x_points[iter1] == x_p && y_points[iter1] == y_p && θ_points[iter1] = θ_p)
                                   connijk_spa[iel][i,j,ip_ang] = iter1  
                                else
                                    iter1 +=1
                                end
                            end
                        else
                            connijk_spa[iel][i,j,ip_ang] = iter
                            push!(x_points,x_p)
                            push!(y_points,y_p)
                            push!(θ_points,θ_p)
                            push!(points,iter)
                            iter += 1
                        end

                    end
                end
            end
        end
    end
    max_entries = iter^2
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()
    sizehint!(I_vec, Int64(round(max_entries*0.0001)))
    sizehint!(J_vec, Int64(round(max_entries*0.0001)))
    sizehint!(V_vec, Int64(round(max_entries*0.0001)))
    ### Done with conforming connectivity
    #handle non-conformity
    if (maximum(ref_level) == 0)
        for iel = 1:nelem
            if (1 in neighbors[iel,:,2]) #element has non-conforming neighbors
                for ineighbor = 1:8
                    if (neighbors[iel,ineighbor,2] == 1)
                        iel1 = neighbors[iel,ineighbor,1] # identify neighbor spatial element number
                        i = 0
                        j = 0
                        i1 = 0
                        j1 = 0
                        ip1 = 0
                        ip = 0
                        #loop through element edges to find matching spatial nodes
                        for igl=1:ngl
                            if (connijk[iel,1,igl] in connijk[iel1,:,:]) # a matching spatial node exists
                                ip = connijk[iel,1,igl]
                                i = 1
                                j = igl
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                            elseif (connijk[iel,igl,1] in connijk[iel1,:,:]) # a matching spatial node exists
                                ip = connijk[iel,igl,1]
                                i = igl
                                j = 1
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                            elseif (connijk[iel,ngl,igl] in connijk[iel1,:,:]) # a matching spatial node exists
                                ip = connijk[iel,ngl,igl]
                                i = ngl
                                j = igl
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                            else
                                ip = connijk[iel,igl,ngl]
                                i = igl
                                j = ngl
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                            end


                            # matched spatial node found
                            # matching spatial nodes must communicate non-conforming angular nodes
                            # loop through angular elements
                            for e_ext=1:nelem_ang[iel]
                                #angular elements are ordered
                                #check for an exact angular element match on neighboring element
                                e_check = 1
                                found = false
                                while (found == false && e_check <= nelem_ang[iel1])
                                    matched = true
                                    iter = 1
                                    while (matched)
                                        ip_ang = connijk_ang[iel][e_ext,iter]
                                        ip_ang1 = connijk_ang[iel1][e_check,iter]
                                        if !(coords_ang[iel][1,ip_ang] == coords_ang[iel1][1,ip_ang1])
                                            matched = false
                                        end
                                        iter +=1
                                    end
                                    if (matched == true)
                                        found = true
                                    end
                                    e_check += 1
                                end
                                # There is no need for nc-treatment for angular elements that are exact matches
                                if (found == false) # no exact angular element matches were found for this element, it's necessary to find NC-DSS target elements
                                    #find this element's end nodes
                                    θmin = coords_ang[1,connijk_ang[e_ext,1]]
                                    θmax = coords_ang[1,connijk_ang[e_ext,nop_ang[iel][e_ext]]]
                                    #look for target elements
                                    for e_ext1 = 1:nelem_ang[iel1]
                                        θmin1 = coords_ang[1,connijk_ang[e_ext1,1]]
                                        θmax1 = coords_ang[1,connijk_ang[e_ext,nop_ang[iel1][e_ext1]]]
                                        if (θmin1 >= θmin && θmax1 <= θmax)
                                            #First situation e_ext the "parent" element and e_ext1 is the "child" element
                                            #Find interpolating points
                                            #All parent nodes send information to all child nodes
                                            for iθ=1:nop_ang[iel][e_ext]
                                                interp_sources[iθ] = coords_ang[1,connijk_ang[e_ext,iθ]]
                                            end
                                            #Find target points
                                            resize(interp_targets,nop_ang[iel1][e_ext1])
                                            resize(L,nop_ang[iel][e_ext], nop_ang[iel1][e_ext1])
                                            for iθ=1:nop_ang[iel1][e_ext1]
                                                interp_targets[iθ] = coords_ang[1,connijk_ang[e_ext1,iθ]]
                                            end
                                            #contruct lagrange interpolator
                                            #find barycentric weights
                                            #BarycentricWeights!(iterp_sources,ω)
                                            #build interpolation matrix
                                            PolynomialInterpolationMatrix!(interp_sources,ω,interp_targets,L)
                                            #Store data for assembly
                                            for iθ=1:nop_ang[iel][e_ext]
                                                ip_ang = connijk_ang[e_ext,iθ]
                                                ip_spa = connijk_spa[iel][i,j,ip_ang]
                                                for jθ=1:nop_ang[iel1][e_ext1]
                                                    jp_ang = connijk_ang[e_ext1,jθ]
                                                    jp_spa = connijk_spa[iel][i1,j1,jp_ang]
                                                    if (ip_spa != jp_spa)
                                                        push!(I_vec, ip_ang)
                                                        push!(J_vec, jp_ang)
                                                        push!(V_vec, L[iθ,jθ])
                                                    end
                                                end
                                            end

                                        elseif (θmin >= θmin1 && θmax <= θmax1)
                                            # Second situation e_ext the "child" element and e_ext1 is the "parent" element
                                            #Find interpolating points
                                            #All child nodes send information to non-extrapolated parent nodes
                                            for iθ=1:nop_ang[iel][e_ext]
                                                interp_sources[iθ] = coords_ang[1,connijk_ang[e_ext,iθ]]
                                            end
                                            #find number of eligible target points
                                            count = 0 
                                            for iθ=1:nop_ang[iel1][e_ext1]
                                                if (coords_ang[1,connijk_ang[e_ext1,iθ]] >= θmin && coords_ang[1,connijk_ang[e_ext1,iθ]] <= θmax)
                                                    count += 1
                                                end
                                            end
                                            resize(interp_targets,count)
                                            resize(L,nop_ang[iel][e_ext], count)
                                            #Find target points
                                            count = 1
                                            for iθ=1:nop_ang[iel1][e_ext1]
                                                if (coords_ang[1,connijk_ang[e_ext1,iθ]] >= θmin && coords_ang[1,connijk_ang[e_ext1,iθ]] <= θmax)
                                                    interp_targets[count] = coords_ang[1,connijk_ang[e_ext1,iθ]]
                                                    count += 1
                                                end
                                            end
                                            #find barycentric weights
                                            BarycentricWeights!(iterp_sources,ω)
                                            #build interpolating polynomials
                                            PolynomialInterpolationMatrix!(interp_sources,ω,interp_targets,L)
                                            #Store data for assembly
                                            for iθ=1:nop_ang[iel][e_ext]
                                                ip_ang = connijk_ang[e_ext,iθ]
                                                ip_spa = connijk_spa[iel][i,j,ip_ang]
                                                count = 1
                                                for jθ=1:nop_ang[iel1][e_ext1]
                                                    if (coords_ang[1,connijk_ang[e_ext1,iθ]] >= θmin && coords_ang[1,connijk_ang[e_ext1,iθ]] <= θmax)
                                                        jp_ang = connijk_ang[e_ext1,jθ]
                                                        jp_spa = connijk_spa[iel][i1,j1,jp_ang]
                                                        if (L[iθ,count] != 1.0)
                                                            push!(I_vec, ip_ang)
                                                            push!(J_vec, jp_ang)
                                                            push!(V_vec, L[iθ,count])
                                                        end
                                                        count +=1
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return sparse(I_vec, J_vec, V_vec)
end

function find_edge_node_match(ngl,iel,ip,connijk)
    found = false
    ip1 = 0
    iter = 1
    while (found == false)
        if (iter <= 5)
            k=1
            j=iter
        elseif (iter > ngl && iter < 2*ngl)
            j = ngl
            k = iter % (ngl-1)
        elseif (iter >= 2*ngl < (3*ngl-1))
            k= ngl
            j = (ngl-1) - (iter % (ngl))
        else
            j = 1
            k = (ngl-1) - ((iter + 1) % ngl)
        end

        if (connijk[iel,k,j] == ip)
            ip1 = connijk[iel,k,j]
            found = true
        end
        iter += 1

    end
    return j, k, ip1
end
