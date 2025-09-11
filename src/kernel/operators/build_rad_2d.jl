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
        npoin_ang_total = 0
        for e=1:nelem
            extra_meshes_coords[e] = extra_mesh[e].extra_coords
            extra_meshes_connijk[e] = extra_mesh[e].extra_connijk
            extra_meshes_extra_Je[e] = extra_mesh[e].extra_metrics.Je
            extra_meshes_extra_npoins[e] = extra_mesh[e].extra_npoin
            extra_meshes_extra_nelems[e] = extra_mesh[e].extra_nelem
            extra_meshes_extra_nops[e] = extra_mesh[e].extra_nop
            npoin_ang_total += mesh.ngl*mesh.ngl*extra_mesh[e].extra_npoin
        end

        @time LHS = sparse_lhs_assembly_2Dby1D_adaptive(ω, Je, mesh.connijk, extra_mesh[1].ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                    extra_meshes_extra_Je,
                                    extra_meshes_coords, extra_meshes_extra_nops, npoin_ang_total, nelem, ngl, extra_meshes_extra_nelem,
                                   dξdx, dξdy, dηdx, dηdy, extra_mesh.extra_npoin, inputs[:rad_HG_g])
        @time M = sparse_mass_assembly_2Dby1D(ω, Je, mesh.connijk, extra_mesh[1].ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_extra_connijk,
                                    extra_meshes_extra_Je,
                                    extra_meshes_coords, extra_meshes_extra_nops, npoin_ang_total, nelem, ngl, extra_meshes_extra_nelem,
                                   extra_meshes_extra_npoin)
        total_ip = size(LHS,1)
        pointwise_interaction = LHS * ones(Float64,total_ip)
        criterion = compute_adaptivity_criterion(pointwise_interaction, nelem, ngl, mesh.connijk, extra_meshes_connijk, extra_meshes_extra_nop, extra_meshes_extra_nelem, extra_meshes_coords)
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
        @info maximum(A), minimum(A), size(A)
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
    @time solution = solve_parallel_lsqr(ip2gip_extra, gip2owner_extra, As, RHS, gnpoin, npoin_ang_total)
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

function sparse_lhs_assembly_2Dby1D_adaptive(ω, Je, connijk, ωθ, x, y, ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang,
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

function sparse_mass_assembly_2Dby1D_adaptive(ω, Je, connijk, ωθ, x, y, ψ, dψ, ψ_ang,
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

function adapt_angular_grid_2Dby1D(criterion,thresholds,LHS,M,ref_levels,nelem,ngl,nelem_ang,nop_ang)
    lgl = basis_structs_ξ_ω!(LGL(), nop, backend)
    #loop through all spatial elements
    for iel = 1:nelem
        #loop through angular elements
        original_e_ext = nelem_ang[e] #save original number of angular elements
        for e_ext = 1:original_e_ext
            #determine if angular element is to be adapted
            if criterion[iel,e_ext] > thresholds
                #adapt this angular element
                ang_ips = zeros(Int, nop_ang[iel][e_ext])
                ang_coords = zeros(Float64, nop_ang[iel][e_ext])
                ang_connijk = zeros(Float64, nelem_ang[e]+1, nop_ang[iel][e_ext])
                ang_connijk .= connijk_ang[iel][1:nelem_ang[e],:]
                ip_adapt = zeros(Int, 2, nop_ang[iel][e_ext])
                ip_taken = zeros(Int, nop_ang[iel][e_ext])
                new_coords = zeros(Float64, 2, nop_ang[iel][e_ext])
                new_ips = zeros(Int, 2, nop_ang[iel][e_ext])
                exact_node = zeros(Int, 2, nop_ang[iel][e_ext])
                non_zeros = [Vector{Int}(undef) for iθ=1:nop_ang[iel][e_ext]]
                min_idx = 100000000000000000000000
                max_idx = 0
                for iθ=1:nop_ang[iel][e_ext]
                    ip_ext = connijk_ang[iel][e_ext,iθ]
                    #first determine original angular ips to reuse and use for adaptivity
                    ang_ips[iθ] = ip_ext
                    #find original coordinates of points on the element
                    ang_coords[iθ] = coords_ang[iel][ip_ext]
                     
                end
                ω = BarycentricWeights(ang_coords)
                #split original element in two, first half keeps e_ext, second takes number  nelem_ang[e] + 1
                θmin = minimum(ang_coords)
                θmax = maximum(ang_coords)
                θhalf = (θmax+θ_min)/2
                #build LGL points for first new element and assign exact corresponding nodes
                for iθ = 1:nop_ang[iel][e_ext]
                    ξ = lgl.ξ[iθ]
                    θ = θ_min*(1.0-ξ)*0.5+θhalf*(1.0 + ξ)*0.5
                    new_coords[1,iθ] = θ
                    for jθ = 1:nop_ang[iel][e_ext]
                        if AlmostEqual(θ, ang_coords[jθ])
                            new_ips[1, iθ] = ang_ips[jθ]
                            ip_taken[jθ] = 1
                            exact_node[1,iθ] = 1
                        end
                    end

                    θ = θhalf*(1.0-ξ)*0.5+θmax*(1.0 + ξ)*0.5
                    new_coords[2,iθ] = θ
                    for jθ = 1:nop_ang[iel][e_ext]
                        if AlmostEqual(θ, ang_coords[jθ])
                            new_ips[2, iθ] = ang_ips[jθ]
                            ip_taken[jθ] = 1
                            exact_node[2,iθ] = 1
                        end
                    end
                end
                iter = npoin_ang[iel]+1
                for iθ = 1:nop_ang[iel][e_ext]
                    if (new_ips[1,iθ] == 0)
                        if 0 in ip_taken
                            for j=1:nop_ang[iel][e_ext]
                                if ip_taken[jθ] == 0 && new_ips[1,iθ] == 0
                                    new_ips[1,iθ] = ang_ips[jθ]
                                    ip_taken[jθ] = 1
                                end
                            end
                        else
                            new_ips[1,iθ] = iter
                            iter += 1
                        end
                    end
                end
                new_ips[2,1] = new_ips[1,nop_ang[iel][e_ext]]
                for iθ = 2:nop_ang[iel][e_ext]
                    if (new_ips[2,iθ] == 0)
                        if 0 in ip_taken
                            for j=1:nop_ang[iel][e_ext]
                                if ip_taken[jθ] == 0 && new_ips[2,iθ] == 0
                                    new_ips[2,iθ] = ang_ips[jθ]
                                    ip_taken[jθ] = 1
                                end
                            end
                        else
                            new_ips[2,iθ] = iter
                            iter += 1
                        end
                    end
                end
                #nodes that existed before use exact same ip, new nodes can use remaining untaken ips, when these run out ips are npoin_ang+1
                #spatial-angular numbering will use the same approach as above.
                #Do interpolations for LHS and M
                #Loop through spatial nodes since new angular grid applies to all nodes on a spatial element
                spa_ang = zeros(Int,2,nop_ang[iel][e_ext])
                for j=1:ngl
                    for i=1:ngl
                        # determine spatial ip
                        ip = connijk[iel,i,j]
                        # find non-zero columns corresponding to row of original points in matrix
                        for iθ=1:nop_ang[iel][e_ext]
                            ip_spa = (ip-1)*(npoin_ang[iel]) + ang_ips[iθ]
                            non_zeros[iθ] = nzrange(LHS',ip_spa)
                            min_idx = min(min_idx, minimum(non_zeros[iθ]))
                            max_idx = max(max_idx, maximum(non_zeros[iθ]))
                        end
                        f = zeros(Float64, max_idx - min_idx + 1, nop_ang[iel][e_ext])
                        f_int_interact = zeros(Float64, nop_ang[iel][e_ext], nop_ang[iel][e_ext])
                        fM = zeros(Float64, nop_ang[iel][e_ext])
                        #store values of LHS of original points to use in interpolation
                        for iθ = 1:nop_ang[iel][e_ext]
                            ip_spa = (ip-1)*(npoin_ang[iel]) + ang_ips[iθ]
                            for icol=min_idx:max_idx
                                f[icol, iθ] = LHS[ip_spa,icol]
                            end
                            fM[iθ] = M[ip_spa,ip_spa]
                            for jθ = 1:nop_ang[iel][e_ext]
                                jp_spa = (ip-1)*(npoin_ang[iel]) + ang_ips[jθ]
                                f_int_interact[iθ,jθ] = LHS[ip_spa, jp_spa]
                            end

                        end
                        #determine spatial angular numbering for new nodes
                        spa_iter = 1
                        for iθ = 1:nop_ang[iel][e_ext]
                            if new_ips[1,iθ] <= npoin_ang[iel]
                                spa_ang[1,iθ] = (ip-1)*(npoin_ang[iel]) + new_ips[1,iθ]
                            else
                                spa_ang[1,iθ] = npoin_ang_total+spa_iter
                                spa_iter += 1
                            end

                            if new_ips[2,iθ] <= npoin_ang[iel]
                                spa_ang[2,iθ] = (ip-1)*(npoin_ang[iel]) + new_ips[2,iθ]
                            else 
                                spa_ang[2,iθ] = npoin_ang_total+spa_iter
                                spa_iter += 1
                            end
                            
                            #Do LHS interpolation
                            #If node already existed, do nothing
                            if exact_node[1,iθ] != 1
                                #Otherwise interpolate
                                #First check if this node number already existed (no new rows or columns added to the matrix)
                                if new_ips[1,iθ] <= npoin_ang[iel]
                                    #spatial angular number corresponds to matrix row
                                    row_number = spa_ang[1,iθ]
                                    #For LHS_ij to be non-zero use only non-zero-columns corresponding to original points
                                    for icol=min_idx:max_idx
                                        new_vals = LagrangeInterpolation(new_coords[1,:],ang_coords,f[icol,:], ω)
                                    end

                                else
                                    #spatial angular number corresponds to matrix row
                                    row_number = spa_ang[1,iθ]
                                    #For LHS_ij to be non-zero use only non-zero-columns corresponding to original points
                                    for icol=min_idx:max_idx
                                        new_vals = LagrangeInterpolation(new_coords[1,iθ],ang_coords,f[icol,:], ω)
                                    end
                                    #Since this node appended at the bottom corner of the matrix (new row) make it contributes to all applicable old rows in this new column
                                    #First determine what values would have corresponded to the interpolation points on the original grid
                                    for jθ = 1:nop_ang[iel][e_ext]
                                        col_val[jθ] = LagrangeInterpolation(new_coords[1,iθ],ang_coords,f_int_interact[jθ,:], ω)
                                    end
                                    #interpolate onto previous rows
                                    for jθ = 1:iθ
                                        rowi = spa_ang[1,jθ] 
                                        new_val = LagrangeInterpolation(new_coords[1,jθ],ang_coords,col_val, ω)
                                        rowi = spa_ang[2,jθ]
                                        new_val = LagrangeInterpolation(new_coords[2,jθ],ang_coords,col_val, ω)
                                    end
                                end
                            end
                            #Same for second half of old element
                            if exact_node[2,iθ] != 1
                                #First check if this node number already existed (no new rows or columns added to the matrix)
                                if new_ips[2,iθ] <= npoin_ang[iel]
                                    #spatial angular number corresponds to matrix row
                                    row_number = spa_ang[2,iθ]
                                    #For LHS_ij to be non-zero use only non-zero-columns corresponding to original points 
                                    for icol=min_idx:max_idx
                                        new_vals = LagrangeInterpolation(new_coords[2,:],ang_coords,f[icol,:], ω)
                                    end
                                else
                                    #spatial angular number corresponds to matrix row
                                    row_number = spa_ang[2,iθ]
                                    #For LHS_ij to be non-zero use only non-zero-columns corresponding to original points
                                    for icol=min_idx:max_idx
                                        new_vals = LagrangeInterpolation(new_coords[2,:],ang_coords,f[icol,:], ω)
                                    end
                                    #Since this node appended at the bottom corner of the matrix (new row) make it contributes to all applicable old row in this new column
                                    #First determine what values would have corresponded to the interpolation points on the original grid
                                    for jθ = 1:nop_ang[iel][e_ext]
                                        col_val[jθ] = LagrangeInterpolation(new_coords[2,iθ],ang_coords,f_int_interact[jθ,:], ω)
                                    end 
                                    #interpolate onto previous rows
                                    for jθ = 1:iθ
                                        rowi = spa_ang[1,jθ]
                                        new_val = LagrangeInterpolation(new_coords[1,jθ],ang_coords,col_val, ω)
                                        rowi = spa_ang[2,jθ]
                                        new_val = LagrangeInterpolation(new_coords[2,jθ],ang_coords,col_val, ω)
                                    end 
                                end
                            end

                        end


                #update extra_mesh where still necessary
            end 
        end

    end

end
