using SparseArrays
function build_radiative_transfer_problem(mesh, inputs, neqs, ngl, dψ, ψ, ω, Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
        nx, ny, nz, elem_to_face,
        extra_mesh, QT::Inexact, SD::NSD_3D, AD::ContGal)

    npoin = mesh.npoin
    nelem = mesh.nelem
    #LHS_ang_spat = Array{Array,3}(undef, Int64(nelem), Int64(ngl^2), Int64(ngl^2))
    #Mass_ang_spat = Array{Array,3}(undef, Int64(nelem), Int64(ngl^2), Int64(ngl^2))
    #=if (inputs[:adaptive_extra_meshes])
        LHS_ang_spat = Array{Array,3}(undef, Int64(nelem), Int64(ngl^2), Int64(ngl^2))
        Mass_ang_spat = Array{Array,3}(undef, Int64(nelem), Int64(ngl^2), Int64(ngl^2))
        npoin_ang_total = 0
        for iel=1:nelem
            npoin_ang_total += extra_mesh[iel].extra_npoin*ngl*ngl
            nelem_ext = extra_mesh[iel].extra_nelem
            ngl_ext = extra_mesh[iel].extra_nop[1]+1
            κ = zeros(TFloat, ngl, ngl)            
            σ = zeros(TFloat, ngl, ngl)
            Φ = zeros(TFloat, extra_mesh[iel].extra_npoin, extra_mesh[iel].extra_npoin)
            E_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh[iel].extra_nop[1]+1, extra_mesh[iel].extra_nop[1]+1)
            P_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh[iel].extra_nop[1]+1, extra_mesh[iel].extra_nop[1]+1)
            S_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh[iel].extra_nop[1]+1, extra_mesh[iel].extra_nop[1]+1)
            user_extinction!(κ, mesh.x, mesh.y, mesh.connijk, iel, extra_mesh, ngl)
            user_scattering!(σ, Φ, mesh.x, mesh.y, mesh.connijk, iel, extra_mesh, ngl)
            _rad_element_extinction_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, E_rad_el, iel, κ, extra_mesh, QT, SD, AD)
            _rad_element_internal_propagation_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, dξdx, dξdy, dηdx, dηdy, P_rad_el, iel, extra_mesh, QT, SD, AD)
            _rad_element_scattering_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, S_rad_el, iel, σ, Φ, extra_mesh, QT, SD, AD)
            for i=1:ngl
                for j =1:ngl
                    LHS_ang_spat[iel,i,j] = zeros(TFloat, extra_mesh[iel].extra_npoin, extra_mesh[iel].extra_npoin)
                end
            end
        end
        LHS_el_ang = E_rad_el + P_rad_el - S_rad_el 
        for iel=1:nelem
            DSS_angular!(LHS_el_ang, LHS_ang_spat, nelem, extra_mesh, ngl, mesh.connijk, iel, SD)
        end

    else
        LHS_ang_spat = zeros(Float64, Int64(nelem), Int64(ngl^2), Int64(ngl^2), extra_mesh.extra_npoin, extra_mesh.extra_npoin)
        Mass_ang_spat = zeros(Float64, Int64(nelem), Int64(ngl^2), Int64(ngl^2), extra_mesh.extra_npoin, extra_mesh.extra_npoin)
        npoin_ang_total = npoin*extra_mesh.extra_npoin
        nelem_ext = extra_mesh.extra_nelem
        E_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh.extra_nop[1]+1, extra_mesh.extra_nop[1]+1)
        P_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh.extra_nop[1]+1, extra_mesh.extra_nop[1]+1)
        S_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh.extra_nop[1]+1, extra_mesh.extra_nop[1]+1)
        M_rad_el = zeros(TFloat, nelem, ngl^2, ngl^2, nelem_ext, extra_mesh.extra_nop[1]+1, extra_mesh.extra_nop[1]+1)
        κ = zeros(TFloat, ngl, ngl)
        σ = zeros(TFloat, ngl, ngl)
        Φ = zeros(TFloat, extra_mesh.extra_npoin, extra_mesh.extra_npoin)
        for iel=1:nelem
            fill!(κ, zero(Float64))
            fill!(σ, zero(Float64))
            fill!(Φ, zero(Float64))
            user_extinction!(κ, mesh.x, mesh.y, mesh.connijk, iel, extra_mesh, ngl)
            
            user_scattering!(σ, Φ, mesh.x, mesh.y, mesh.connijk, iel, extra_mesh.extra_nelem, extra_mesh.extra_nop, extra_mesh.extra_connijk,
                             extra_mesh.extra_coords, ngl)
            
            _rad_element_extinction_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, E_rad_el, iel, κ, 
                                                   extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_metrics.Je, extra_mesh.ωθ, 
                                                   extra_mesh.ψ, QT, SD, AD)
            
            _rad_element_internal_propagation_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, dξdx, dξdy, dηdx, dηdy, 
                                                             P_rad_el, iel, 
                                                             extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_metrics.Je, extra_mesh.ωθ,
                                                             extra_mesh.ψ, extra_mesh.extra_connijk, extra_mesh.extra_coords, QT, SD, AD)
            
            _rad_element_scattering_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, S_rad_el, iel, σ, Φ, 
                                                   extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_metrics.Je, extra_mesh.ωθ,
                                                   extra_mesh.ψ, extra_mesh.extra_connijk, QT, SD, AD)
            
            _rad_element_mass_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, M_rad_el, iel, 
                                             extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_metrics.Je, extra_mesh.ωθ,
                                                   extra_mesh.ψ, QT, SD, AD)
        end
        
        LHS_el_ang = E_rad_el -S_rad_el + P_rad_el #- S_rad_el
        for iel=1:nelem
            DSS_angular!(LHS_el_ang, LHS_ang_spat, nelem, extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_connijk, ngl, mesh.connijk, iel, SD)
            DSS_angular!(M_rad_el, Mass_ang_spat, nelem, extra_mesh.extra_nop, extra_mesh.extra_nelem, extra_mesh.extra_connijk, ngl, mesh.connijk, iel, SD)
        end
        LHS_sp = zeros(TFloat, npoin, npoin, extra_mesh.extra_npoin, extra_mesh.extra_npoin)
        M_sp = zeros(TFloat, npoin, npoin, extra_mesh.extra_npoin, extra_mesh.extra_npoin)
        DSS_spatial!(LHS_sp, LHS_ang_spat, nelem, extra_mesh.extra_npoin, ngl, mesh.connijk, inputs, SD)
        DSS_spatial!(M_sp, Mass_ang_spat, nelem, extra_mesh.extra_npoin, ngl, mesh.connijk, inputs, SD)
        #@time LHS_sp = DSS_spatial(LHS_ang_spat, npoin, nelem, extra_mesh.extra_npoin, ngl, mesh.connijk, inputs, SD)
        #@time M_sp = DSS_spatial(Mass_ang_spat, npoin, nelem, extra_mesh.extra_npoin, ngl, mesh.connijk, inputs, SD)
    end=#
    #LHS = zeros(TFloat, npoin_ang_total, npoin_ang_total)
    #M = zeros(TFloat, npoin_ang_total, npoin_ang_total)
    #@time LHS = Map_rad_global(LHS_sp, npoin_ang_total, npoin, extra_mesh.extra_npoin, nelem, mesh.connijk, ngl, inputs, SD)
    #@info "assembled LHS"
    #@time M = Map_rad_global(M_sp, npoin_ang_total, npoin, extra_mesh.extra_npoin, nelem, mesh.connijk, ngl, inputs, SD)
    #@info "assembled Mass matrix"
    #M = sparse(M)
    @info extra_mesh.extra_coords[1,:]
    @info extra_mesh.extra_coords[2,:]
    npoin_ang_total = npoin*extra_mesh.extra_npoin
    @info npoin_ang_total, extra_mesh.extra_npoin, npoin
    @time LHS = sparse_lhs_assembly_3Dby2D(ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ, 
                                           mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh.ψ, extra_mesh.extra_connijk, 
                                    extra_mesh.extra_metrics.Je, 
                                    extra_mesh.extra_coords, extra_mesh.extra_nop, npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem,
                                    dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, extra_mesh.extra_npoin, inputs[:rad_HG_g])
    @info "assembled LHS"
    @time M = sparse_mass_assembly_3Dby2D(ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ, mesh.x, mesh.y, ψ, dψ, extra_mesh.ψ, extra_mesh.extra_connijk,
                                    extra_mesh.extra_metrics.Je,
                                    extra_mesh.extra_coords, extra_mesh.extra_nop, npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem,
                                   extra_mesh.extra_npoin)
    @info "assembled Mass matrix"
    @info nnz(M), nnz(LHS), npoin_ang_total^2, nnz(M)/npoin_ang_total^2, nnz(LHS)/npoin_ang_total^2
    # inexact integration makes M diagonal, build the sparse inverse to save space
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
    @info maximum(A), minimum(A)
    M_inv = nothing
    LHS = nothing
    GC.gc()
    RHS = zeros(TFloat, npoin_ang_total,1)
    ref = zeros(TFloat, npoin_ang_total,1)
        
    @time for iel=1:nelem
        for i = 1:ngl
            for j = 1:ngl
                for k = 1:ngl
                    ip = mesh.connijk[iel,i,j,k]
                    x = mesh.x[ip]
                    y = mesh.y[ip]
                    z = mesh.z[ip]
                    if (inputs[:adaptive_extra_meshes])

                    else

                        for e_ext = 1:extra_mesh.extra_nelem
                            for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                                for iϕ = 1:extra_mesh.extra_nop[e_ext]+1
                                    ip_ext = extra_mesh.extra_connijk[e_ext,iϕ,iθ]
                                    θ = extra_mesh.extra_coords[1,ip_ext]
                                    ϕ = extra_mesh.extra_coords[2,ip_ext]
                                    ip_g = (ip-1) * extra_mesh.extra_npoin + ip_ext
                                    κip = 10*exp(-((x-3/2)/3)^2)*exp(-y/2)
                                    σip = 0.1*κip
                                    gip = x^2#exp(-((1. / 3) * (x - (3 / 3.)))^2)#exp(-((x-3/3)/3)^2)
                                    dgip = 2*x#-(2. / 3^2) * (x - (3 / 3.)) * gip
                                    hip = y^3#exp(-4. * (2 - y) / 2)#exp(-4*(2-y)/2)
                                    fip = z^4
                                    dfip = 4*z^3
                                    dhip = 3*y^2#(4. / 2) * hip
                                    sip = cos(θ)#exp(-((96 / (2. * π)) * (θ - (7. * π / 5.)))^2)#exp(-((96/(2*π))*(θ-7*π/5))^2)
                                    bip = sin(ϕ)
                                    uip = gip*hip*fip*sip*bip
                                    ref[ip_g] = uip
                                    propip = 0.0#(cos(θ)*dgip*hip+sin(θ)*gip*dhip)*sip#cos(θ)*hip*sip*gip*((-2*x)/9 + 2/9) + sin(θ)*gip*sip*2*hip
                                    if (ip in mesh.poin_in_bdy_face)
                                        
                                        applied = false
                                        iface = elem_to_face[iel,i,j,k,1]
                                        face_i = elem_to_face[iel,i,j,k,2]
                                        face_j = elem_to_face[iel,i,j,k,3]
                                        prodx = nx[iface,face_i,face_j]*sin(θ)*cos(ϕ)
                                        prody = ny[iface,face_i,face_j]*sin(θ)*sin(ϕ)
                                        prodz = nz[iface,face_i,face_j]*cos(θ)
                                        xb = (x == 3.0 || x == 0.0)
                                        yb = (y == 3.0 || y == 0.0)
                                        zb = (z == 2.0 || z == 0.0)
                                        #=if (xb && yb)
                                                if (y == 0.0)
                                                    ny[iface,face_i,face_j] = 1.0
                                                else
                                                    ny[iface,face_i,face_j] = -1.0
                                                end
                                                if (x == 0.0)
                                                    nx[iface,face_i,face_j] = 1.0
                                                else
                                                    nx[iface,face_i,face_j] = -1.0
                                                end
                                            prodx = nx[iface,face_i,face_j]*sin(θ)*cos(ϕ)
                                            prody = ny[iface,face_i,face_j]*sin(θ)*sin(ϕ)
                                        end
                                        if (xb && zb)
                                                if (z == 0.0)
                                                    nz[iface,face_i,face_j] = 1.0
                                                else
                                                    nz[iface,face_i,face_j] = -1.0
                                                end
                                                if (x == 0.0)
                                                    nx[iface,face_i,face_j] = 1.0
                                                else
                                                    nx[iface,face_i,face_j] = -1.0
                                                end
                                            prodx = nx[iface,face_i,face_j]*sin(θ)*cos(ϕ)
                                            prodz = nz[iface,face_i,face_j]*cos(θ)
                                        end
                                        if (yb && zb)
                                                if (z == 0.0)
                                                    nz[iface,face_i,face_j] = 1.0
                                                else
                                                    nz[iface,face_i,face_j] = -1.0
                                                end
                                                if (y == 0.0)
                                                    ny[iface,face_i,face_j] = 1.0
                                                else
                                                    ny[iface,face_i,face_j] = -1.0
                                                end
                                            prody = ny[iface,face_i,face_j]*sin(θ)*sin(ϕ)
                                            prodz = nz[iface,face_i,face_j]*cos(θ)
                                        end=#

                                        #@info nx[iface,face_i,face_j], ny[iface,face_i,face_j], nz[iface,face_i,face_j], x, y, z
                                        if (prodx + prody + prodz >= -1e-13)
                                            RHS[ip_g] = user_rad_bc(x,y,z,θ,ϕ)#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                            A[ip_g,:] .= 0.0
                                            A[ip_g,ip_g] = 1.0
                                            applied = true
                                        end
                                        if (applied == false)
                                            RHS[ip_g] = user_rhs_sphere(x,y,z,θ,ϕ)#(-gip*hip*(user_f!(x,y,θ))*σip + κip*uip +  propip)
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
                                        RHS[ip_g] = user_rhs_sphere(x,y,z,θ,ϕ)#(-gip*hip*(user_f!(x,y,θ))*σip + κip*uip +  propip) 
                                    end
                                end
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
    @time solution = As \ RHS
    
    @info "done radiation solved"
    @info "dof", npoin_ang_total
    @info "absolute errors, inf, L1, L2", maximum(abs.(solution - ref)), sum(abs.(solution-ref))/npoin_ang_total, sqrt(sum((solution-ref).^2))/npoin_ang_total
    @info "relative errors, inf, L1 ,L2", maximum(abs.(solution - ref))/maximum(abs.(ref)), sum(abs.(solution-ref))/sum(abs.(ref)), sqrt(sum((solution-ref).^2))/sqrt(sum((ref).^2))
    @info "norms", maximum(abs.(solution - ref)), maximum(abs.(ref)), sum(abs.(solution-ref)), sum(abs.(ref)), sqrt(sum((solution-ref).^2)), sqrt(sum((ref).^2))
    A = nothing
    RHS = nothing
    GC.gc()
    @info "integrating solution and reference in angle"
    int_sol = zeros(TFloat, npoin,1)
    int_ref = zeros(TFloat, npoin,1)
    L2_err = 0.0
    L2_ref = 0.0
    
    #=for iel=1:nelem
        for i = 1:ngl
            for j = 1:ngl
                for k = 1:ngl
                    ip = mesh.connijk[iel,i,j,k]
                    x = mesh.x[ip]
                    y = mesh.y[ip]
                    z = mesh.z[ip]
                    for e_ext = 1:extra_mesh.extra_nelem
                        for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                            for iϕ = 1:extra_mesh.extra_nop[e_ext]+1
                                ip_ext = extra_mesh.extra_connijk[e_ext,iϕ,iθ]
                                θ = extra_mesh.extra_coords[1,ip_ext]
                                ϕ = extra_mesh.extra_coords[2,ip_ext]
                                ip_g = (ip-1) * extra_mesh.extra_npoin + ip_ext
                                if (abs(solution[ip_g] - 27) > 1e-5)
                                    @info solution[ip_g], x, y, z, θ, ϕ
                                end
                            end
                        end
                    end
                end
            end
        end
    end=#
    #=for iel=1:nelem
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
    @info "new L2 norms", sqrt(L2_ref), sqrt(L2_err), sqrt(L2_err/L2_ref)
    plot_triangulation(NSD_2D(), mesh, int_ref[:], "ref",  inputs[:output_dir], inputs; iout=1, nvar=1)
    plot_triangulation(NSD_2D(), mesh, int_sol[:], "sol",  inputs[:output_dir], inputs; iout=2, nvar=1)
=#
end

function _rad_element_mass_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, M_rad_el, iel, 
        nop_ang, nelem_ang, Je_ang, ωθ, ψ_ang, QT::Inexact, SD::NSD_2D, AD::ContGal)
    
    for j=1:ngl
        for i=1:ngl
            ωJac = ω[i]*ω[j]*Je[iel,i,j]
            J = i +(j-1)*(ngl) 
            for e_ext = 1:nelem_ang
                for iθ = 1:nop_ang[e_ext]+1
                    J_ext = iθ
                    ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]
                            
                    for n=1:ngl
                        for m=1:ngl
                            I = m + (n-1)*(ngl)
                            for jθ = 1:nop_ang[e_ext]+1
                                I_ext = jθ
                                M_rad_el[iel,I,J,e_ext,I_ext,J_ext] += ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]*ψ_ang[iθ,jθ]
                            end
                        end
                    end
                end
            end
        end
    end
end

function _rad_element_extinction_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, E_rad_el, iel, κ, nop_ang, nelem_ang, Je_ang, ωθ, ψ_ang, QT::Inexact, SD::NSD_2D, AD::ContGal)

    for j=1:ngl
        for i=1:ngl
            ωJac = ω[i]*ω[j]*Je[iel,i,j]
            J = i +(j-1)*(ngl)
            for e_ext = 1:nelem_ang
                for iθ = 1:nop_ang[e_ext]+1
                    J_ext = iθ 
                    ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]

                    for n=1:ngl
                        for m=1:ngl
                            #ωJac = ω[m]*ω[n]*Je[iel,m,n]
                            I = m + (n-1)*(ngl)
                            for jθ = 1:nop_ang[e_ext]+1
                                I_ext = jθ
                                E_rad_el[iel,I,J,e_ext,I_ext,J_ext] += κ[i,j]*ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]*ψ_ang[iθ,jθ]
                            end
                        end
                    end
                end
            end
        end
    end
end

function _rad_element_internal_propagation_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, dξdx, dξdy, dηdx, dηdy, P_rad_el, iel,
        nop_ang, nelem_ang, Je_ang, ωθ, ψ_ang, connijk_ang, coords_ang, QT::Inexact, SD::NSD_2D, AD::ContGal)

    for j=1:ngl
        for i=1:ngl
            ωJac = ω[i]*ω[j]*Je[iel,i,j]
            I = i +(j-1)*(ngl)
            dξdx_ij = dξdx[iel,i,j]
            dξdy_ij = dξdy[iel,i,j]
            dηdx_ij = dηdx[iel,i,j]
            dηdy_ij = dηdy[iel,i,j]

            for e_ext = 1:nelem_ang
                for iθ = 1:nop_ang[e_ext]+1
                    I_ext = iθ
                    ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]
                    #=ψdψdξ = 0.0
                    ψdψdη = 0.0 
                    for n=1:ngl
                        for m=1:ngl
                            ψdψdξ += ψ[n,j]*dψ[m,i]
                            ψdψdη += ψ[m,i]*dψ[n,j]
                            @info m, n, ψdψdξ, ψdψdη, ψ[n,j]*dψ[m,i], ψ[m,i]*dψ[n,j]
                        end
                    end
                    ddx = (ψdψdξ*dξdx_ij + ψdψdη * dηdx_ij)
                    ddy = (ψdψdξ*dξdy_ij + ψdψdη * dηdy_ij)
                    =#
                    for n=1:ngl
                        for m=1:ngl
                            #ωJac = ω[m]*ω[n]*Je[iel,m,n]
                            J = m + (n-1)*(ngl)
                            for jθ = 1:nop_ang[e_ext]+1
                                J_ext = jθ
                                i_angle = connijk_ang[e_ext,iθ]

                                θ = coords_ang[1,i_angle]
                                P_rad_el[iel,I,J,e_ext,I_ext,J_ext] += (ωJac*ωJac_rad)*ψ_ang[iθ,jθ]*(ψ[n,j]*dψ[m,i]*dξdx_ij*cos(θ) + ψ[m,i]*dψ[n,j]*dηdy_ij*sin(θ))
                            end
                            #=for jθ = 1:nop_ang[e_ext]+1
                                I_ext = jθ
                                i_angle = connijk_ang[e_ext,iθ]

                                θ = coords_ang[1,i_angle]
                                #s = [cos(θ), sin(θ)]
                                for iq = 1:ngl
                                    ψdψdξ += ψ[j,n]*dψ[iq,i]
                                    ψdψdη += ψ[i,m]*dψ[iq,j]
                                end
                                ddx = (ψdψdξ*dξdx_ij + ψdψdη * dηdx_ij)*cos(θ)
                                ddy = (ψdψdξ*dξdy_ij + ψdψdη * dηdy_ij)*sin(θ)
                                
                                P_rad_el[iel,I,J,e_ext,I_ext,J_ext] += (ωJac*ωJac_rad)*(ddx + ddy)*(ψ_ang[iθ,jθ])
                                #P_rad_el[iel,I,J,e_ext,I_ext,J_ext] += (ωJac*ωJac_rad)*(ψ[n,j]*dψ[m,i]*cos(θ)*dξdx_ij+ψ[m,i]*dψ[n,j]*sin(θ)*dηdy_ij)*(ψ_ang[iθ,jθ])
                                #=if (ψ[i,m]*ψ[j,n]*extra_mesh.ψ[iθ,jθ] != 0)
                                    @info P_rad_el[iel,I,J,e_ext,I_ext,J_ext], dψ[k,n]*s[2]*dηdy_ij+dψ[k,m]*s[1]*dξdx_ij, ψ[i,m]*ψ[j,n]*extra_mesh.ψ[iθ,jθ], ωJac*ωJac_rad
                                end=#
                            end=#
                        end
                    end
                end
            end
        end
    end
end


function _rad_element_scattering_matrix_2Dby1D!(neqs, ngl, dψ, ψ, ω, Je, S_rad_el, iel, σ, Φ, nop_ang, nelem_ang, Je_ang, ωθ, ψ_ang, connijk_ang, QT::Inexact, SD::NSD_2D, AD::ContGal)

    for j=1:ngl
        for i=1:ngl
            ωJac = ω[i]*ω[j]*Je[iel,i,j]
            I = i +(j-1)*(ngl)
            for e_ext = 1:nelem_ang
                for iθ = 1:nop_ang[e_ext]+1
                    I_ext = iθ 
                    ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]
                    jpθ = connijk_ang[e_ext,iθ]
                    for n=1:ngl
                        for m=1:ngl
                            J = m + (n-1)*(ngl)
                            for jθ = 1:nop_ang[e_ext]+1
                                J_ext = jθ
                                intϕ = 0.0
                                #ωJac_rad = ωθ[jθ]*Je_ang[e_ext,jθ]
                                for e_ext_scatter = 1:nelem_ang
                                    for kθ = 1:nop_ang[e_ext]+1
                                        ipθ = connijk_ang[e_ext_scatter,kθ]
                                        ωJac_rad_scatter = ωθ[kθ]*Je_ang[e_ext_scatter,kθ]
                                        intϕ +=   ωJac_rad_scatter*Φ[jpθ,ipθ]             

                                        #S_rad_el[iel,I,J,e_ext,I_ext,J_ext] += ψ_ang[iθ,jθ]*ωJac*ωJac_rad*ωJac_rad_scatter*ψ[j,n]*ψ[i,m]*σ[m,n]*Φ[jpθ,ipθ]/div
                                    end
                                end
                                S_rad_el[iel,I,J,e_ext,I_ext,J_ext] += ψ_ang[iθ,jθ]*intϕ * ψ[i,m] * ψ[j,n] * ωJac*ωJac_rad*σ[m,n]
                                #=for e_ext_scatter = 1:nelem_ang
                                    for kθ = 1:nop_ang[e_ext]+1
                                        ipθ = connijk_ang[e_ext_scatter,kθ]
                                        ωJac_rad_scatter = ωθ[kθ]*Je_ang[e_ext_scatter,kθ]
                                        div = 1
                                        if (kθ == 1 || kθ == 1:nop_ang[e_ext]+1)
                                            div = 2
                                        end
                                                        
                                        S_rad_el[iel,I,J,e_ext,I_ext,J_ext] += ψ_ang[iθ,jθ]*ωJac*ωJac_rad*ωJac_rad_scatter*ψ[j,n]*ψ[i,m]*σ[m,n]*Φ[jpθ,ipθ]/div
                                    end
                                end=#
                            end
                        end
                    end
                end
            end
        end
    end
end

function DSS_angular!(M_el_ang, M, nelem, nop_ang, nelem_ang, connijk_ang, ngl, connijk, iel, ::NSD_2D)
    
    for j=1:ngl
        for i=1:ngl
            I = i +(j-1)*(ngl)
            for e_ext = 1:nelem_ang
                for iθ = 1:nop_ang[e_ext]+1
                    I_ext = iθ 
                    IP_ext = connijk_ang[e_ext,iθ]
                    for n=1:ngl
                        for m=1:ngl
                            J  = m + (n-1)*(ngl)
                            for jθ = 1:nop_ang[e_ext]+1
                                J_ext = jθ 
                                JP_ext = connijk_ang[e_ext,jθ]
                                #@info M[iel, I, J]#, M_el_ang[iel, I, J, e_ext, I_ext, J_ext]
                                M[iel, I, J, IP_ext, JP_ext] += M_el_ang[iel, I, J, e_ext, I_ext, J_ext]
                            end
                        end
                    end
                end
            end
        end
    end
end

function DSS_spatial!(M_sp, M_ang_spat, nelem, npoin_ang, ngl, connijk, inputs, ::NSD_2D)

    for iel = 1:nelem
        for j = 1:ngl
            for i=1:ngl
                I = i +(j-1)*(ngl)
                ip = connijk[iel, i, j]
                if (inputs[:adaptive_extra_meshes])

                else
                    for ip_ext = 1:npoin_ang

                        for n=1:ngl
                            for m = 1:ngl
                                J  = m + (n-1)*(ngl)
                                jp = connijk[iel, m, n]

                                for jp_ext = 1:npoin_ang
                                    M_sp[ip, jp, ip_ext, jp_ext] += M_ang_spat[iel, I, J, ip_ext, jp_ext]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

#=function DSS_spatial(M_ang_spat, npoin, nelem, npoin_ang, ngl, connijk, inputs, ::NSD_2D)
    
    max_entries = npoin*npoin*npoin_ang*npoin_ang
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    K_vec = Vector{Int}()
    L_vec = Vector{Int}()
    V_vec = Vector{Float64}()

        # Reserve space to avoid frequent reallocations
    sizehint!(I_vec, max_entries)
    sizehint!(J_vec, max_entries)
    sizehint!(K_vec, max_entries)
    sizehint!(L_vec, max_entries)
    sizehint!(V_vec, max_entries)

    for iel = 1:nelem
        for j = 1:ngl
            for i=1:ngl
                I = i +(j-1)*(ngl)
                ip = connijk[iel, i, j]
                if (inputs[:adaptive_extra_meshes])

                else
                    for ip_ext = 1:npoin_ang

                        for n=1:ngl
                            for m = 1:ngl
                                J  = m + (n-1)*(ngl)
                                jp = connijk[iel, m, n]

                                for jp_ext = 1:npoin_ang
                                    #M_sp[ip, jp, ip_ext, jp_ext] += M_ang_spat[iel, I, J, ip_ext, jp_ext]
                                    val = M_ang_spat[iel, I, J, ip_ext, jp_ext]
                                    if abs(val) > eps(Float64)  # Skip near-zero entries
                                        push!(I_vec, ip)
                                        push!(J_vec, jp)
                                        push!(K_vec, ip_ext)
                                        push!(L_vec, jp_ext)
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
    return sparse(I_vec, J_vec, K_vec, L_vec, V_vec)
end=#


function Map_rad_global(M_sp, npoin_ang_total, npoin, npoin_ang, nelem, connijk, ngl, inputs, ::NSD_2D)

    if (inputs[:adaptive_extra_meshes])

    else
        max_entries = npoin_ang_total^2
        I_vec = Vector{Int}()
        J_vec = Vector{Int}()
        V_vec = Vector{Float64}()
    
        # Reserve space to avoid frequent reallocations
        sizehint!(I_vec, npoin_ang_total)
        sizehint!(J_vec, npoin_ang_total)
        sizehint!(V_vec, max_entries)

        step_ip = 0 
        for ip=1:npoin
            step_jp = 0
            for jp =1:npoin
                npoin_ext_ip = npoin_ang
                npoin_ext_jp = npoin_ang
                for ip_ext = 1:npoin_ext_ip
                    idx_ip = (ip-1)*(npoin_ext_ip) + ip_ext
                    for jp_ext = 1:npoin_ext_jp
                        idx_jp = (jp-1)*(npoin_ext_jp) + jp_ext
                        #M[idx_ip, idx_jp] = M_sp[ip, jp, ip_ext, jp_ext]
                        val = M_sp[ip, jp, ip_ext, jp_ext]
                        if abs(val) > eps(Float64)  # Skip near-zero entries
                            push!(I_vec, idx_ip)
                            push!(J_vec, idx_jp)
                            push!(V_vec, val)
                        end
                    end
                end
                step_jp +=1
            end
            step_ip +=1
        end
    end
    # Create sparse matrix and sum duplicate entries automatically
    return sparse(I_vec, J_vec, V_vec)
end
#####THIS function is a work in progress until this comment is removed
function sparse_lhs_assembly_3Dby2D(ω, Je, connijk, ωθ, ωϕ, x, y, z,
                                    ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang,
                                   dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, npoin_ang, rad_HG_g)

    max_entries = npoin_ang_total^2
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()

    sizehint!(I_vec, Int64(round(max_entries*0.0009)))
    sizehint!(J_vec, Int64(round(max_entries*0.0009)))
    sizehint!(V_vec, Int64(round(max_entries*0.0009)))
    HG, error = quadgk(v -> (1-rad_HG_g^2)/((1+rad_HG_g^2-2*rad_HG_g*cos(v))^(3/2)), 0, 2*π, rtol=1e-13, atol = 1e-13)    
    for iel=1:nelem
        for k=1:ngl
            for j=1:ngl
                for i=1:ngl
                    ip = connijk[iel,i,j,k]
                    ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                    dξdx_ij = dξdx[iel,i,j,k]
                    dξdy_ij = dξdy[iel,i,j,k]
                    dξdz_ij = dξdz[iel,i,j,k]
                    dηdx_ij = dηdx[iel,i,j,k]
                    dηdy_ij = dηdy[iel,i,j,k]
                    dηdz_ij = dηdz[iel,i,j,k]
                    dζdx_ij = dζdx[iel,i,j,k]
                    dζdy_ij = dζdy[iel,i,j,k]
                    dζdz_ij = dζdz[iel,i,j,k]
                    κ = user_extinction(x[ip],y[ip],z[ip])
                    σ = user_scattering_coef(x[ip],y[ip],z[ip])
                    for e_ext = 1:nelem_ang
                        for jθ = 1:nop_ang[e_ext]+1
                            for iθ = 1:nop_ang[e_ext]+1
                                ip_ext = connijk_ang[e_ext,iθ,jθ]
                                sum = 0.0
                                ωJac_rad = ωθ[iθ]*ωϕ[jθ]*Je_ang[e_ext,iθ,jθ]
                                #@info coords_ang[1,ip_ext], coords_ang[2,ip_ext], e_ext, iθ, jθ
                                for o=1:ngl
                                    for n=1:ngl
                                        for m=1:ngl

                                            jp = connijk[iel,m,n,o]
                                            for lθ = 1:nop_ang[e_ext]+1
                                                for kθ = 1:nop_ang[e_ext]+1
                                                    jp_ext = connijk_ang[e_ext,kθ,lθ]
                                                    extinction = κ*ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]*ψ[k,o]*ψ_ang[iθ,kθ]*ψ_ang[jθ,lθ]
                                                    i_angle = connijk_ang[e_ext,iθ,jθ]

                                                    θ = coords_ang[1,i_angle]
                                                    ϕ = coords_ang[2,i_angle]
                                                     
                                                    propx = ψ[o,k]*ψ[n,j]*dψ[m,i]*dξdx_ij*sin(θ)*cos(ϕ)
                                                    propx += ψ[o,k]*ψ[m,i]*dψ[n,j]*dηdx_ij*sin(θ)*cos(ϕ)
                                                    propx += ψ[n,j]*ψ[m,i]*dψ[o,k]*dζdx_ij*sin(θ)*cos(ϕ) 
                                                    propy = ψ[o,k]*ψ[n,j]*dψ[m,i]*dξdy_ij*sin(θ)*sin(ϕ)
                                                    propy += ψ[o,k]*ψ[m,i]*dψ[n,j]*dηdy_ij*sin(θ)*sin(ϕ)
                                                    propy += ψ[n,j]*ψ[m,i]*dψ[o,k]*dζdy_ij*sin(θ)*sin(ϕ)
                                                    propz = ψ[o,k]*ψ[n,j]*dψ[m,i]*dξdz_ij*cos(θ)
                                                    propz += ψ[o,k]*ψ[m,i]*dψ[n,j]*dηdz_ij*cos(θ)
                                                    propz += ψ[n,j]*ψ[m,i]*dψ[o,k]*dζdz_ij*cos(θ)
                                                    propagation = (ωJac*ωJac_rad)*ψ_ang[iθ,kθ]*ψ_ang[jθ,lθ]*(propx+propy+propz)
                                                    intΦ = 0.0 
                                                    sum += propagation#propy*(ωJac*ωJac_rad)*ψ_ang[iθ,kθ]*ψ_ang[jθ,lθ]
                                                    for e_ext_scatter = 1:nelem_ang
                                                        for nθ = 1:nop_ang[e_ext]+1
                                                            for mθ = 1:nop_ang[e_ext]+1
                                                                div = 1
                                                                #=if (mθ == nop_ang[e_ext]+1 || mθ == 1) && (nθ == nop_ang[e_ext]+1 || nθ == 1)
                                                                    div = 3
                                                                #elseif (mθ == nop_ang[e_ext]+1 || mθ == 1) || (nθ == nop_ang[e_ext]+1 || nθ == 1)
                                                                 #   div = 2
                                                                end=#
                                                                ipθ = connijk_ang[e_ext_scatter,mθ,nθ]
                                                                θ1 = coords_ang[1,ipθ]
                                                                ϕ1 = coords_ang[2,ipθ]

                                                                #=if (div > 1) && (θ1 < eps(Float64) || abs(θ1 - π) < eps(Float64))
                                                                    div = div/2
                                                                end
                                                                if (div > 1) && (ϕ1 < eps(Float64) && nθ == nop_ang[e_ext]+1)
                                                                    div = div/2
                                                                end=#
                                                                Φ = user_scattering_functions(θ,θ1,ϕ,ϕ1,HG)
                                                                ωJac_rad_scatter = ωθ[mθ]*ωϕ[nθ]*Je_ang[e_ext_scatter,mθ,nθ]
                                                                intΦ +=   ωJac_rad_scatter*Φ/div
                                                                #@info intΦ, Je_ang[e_ext_scatter,mθ,nθ], θ1, ϕ1, e_ext_scatter, mθ, nθ
                                                            end
                                                        end
                                                    end
                                                    @info intΦ
                                                    scattering = ψ_ang[iθ,kθ] * ψ_ang[jθ,lθ] * intΦ * ψ[k,o] * ψ[i,m] * ψ[j,n] * ωJac*ωJac_rad*σ
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
                end
            end
        end
    end
    return sparse(I_vec, J_vec, V_vec)
end

function sparse_mass_assembly_3Dby2D(ω, Je, connijk, ωθ, ωϕ, x, y, ψ, dψ, ψ_ang, 
        connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang, npoin_ang)
  
    max_entries = npoin_ang_total^2
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()

    sizehint!(I_vec, Int64(round(max_entries*0.0001)))
    sizehint!(J_vec, Int64(round(max_entries*0.0001)))
    sizehint!(V_vec, Int64(round(max_entries*0.0001)))

    for iel=1:nelem
        for k=1:ngl
            for j=1:ngl
                for i=1:ngl
                    ip = connijk[iel,i,j,k]
                    ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                    for e_ext = 1:nelem_ang
                        for jθ = 1:nop_ang[e_ext]+1
                            for iθ = 1:nop_ang[e_ext]+1
                                ωJac_rad = ωθ[iθ]*ωϕ[jθ]*Je_ang[e_ext,iθ,jθ]
                                ip_ext = connijk_ang[e_ext,iθ,jθ]
                                for o=1:ngl
                                    for n=1:ngl
                                        for m=1:ngl
                                            jp = connijk[iel,m,n,o]
                                            for lθ = 1:nop_ang[e_ext]+1
                                                for kθ = 1:nop_ang[e_ext]+1
                                                    jp_ext = connijk_ang[e_ext,kθ,lθ]

                                                    val = ωJac*ωJac_rad*ψ[k,o]*ψ[j,n]*ψ[i,m]*ψ_ang[iθ,kθ]*ψ_ang[jθ,lθ]
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
                end
            end
        end
    end
    return sparse(I_vec, J_vec, V_vec)
end
