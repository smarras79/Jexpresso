function compute_viscosity!(μ::Vector{Float64}, ::NSD_1D, PT::ShallowWater, q, q1, q2, rhs, Δt, mesh,metrics)

    #compute domain averages
    H_avg = 0.0
    Hu_avg = 0.0
    for ie=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[i,ie]
            H_avg += q[ip,1]
            Hu_avg += q[ip,2]
        end
    end
    H_avg = H_avg / (mesh.nelem*mesh.ngl)
    Hu_avg = Hu_avg / (mesh.nelem*mesh.ngl)

    #Get denominator infinity norms

    Hdiff = zeros(mesh.ngl,mesh.nelem)
    Hudiff = zeros(mesh.ngl,mesh.nelem)
    for ie=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[i,ie]
            Hdiff[i,ie] = abs(q[ip,1] - H_avg)
            Hudiff[i,ie] = abs(q[ip,2] - Hu_avg)
        end
    end
    denom1 = maximum(Hdiff)  + 1.0e-16
    denom2 = maximum(Hudiff) + 1.0e-16
    #@info denom1, denom2
    #Get Numerator inifinity norms, mu_max infinity norm
    #mu_SGS = zeros(mesh.nelem,1)

    for ie =1:mesh.nelem
        Δ = mesh.Δx[ie]/mesh.ngl
        RH = zeros(mesh.ngl)
        RHu = zeros(mesh.ngl)
        ughs = zeros(mesh.ngl)
        for i=1:mesh.ngl
            ip = mesh.conn[i,ie]
            x = mesh.x[ip]
            #@info "mass", (3*q[ip,1]-4*q1[ip,1]+q2[ip,1])/(2*Δt),q[ip,1],q1[ip,1],q2[ip,1]
            #@info "velocity", (3*q[ip,2]-4*q1[ip,2]+q2[ip,2])/(2*Δt),q[ip,2],q1[ip,2],q2[ip,2]      
            RH[i] = abs((q[ip,1] - q1[ip,1])/Δt + rhs[ip,1]) #abs((3*q[ip,1]-4*q1[ip,1]+q2[ip,1])/(2*Δt)+rhs[ip,1])#rhs[ip,1] #abs((q[ip,1] - q1[ip,1])/Δt + rhs[ip,1])
            RHu[i] = abs((q[ip,2] - q1[ip,2])/Δt + rhs[ip,2])#abs((3*q[ip,2]-4*q1[ip,2]+q2[ip,2])/(2*Δt)+rhs[ip,2])#rhs[ip,2] #(q[ip,2] - q1[ip,2])/Δt + rhs[ip,2]
            ughs[i] = abs(q[ip,2]/q[ip,1])+sqrt(9.81*(max(q[ip,1],0.001)))
        end
        numer1 = maximum(RH)
        numer2 = maximum(RHu)
        #@info numer1, numer2
        mu_res = Δ^2 * max(numer1/denom1, numer2/denom2)
        mu_max = 0.5*Δ*maximum(ughs)
        μ[ie] = max(0.0,min(mu_max,mu_res))
    end 

    #return mu_SGS

end

function compute_viscosity!(μ::Vector{Float64}, ::NSD_1D, PT::AdvDiff, q, q1, q2, rhs, Δt, mesh,metrics)

    #compute domain averages
    ρ_avg  = 0.0
    #ρu_avg = 0.0
    #ρE_avg = 0.0
    for e=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[i,e]
            ρ_avg  += q[ip,1]
            #ρu_avg += q[ip,2]
            #ρu_avg += q[ip,3]
        end
    end
    ρ_avg = ρ_avg / (mesh.nelem*mesh.ngl)

    #Get denominator infinity norms
    ρdiff  = zeros(mesh.ngl,mesh.nelem)
    for ie=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[i,ie]
            ρdiff[i,ie]  = abs(q[ip,1] - ρ_avg)
        end
    end
    denom1 = maximum(ρdiff)  + 1.0e-16
    #denom2 = maximum(ρudiff) + 1.0e-16
   
    #Get Numerator inifinity norms, mu_max infinity norm
    #mu_SGS = zeros(mesh.nelem,1)
                           
    for ie =1:mesh.nelem
        Δ = mesh.Δx[ie]/mesh.ngl
        Rρ = zeros(mesh.ngl)
        #Rρu = zeros(mesh.ngl)
        for i=1:mesh.ngl
            ip = mesh.conn[i,ie]
            x = mesh.x[ip]
            Rρ[i] = abs((q[ip,1] - q1[ip,1])/Δt - rhs[ip,1])
            Rρu[i] = (q[ip,2] - q1[ip,2])/Δt - rhs[ip,2]
        end
        numer1 = maximum(Rρ)
        #numer2 = maximum(Rρ)
        #@info numer1, numer2
        mu_res = Δ^2*numer1/denom1
        #mu_max = 0.5*Δ*maximum(ughs)
        μ[ie] = max(0.0, mu_res)
    end
                          
end


function compute_viscosity!(μ_SGS, ::NSD_1D, PT::CompEuler, q, q1, q2, rhs, Δt, mesh, metrics, T)

    #compute domain averages
    ρ_avg  = 0.0
    ρu_avg = 0.0
    ρE_avg = 0.0
    for e=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[i,e]
            ρ_avg  += q[ip,1]
            ρu_avg += q[ip,2]
            ρE_avg += q[ip,3]
        end
    end
    ρ_avg  = ρ_avg  / (mesh.nelem*mesh.ngl)
    ρu_avg = ρu_avg / (mesh.nelem*mesh.ngl)
    ρE_avg = ρE_avg / (mesh.nelem*mesh.ngl)

    #Get denominator infinity norms
    ρdiff  = zeros(T, mesh.ngl,mesh.nelem)
    ρdiff  = zeros(T, mesh.ngl,mesh.nelem)
    ρudiff = zeros(T, mesh.ngl,mesh.nelem)
    ρEdiff = zeros(T, mesh.ngl,mesh.nelem)
    for e=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[i,e]
            
            ρdiff[i,e]  = abs(q[ip,1] - ρ_avg)
            ρudiff[i,e] = abs(q[ip,2] - ρu_avg)
            ρEdiff[i,e] = abs(q[ip,3] - ρE_avg)
        end
    end
    denom1 = maximum(ρdiff)  + 1.0e-16
    denom2 = maximum(ρudiff) + 1.0e-16
    denom3 = maximum(ρEdiff) + 1.0e-16
        
    ρ   = @MVector zeros(Float64, mesh.ngl)
    u   = @MVector zeros(Float64, mesh.ngl)
    T   = @MVector zeros(Float64, mesh.ngl)
    
    Rρ  = @MVector zeros(Float64, mesh.ngl)
    Rρu = @MVector zeros(Float64, mesh.ngl)
    RρE = @MVector zeros(Float64, mesh.ngl)
    for ie =1:mesh.nelem
        Δ = mesh.Δx[ie]/mesh.ngl
        for i=1:mesh.ngl
            ip = mesh.conn[i,ie]
            
             Rρ[i] = abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])/(2*Δt) + rhs[ip,1])
            Rρu[i] = abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])/(2*Δt) + rhs[ip,2])
            RρE[i] = abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])/(2*Δt) + rhs[ip,3])
            
            ρ[i] = q[ip,1]
            u[i] = q[ip,2]/ρ[i]
            e    = q[ip,3]/ρ[i]
            T[i] = e - 0.5*u[i]^2
        end
        γ = 1.4
        C1 = 1.0
        C2 = 0.5
        numer1 = maximum(Rρ)
        numer2 = maximum(Rρu)
        numer3 = maximum(RρE)
        #@info numer1, numer2
        μ_res = C1*Δ^2*denom1*max(numer1/denom1, numer2/denom2, numer3/denom3)
        μ_max = C2*Δ*maximum(ρ)*maximum(abs.(u) .+ sqrt.(γ*T))
        μ_SGS[ie] = max(0.0, min(μ_max, μ_res))
    end
    
end

function compute_viscosity!(μ::Vector{Float64}, ::NSD_2D, PT::CompEuler, q, q1, q2, rhs, Δt, mesh, metrics)
    
    PhysConst = PhysicalConst{Float64}()
    #compute domain averages
    ρ_avg  = 0.0
    ρu_avg = 0.0
    ρv_avg = 0.0
    ρE_avg = 0.0
    for e=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,e]
                ρ_avg  += q[ip,1]
                ρu_avg += q[ip,2]
                ρv_avg += q[ip,3]
                ρE_avg += q[ip,4]
            end
        end
    end
    ρ_avg  = ρ_avg  / (mesh.nelem*mesh.ngl*mesh.ngl)
    ρu_avg = ρu_avg / (mesh.nelem*mesh.ngl*mesh.ngl)
    ρv_avg = ρv_avg / (mesh.nelem*mesh.ngl*mesh.ngl)
    ρE_avg = ρE_avg / (mesh.nelem*mesh.ngl*mesh.ngl)

    #Get denominator infinity norms
    ρdiff  = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    ρudiff = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    ρvdiff = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    ρEdiff = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    for e=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,e]
            
                ρdiff[i,j,e]  = abs(q[ip,1] - ρ_avg)
                ρudiff[i,j,e] = abs(q[ip,2] - ρu_avg)
                ρvdiff[i,j,e] = abs(q[ip,3] - ρv_avg)
                ρEdiff[i,j,e] = abs(q[ip,4] - ρE_avg)
            end
        end
    end
    denom1 = maximum(ρdiff)  + 1.0e-16
    denom2 = maximum(ρudiff) + 1.0e-16
    denom3 = maximum(ρvdiff) + 1.0e-16
    denom4 = maximum(ρEdiff) + 1.0e-16
    #@info denom1 denom2 denom3
    
    #Get Numerator inifinity norms, μ_max infinity norm
    μ_SGS = zeros(mesh.nelem,1)
    
    for ie =1:mesh.nelem
        Δ   = (maximum(mesh.x) - minimum(mesh.x))/(25*mesh.nop) #temporary. make general ASAP
        ρ_bar   = 0.0#zeros(mesh.ngl,mesh.ngl)
        u   = zeros(mesh.ngl,mesh.ngl)
        v   = zeros(mesh.ngl,mesh.ngl)
        T   = zeros(mesh.ngl,mesh.ngl)
        e   = zeros(mesh.ngl,mesh.ngl)
        p_bar = 0.0
        
        Rρ  = zeros(mesh.ngl,mesh.ngl)
        Rρu = zeros(mesh.ngl,mesh.ngl)
        Rρv = zeros(mesh.ngl,mesh.ngl)
        RρE = zeros(mesh.ngl,mesh.ngl)
        for j=1:mesh.ngl
            for i=1:mesh.ngl
                ip = mesh.connijk[i,j,ie]
                
                #Rρ[i] = abs((q[ip,1] - q1[ip,1])/Δt + rhs[ip,1]) #abs((3*q[ip,1]-4*q1[ip,1]+q2[ip,1])/(2*Δt)+rhs[ip,1])#rhs[ip,1] #abs((q[ip,1] - q1[ip,1])/Δt + rhs[ip,1])
                #Rρu[i] = abs((q[ip,2] - q1[ip,2])/Δt + rhs[ip,2])#abs((3*q[ip,2]-4*q1[ip,2]+q2[ip,2])/(2*Δt)+rhs[ip,2])#rhs[ip,2] #(q[ip,2] - q1[ip,2])/Δt + rhs[ip,2]
                #RρE[i] = abs((q[ip,3] - q1[ip,3])/Δt + rhs[ip,3])#abs((3*q[ip,3]-4*q1[ip,3]+q2[ip,3])/(2*Δt)+rhs[ip,2])#rhs[ip,3] #(q[ip,2] - q1[ip,2])/Δt + rhs[ip,2]

                Rρ[i,j] = abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])/(2*Δt) + rhs[ip,1])
                Rρu[i,j] = abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])/(2*Δt) + rhs[ip,2])
                Rρv[i,j] = abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])/(2*Δt) + rhs[ip,3])
                RρE[i,j] = abs((3*q[ip,4] - 4*q1[ip,4] + q2[ip,4])/(2*Δt) + rhs[ip,4])
                
                ρ_bar += q[ip,1]
                u[i,j] = q[ip,2]/q[ip,1]
                v[i,j] = q[ip,3]/q[ip,1]
                e[i,j] = q[ip,4]/q[ip,1]
                #T[i,j] = e[i] - 0.5*(u[i,j]^2 + v[i,j]^2)
                p_bar += perfectGasLaw_ρθtoP(PhysConst, ρ= q[ip,1], θ=e[i,j]) 
            end
        end
        ρ_bar = ρ_bar/mesh.ngl^2
        p_bar = p_bar/mesh.ngl^2 
        γ = 1.4
        C1 = 1.0
        C2 = 0.5
        numer1 = maximum(Rρ)
        numer2 = maximum(Rρu)
        numer3 = maximum(Rρv)
        numer4 = maximum(RρE)
        μ_res = C1*Δ^2*denom1*max(numer1/denom1, numer2/denom2, numer3/denom3, numer4/denom4)
        μ_max = C2*Δ*maximum(sqrt.(u.*u .+ v.*v) .+ sqrt(γ*p_bar/ρ_bar))
        μ_SGS[ie] = max(0.0, min(μ_max, μ_res))
    end

    return μ_SGS

end
