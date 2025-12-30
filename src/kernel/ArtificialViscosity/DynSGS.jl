function compute_viscosity!(μ::Vector{Float64}, ::NSD_1D, PT::ShallowWater, q, q1, q2, rhs, Δt, mesh,metrics)

    #compute domain averages
    H_avg = 0.0
    Hu_avg = 0.0
    @inbounds for ie=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[ie,i]
            H_avg += q[ip,1]
            Hu_avg += q[ip,2]
        end
    end
    H_avg = H_avg / (mesh.nelem*mesh.ngl)
    Hu_avg = Hu_avg / (mesh.nelem*mesh.ngl)

    #Get denominator infinity norms
    Hdiff = zeros(mesh.ngl,mesh.nelem)
    Hudiff = zeros(mesh.ngl,mesh.nelem)
    @inbounds for ie=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[ie,i]
            Hdiff[i,ie] = abs(q[ip,1] - H_avg)
            Hudiff[i,ie] = abs(q[ip,2] - Hu_avg)
        end
    end
    denom1 = maximum(Hdiff)  + 1.0e-16
    denom2 = maximum(Hudiff) + 1.0e-16

    #Get Numerator inifinity norms, mu_max infinity norm
    # Pre-allocate temporary arrays outside loop for performance - CRITICAL FIX
    RH = zeros(mesh.ngl)
    RHu = zeros(mesh.ngl)
    ughs = zeros(mesh.ngl)

    @inbounds for ie =1:mesh.nelem
        Δ = mesh.Δx[ie]/mesh.ngl
        fill!(RH, 0.0)
        fill!(RHu, 0.0)
        fill!(ughs, 0.0)

        @simd for i=1:mesh.ngl
            ip = mesh.conn[ie,i]
            RH[i] = abs((q[ip,1] - q1[ip,1])/Δt + rhs[ip,1])
            RHu[i] = abs((q[ip,2] - q1[ip,2])/Δt + rhs[ip,2])
            ughs[i] = abs(q[ip,2]/q[ip,1])+sqrt(9.81*(max(q[ip,1],0.001)))
        end
        numer1 = maximum(RH)
        numer2 = maximum(RHu)
        mu_res = Δ^2 * max(numer1/denom1, numer2/denom2)
        mu_max = 0.5*Δ*maximum(ughs)
        μ[ie] = max(0.0,min(mu_max,mu_res))
    end

end

function compute_viscosity!(μ::Vector{Float64}, ::NSD_1D, PT::AdvDiff, q, q1, q2, rhs, Δt, mesh,metrics)

    #compute domain averages
    ρ_avg  = 0.0
    for e=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[i,e]
            ρ_avg  += q[ip,1]
        end
    end
    ρ_avg = ρ_avg / (mesh.nelem*mesh.ngl)

    #Get denominator infinity norms
    ρdiff  = zeros(mesh.ngl,mesh.nelem)
    for ie=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[ie,i]
            ρdiff[i,ie]  = abs(q[ip,1] - ρ_avg)
        end
    end
    denom1 = maximum(ρdiff)  + 1.0e-16

    #Get Numerator inifinity norms, mu_max infinity norm
    # Pre-allocate temporary arrays outside loop for performance
    Rρ = zeros(mesh.ngl)

    @inbounds for ie =1:mesh.nelem
        Δ = mesh.Δx[ie]/mesh.ngl
        fill!(Rρ, 0.0)  # Reuse pre-allocated array
        @simd for i=1:mesh.ngl
            ip = mesh.conn[ie,i]
            Rρ[i] = abs((q[ip,1] - q1[ip,1])/Δt - rhs[ip,1])
        end
        numer1 = maximum(Rρ)
        mu_res = Δ^2*numer1/denom1
        μ[ie] = max(0.0, mu_res)
    end

end

function compute_viscosity!(μ_SGS, ::NSD_1D, PT::CompEuler, q, q1, q2, rhs, Δt, mesh, metrics, TT)

    #compute domain averages
    ρ_avg  = zero(TT)
    ρu_avg = zero(TT)
    ρE_avg = zero(TT)
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
    # Fixed: removed duplicate allocation
    ρdiff  = zeros(TT, mesh.ngl,mesh.nelem)
    ρudiff = zeros(TT, mesh.ngl,mesh.nelem)
    ρEdiff = zeros(TT, mesh.ngl,mesh.nelem)
    for e=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.conn[i,e]

            ρdiff[i,e]  = abs(q[ip,1] - ρ_avg)
            ρudiff[i,e] = abs(q[ip,2] - ρu_avg)
            ρEdiff[i,e] = abs(q[ip,3] - ρE_avg)
        end
    end
    denom1 = maximum(ρdiff)  + TT(1.0e-16)
    denom2 = maximum(ρudiff) + TT(1.0e-16)
    denom3 = maximum(ρEdiff) + TT(1.0e-16)

    # Fixed: use type parameter TT instead of hardcoded Float64
    # Renamed T to Temp to avoid shadowing type parameter
    ρ    = @MVector zeros(TT, mesh.ngl)
    u    = @MVector zeros(TT, mesh.ngl)
    Temp = @MVector zeros(TT, mesh.ngl)

    Rρ  = @MVector zeros(TT, mesh.ngl)
    Rρu = @MVector zeros(TT, mesh.ngl)
    RρE = @MVector zeros(TT, mesh.ngl)

    γ = TT(1.4)
    C1 = TT(1.0)
    C2 = TT(0.5)

    @inbounds for ie =1:mesh.nelem
        Δ = mesh.Δx[ie]/mesh.ngl
        @simd for i=1:mesh.ngl
            ip = mesh.conn[ie,i]

            Rρ[i] = abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])/(2*Δt) + rhs[ip,1])
            Rρu[i] = abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])/(2*Δt) + rhs[ip,2])
            RρE[i] = abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])/(2*Δt) + rhs[ip,3])

            ρ[i] = q[ip,1]
            u[i] = q[ip,2]/ρ[i]
            e    = q[ip,3]/ρ[i]
            Temp[i] = e - TT(0.5)*u[i]^2
        end
        numer1 = maximum(Rρ)
        numer2 = maximum(Rρu)
        numer3 = maximum(RρE)
        μ_res = C1*Δ^2*denom1*max(numer1/denom1, numer2/denom2, numer3/denom3)
        μ_max = C2*Δ*maximum(ρ)*maximum(abs.(u) .+ sqrt.(γ*Temp))
        μ_SGS[ie] = max(zero(TT), min(μ_max, μ_res))
    end

end

function compute_viscosity!(μ::Vector{Float64}, ::NSD_2D, PT::CompEuler, q, q1, q2, rhs, Δt, mesh, metrics)

    PhysConst = PhysicalConst{Float64}()
    #compute domain averages
    ρ_avg  = 0.0
    ρu_avg = 0.0
    ρv_avg = 0.0
    ρE_avg = 0.0
    @inbounds for e=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[e,i,j]
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

    #Get denominator infinity norms - Fixed: correct index order for column-major
    ρdiff  = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    ρudiff = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    ρvdiff = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    ρEdiff = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    @inbounds for e=1:mesh.nelem
        for j=1:mesh.ngl
            for i=1:mesh.ngl
                ip = mesh.connijk[e,i,j]

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

    #Get Numerator inifinity norms, μ_max infinity norm
    μ_SGS = zeros(mesh.nelem,1)

    # Pre-allocate arrays outside loop for performance - CRITICAL FIX
    u   = zeros(mesh.ngl,mesh.ngl)
    v   = zeros(mesh.ngl,mesh.ngl)
    T   = zeros(mesh.ngl,mesh.ngl)
    e   = zeros(mesh.ngl,mesh.ngl)

    Rρ  = zeros(mesh.ngl,mesh.ngl)
    Rρu = zeros(mesh.ngl,mesh.ngl)
    Rρv = zeros(mesh.ngl,mesh.ngl)
    RρE = zeros(mesh.ngl,mesh.ngl)

    γ = 1.4
    C1 = 1.0
    C2 = 0.5

    @inbounds for ie =1:mesh.nelem
        # Fixed: Use element-local grid spacing if available, otherwise compute properly
        # TODO: Replace this with proper element size from metrics
        if hasfield(typeof(mesh), :Δeffective_2d) && length(mesh.Δeffective_2d) > 0
            Δ = mesh.Δeffective_2d[ie]
        else
            # Fallback: compute from element extent (not ideal but better than global)
            Δ = sqrt(metrics.Je[ie,1,1])  # Approximate from Jacobian
        end

        ρ_bar   = 0.0
        p_bar   = 0.0

        # Reuse pre-allocated arrays
        fill!(u, 0.0)
        fill!(v, 0.0)
        fill!(e, 0.0)
        fill!(Rρ, 0.0)
        fill!(Rρu, 0.0)
        fill!(Rρv, 0.0)
        fill!(RρE, 0.0)

        for j=1:mesh.ngl
            @simd for i=1:mesh.ngl
                ip = mesh.connijk[i,j,ie]

                Rρ[i,j] = abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])/(2*Δt) + rhs[ip,1])
                Rρu[i,j] = abs((3*q[ip,2] - 4*q1[ip,2] + q2[ip,2])/(2*Δt) + rhs[ip,2])
                Rρv[i,j] = abs((3*q[ip,3] - 4*q1[ip,3] + q2[ip,3])/(2*Δt) + rhs[ip,3])
                RρE[i,j] = abs((3*q[ip,4] - 4*q1[ip,4] + q2[ip,4])/(2*Δt) + rhs[ip,4])

                ρ_bar += q[ip,1]
                u[i,j] = q[ip,2]/q[ip,1]
                v[i,j] = q[ip,3]/q[ip,1]
                e[i,j] = q[ip,4]/q[ip,1]
                p_bar += perfectGasLaw_ρθtoP(PhysConst, ρ= q[ip,1], θ=e[i,j])
            end
        end
        ρ_bar = ρ_bar/mesh.ngl^2
        p_bar = p_bar/mesh.ngl^2

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
