using Random



function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

            """
    @info " Initialize fields for 3D CompEuler with θ equation ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ("ρ", "ρu", "ρv", "ρw", "ρθ")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    
    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()
        if inputs[:lrestart] == true
            #
            # READ RESTART HDF5:
            #
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars))
            PhysConst = PhysicalConst{Float64}()
        
            for ip=1:mesh.npoin
                ρ  = q.qn[ip,1]
                ρθ = q.qn[ip,5]
                θ  = ρθ/ρ
                P = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
                q.qn[ip,end] = P
            
                ρe  = q.qe[ip,1]
                ρθe = q.qe[ip,5]
                θe  = ρθe/ρ
                Pe = perfectGasLaw_ρθtoP(PhysConst, ρ=ρe, θ=θe)
                q.qe[ip,end] = Pe
            end
        
        else
            #
            # INITIAL STATE from scratch:
            #
            new_param_set = create_updated_TD_Parameters(FT(101325.0))
            for ip = 1:mesh.npoin
            
                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]
            
            
                # @info param_set
                Δθ = 0.0 #K
    
                θ = 0.0
                θref = 0.0
                p    = 0.0
                pref = 0.0
                ρ    = 0.0
                ρref = 0.0
                
                u = 0.0
                v = 0.0
                w = 0.0
                ρref, u, v, w, θref, Δθ, pref = initialize_bomex!(z, new_param_set)
                # @info ρref, u, v, w, θref, Δθ, pref
                # ρref = params.ρ
                # u = params.u
                # v = params.v
                # w = params.w
                # θref = params.θ
                # Δθ = params.Δθ #K
                # pref = params.P
                ρ = ρref
                θ = θref + Δθ
                p = pref
            
                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*w - ρref*w
                    q.qn[ip,5] = ρ*θ - ρref*θref
                    q.qn[ip,end] = p
                
                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = u
                    q.qe[ip,3] = v
                    q.qe[ip,4] = w
                    q.qe[ip,5] = ρref*θref
                    q.qe[ip,end] = pref
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*w
                    q.qn[ip,5] = ρ*θ
                    q.qn[ip,end] = p

                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = ρref*u
                    q.qe[ip,3] = ρref*v
                    q.qe[ip,4] = ρref*w
                    q.qe[ip,5] = ρref*θref
                    q.qe[ip,end] = pref
                end
                #end
            end
        end
    
        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[:,2] .= q.qn[:,2]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,3] .= q.qn[:,3]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,4] .= q.qn[:,4]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,5] .= q.qn[:,5]./(q.qn[:,1] + q.qe[:,1])
            
                #Store initial background state for plotting and analysis of pertuebations
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]
                q.qn[:,5] .= q.qn[:,5]./q.qn[:,1]

                #Store initial background state for plotting and analysis of pertuebations
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            end
        end

        if (inputs[:lwrite_initial] == true)
            outvarsref = ("rho_ref", "uρ_ref", "vρ_ref", "wρ_ref", "theta_ref", "p_ref")    
            write_vtk_ref(SD, mesh, q.qn.-q.qe, "initial_state", inputs[:output_dir]; nvar=length(q.qn[1,:]), outvarsref=outvarsref)
        
            outvarsref = ("rho_ref", "uρ_ref", "vρ_ref", "wρ_ref", "theta_ref", "p_ref")    
            write_vtk_ref(SD, mesh, q.qe, "REFERENCE_state", inputs[:output_dir]; nvar=length(q.qe[1,:]), outvarsref=outvarsref)
        end
    else

        if (inputs[:SOL_VARS_TYPE] == PERT())
            lpert = true
        else
            lpert = false
        end
        PhysConst = PhysicalConst{TFloat}()
        FT = TFloat
        new_param_set = create_updated_TD_Parameters(FT(101325.0))

        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, mesh.z, PhysConst, new_param_set, lpert; ndrange = (mesh.npoin))
    end
    @info " Initialize fields for 3D CompEuler with θ equation ........................ DONE "
    
    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, z, PhysConst, param_set, lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    x = x[ip]
    y = y[ip]
    z = z[ip]


    Δθ = T(0.0) #K
    
    θ = T(0.0)
    θref = T(0.0)
    p    = T(0.0)
    pref = T(0.0)
    ρ    = T(0.0)
    ρref = T(0.0)
    
    u = T(0.0)
    v = T(0.0)
    w = T(0.0)
    ρref, u, v, w, θref, Δθ, pref = initialize_bomex!(z, param_set)

    # u = 0
    ρ = ρref
    θ = θref + Δθ
    p = pref

    if (lpert)
        qn[ip,1] = ρ - ρref
        qn[ip,2] = ρ*u - ρref*u
        qn[ip,3] = ρ*v - ρref*v
        qn[ip,4] = ρ*w - ρref*w
        qn[ip,5] = ρ*θ - ρref*θref
        qn[ip,end] = p
    else
        qn[ip,1] = ρ
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*w
        qn[ip,5] = ρ*θ
        qn[ip,end] = p
    end

                    #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = ρref
    qe[ip,2] = ρref*u
    qe[ip,3] = ρref*v
    qe[ip,4] = ρref*w
    qe[ip,5] = ρref*θref
    qe[ip,end] = pref

end

function initialize_bomex!(z, param_set)
    FT = eltype(z)
    @inbounds begin
        P_sfc::FT = 1.015e5 # Surface air pressure
        qg::FT = 22.45e-3 # Total moisture at surface
        q_pt_sfc = TD.PhasePartition(qg) # Surface moisture partitioning
        Rm_sfc = TD.gas_constant_air(param_set, q_pt_sfc) # Moist gas constant
        T_sfc = FT(300.4) # Surface temperature
        _grav = FT(TP.grav(param_set))
        H = Rm_sfc * T_sfc / _grav
        # z = FT(500.0)
        P = P_sfc * exp(-z / H)
        
        zlv::FT = 700
        zl1::FT = 520
        zl2::FT = 1480
        zl3::FT = 2000
        zl4::FT = 3000
        u = FT(0.0)
        # if z <= zlv
        #     u = -8.75
        # else
        #     u = -8.75 + (z - zlv) * (-4.61 + 8.75) / (zl4 - zlv)
        # end
        v = FT(0.0)
        w = FT(0.0)

        # Assign piecewise quantities to θ_liq and q_tot
        θ_liq::FT = 0
        q_tot::FT = 0

        # Piecewise functions for potential temperature and total moisture
        if FT(0) <= z <= zl1
            # Well mixed layer
            θ_liq = 298.7
            q_tot = 17.0 + (z / zl1) * (16.3 - 17.0)
        elseif z > zl1 && z <= zl2
            # Conditionally unstable layer
            θ_liq = 298.7 + (z - zl1) * (302.4 - 298.7) / (zl2 - zl1)
            q_tot = 16.3 + (z - zl1) * (10.7 - 16.3) / (zl2 - zl1)
        elseif z > zl2 && z <= zl3
            # Absolutely stable inversion
            θ_liq = 302.4 + (z - zl2) * (308.2 - 302.4) / (zl3 - zl2)
            q_tot = 10.7 + (z - zl2) * (4.2 - 10.7) / (zl3 - zl2)
        else
            θ_liq = 308.2 + (z - zl3) * (311.85 - 308.2) / (zl4 - zl3)
            q_tot = 4.2 + (z - zl3) * (3.0 - 4.2) / (zl4 - zl3)
        end

        # Convert total specific humidity to kg/kg
        q_tot /= 1000


        # Establish thermodynamic state and moist phase partitioning
        TS = TD.PhaseEquil_pθq(param_set, P, θ_liq, q_tot)
        T = TD.air_temperature(param_set, TS)
        ρ = TD.air_density(param_set, TS)
        q_pt = TD.PhasePartition(param_set, TS)
        θ = TD.virtual_pottemp(param_set, TS)
        bhasCondense = TD.has_condensate(param_set, TS)
        # @info TS
        Δθ = FT(0.0)
        if z < FT(400)
            Δθ = FT((rand()-0.5) * 1.0)
        end
    end
    return ρ, u, v, w, θ, Δθ, P
    # if bhasCondense
    #     @info z, ρ, u, T, θ_liq, θ, q_tot, q_pt.liq, P, bhasCondense
    # end
end

