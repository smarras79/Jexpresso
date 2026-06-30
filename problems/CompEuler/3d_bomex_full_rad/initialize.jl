using Random

function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 3D BOMEX with full RT coupling ........................ "
    end

    #---------------------------------------------------------------------------------
    # Variables: static energy formulation + two moisture variables
    #   ρhl = ρ × (cp T + g z)  -- liquid-water static energy (no latent heat)
    #   ρqt = ρ × total water mixing ratio
    #   ρqp = ρ × precipitating water mixing ratio
    #---------------------------------------------------------------------------------
    qvars    = ("ρ", "ρu", "ρv", "ρw", "hl", "ρqt", "ρqp")
    qoutvars = ["ρ", "ρu", "ρv", "ρw", "hl", "ρqt", "ρqp",
                "T", "qn", "qc", "qi", "qr", "qs", "qg", "qsatt"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat,
                 inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    if inputs[:backend] == CPU()
        PhysConst = PhysicalConst{Float64}()
        new_param_set = create_updated_TD_Parameters(
            PhysConst.potential_temperature_reference_pressure)

        if inputs[:lrestart] == true
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path],
                                     inputs, mesh.npoin, HDF5(); nvar=length(qvars))
        else
            for ip = 1:mesh.npoin
                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]

                ρref, u, v, w, e_tot_ref, Δe, P, qt_ref, Δqt, ql_ref =
                    initialize_bomex!(z, new_param_set)

                # Static energy: hl = cp T + g z
                # Recover T from total energy: e_tot = cp T + g z + ½(u²+v²+w²)
                # For the background state (u=v=w=0): e_tot_ref = cp T_ref + g z
                hl_ref = e_tot_ref  # background: KE=0, so e_tot = hl for ref state

                # Perturbed state
                ρ  = ρref
                hl = hl_ref + Δe
                qt = qt_ref + Δqt
                ql = ql_ref

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip, 1] = ρ - ρref
                    q.qn[ip, 2] = ρ*u - ρref*u
                    q.qn[ip, 3] = ρ*v - ρref*v
                    q.qn[ip, 4] = ρ*w - ρref*w
                    q.qn[ip, 5] = ρ*hl - ρref*hl_ref
                    q.qn[ip, 6] = ρ*qt - ρref*qt_ref
                    q.qn[ip, 7] = 0.0
                    q.qn[ip, end] = P

                    q.qe[ip, 1] = ρref
                    q.qe[ip, 2] = ρref*u
                    q.qe[ip, 3] = ρref*v
                    q.qe[ip, 4] = ρref*w
                    q.qe[ip, 5] = ρref*hl_ref
                    q.qe[ip, 6] = ρref*qt_ref
                    q.qe[ip, 7] = 0.0
                    q.qe[ip, end] = P
                else
                    q.qn[ip, 1] = ρ
                    q.qn[ip, 2] = ρ*u
                    q.qn[ip, 3] = ρ*v
                    q.qn[ip, 4] = ρ*w
                    q.qn[ip, 5] = ρ*hl
                    q.qn[ip, 6] = ρ*qt
                    q.qn[ip, 7] = 0.0
                    q.qn[ip, end] = P

                    q.qe[ip, 1] = ρref
                    q.qe[ip, 2] = ρref*u
                    q.qe[ip, 3] = ρref*v
                    q.qe[ip, 4] = ρref*w
                    q.qe[ip, 5] = ρref*hl_ref
                    q.qe[ip, 6] = ρref*qt_ref
                    q.qe[ip, 7] = 0.0
                    q.qe[ip, end] = P
                end
            end
        end

        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                ρtot = q.qn[:, 1] .+ q.qe[:, 1]
                q.qn[:, 2] ./= ρtot
                q. qn[:, 3] ./= ρtot
                q.qn[:, 4] ./= ρtot
                q.qn[:, 5] ./= ρtot
                q.qe[:, 5] ./= q.qe[:, 1]
            else
                q.qn[:, 2] ./= q.qn[:, 1]
                q.qn[:, 3] ./= q.qn[:, 1]
                q.qn[:, 4] ./= q.qn[:, 1]
                q.qn[:, 5] ./= q.qn[:, 1]
                q.qn[:, 6] ./= q.qn[:, 1]
                q.qn[:, 7] ./= q.qn[:, 1]
                q.qe[:, 5] ./= q.qe[:, 1]
            end
        end
    end

    if rank == 0
        @info " Initialize fields for 3D BOMEX with full RT coupling ........................ DONE "
    end

    return q
end

function initialize_bomex!(z, param_set)
    PhysConst = PhysicalConst{Float64}()
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
        # Compute energy contributions
       
        e_kin::FT = FT(1 // 2) * (u^2 + v^2 + w^2)
        e_pot::FT = _grav * z
        e_tot = TD.total_energy(param_set, TS, e_kin, e_pot)
        hl_tot = T * PhysConst.cp  + PhysConst.g * z - PhysConst.Lc * q_pt.liq
        #@info hl_tot, e_tot, T, q_pt.liq
        bhasCondense = TD.has_condensate(param_set, TS)
        # println(TS)
        Δθ = FT(0.0)
        Δqt = FT(0.0)
        if z < FT(400)
            Δθ = FT((rand()-0.5) * e_tot / 100)
            Δqt = FT((rand()-0.5) * q_tot / 100)
        end
    end
    return ρ, u, v, w, hl_tot, Δθ, P, q_tot, Δqt, q_pt.liq
    # if bhasCondense
    #     println(z, ρ, u, T, θ_liq, θ, q_tot, q_pt.liq, P, bhasCondense)
    # end
end