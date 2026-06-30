using Random

function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 3D Lasher-Trapp 2001 with full RT coupling ............... "
    end

    #---------------------------------------------------------------------------------
    # Variables: static energy + two moisture scalars (SAM microphysics layout)
    #   hl  = liquid-water static energy = cp T + g z  (per unit mass)
    #   qt  = total water mixing ratio
    #   qp  = precipitating water (rain + snow + graupel, diagnosed by SAM micro)
    #---------------------------------------------------------------------------------
    qvars    = ("ρ", "ρu", "ρv", "ρw", "hl", "ρqt", "ρqp")
    qoutvars = ["ρ", "ρu", "ρv", "ρw", "hl", "ρqt", "ρqp",
                "T", "qn", "qc", "qi", "qr", "qs", "qg", "qsatt"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat,
                 inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    if inputs[:backend] == CPU()
        PhysConst = PhysicalConst{Float64}()

        if inputs[:lrestart] == true
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path],
                                     inputs, mesh.npoin, HDF5(); nvar=length(qvars))
        else
            # Sounding: 22 Jul 1995 1545 UTC (Lasher-Trapp et al. 2001, Fig. 1)
            # Format: z(m)  P(Pa)  T(°C)  qv(g/kg)  u(m/s)  v(m/s)
            data    = read_sounding(inputs[:sounding_file])
            # interpolate_sounding returns columns 2:end of the sounding at each mesh node
            background = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.z, data)
            # background[:,1] = P (Pa)
            # background[:,2] = T (°C)
            # background[:,3] = qv (g/kg)
            # background[:,4] = u (m/s)
            # background[:,5] = v (m/s)

            for ip = 1:mesh.npoin
                z = mesh.z[ip]

                pref   = background[ip, 1]
                T_ref  = background[ip, 2] + 273.15   # °C → K
                qv_ref = background[ip, 3] / 1000     # g/kg → kg/kg
                u_ref  = background[ip, 4]
                v_ref  = background[ip, 5]

                Tv_ref = T_ref * (1 + 0.61 * qv_ref)
                ρref   = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv_ref, Press=pref)
                hl_ref = PhysConst.cp * T_ref + PhysConst.g * z
                pref_m = ρref * PhysConst.Rair * T_ref + ρref * qv_ref * PhysConst.Rvap * T_ref

                ρ  = ρref
                hl = hl_ref

                # Random velocity perturbations below 1 km (Lasher-Trapp et al. 2001)
                Δu = 0.0#(z <= 1000.0) ? (rand() - 0.5) * 1.0 : 0.0   # uniform in [-0.5, 0.5] m/s
                Δv = 0.0#(z <= 1000.0) ? (rand() - 0.5) * 1.0 : 0.0

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip, 1] = ρ - ρref
                    q.qn[ip, 2] = ρ * (u_ref + Δu) - ρref * u_ref
                    q.qn[ip, 3] = ρ * (v_ref + Δv) - ρref * v_ref
                    q.qn[ip, 4] = 0.0
                    q.qn[ip, 5] = ρ * hl - ρref * hl_ref
                    q.qn[ip, 6] = ρ * qv_ref - ρref * qv_ref
                    q.qn[ip, 7] = 0.0
                    q.qn[ip, end] = pref_m

                    q.qe[ip, 1] = ρref
                    q.qe[ip, 2] = ρref * u_ref
                    q.qe[ip, 3] = ρref * v_ref
                    q.qe[ip, 4] = 0.0
                    q.qe[ip, 5] = ρref * hl_ref
                    q.qe[ip, 6] = ρref * qv_ref
                    q.qe[ip, 7] = 0.0
                    q.qe[ip, end] = pref_m
                else
                    q.qn[ip, 1] = ρ
                    q.qn[ip, 2] = ρ * (u_ref + Δu)
                    q.qn[ip, 3] = ρ * (v_ref + Δv)
                    q.qn[ip, 4] = 0.0
                    q.qn[ip, 5] = ρ * hl
                    q.qn[ip, 6] = ρ * qv_ref
                    q.qn[ip, 7] = 0.0
                    q.qn[ip, end] = pref_m

                    q.qe[ip, 1] = ρref
                    q.qe[ip, 2] = ρref * u_ref
                    q.qe[ip, 3] = ρref * v_ref
                    q.qe[ip, 4] = 0.0
                    q.qe[ip, 5] = ρref * hl_ref
                    q.qe[ip, 6] = ρref * qv_ref
                    q.qe[ip, 7] = 0.0
                    q.qe[ip, end] = pref_m
                end
            end
        end

        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                ρtot = q.qn[:, 1] .+ q.qe[:, 1]
                q.qn[:, 2] ./= ρtot
                q.qn[:, 3] ./= ρtot
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
        @info " Initialize fields for 3D Lasher-Trapp 2001 with full RT coupling ............... DONE "
    end

    return q
end
