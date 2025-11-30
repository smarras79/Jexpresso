"""
    ERA5 LES - Initialize

    Initialize simulation fields using ERA5 reanalysis data.

    Author: Jexpresso Development Team
    Date: 2025-11-30
"""

function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """
    Initialize 3D LES simulation with ERA5 data for real weather forecasting.
    """
    @info " Initialize fields for ERA5 LES case ........................ "

    #---------------------------------------------------------------------------------
    # Solution variables for compressible Euler with θ equation and moisture
    #---------------------------------------------------------------------------------
    qvars = ["ρ", "ρu", "ρv", "ρw", "ρθ", "ρqt", "ρql"]
    qoutvars = ["ρ", "u", "v", "w", "θ", "qt", "ql", "P"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    #---------------------------------------------------------------------------------
    # Load ERA5 data if enabled
    #---------------------------------------------------------------------------------
    if inputs[:lERA5]
        @info " Loading ERA5 initial conditions..."

        # Allocate ERA5 fields structure
        era5_fields = allocate_ERA5Fields(mesh.npoin, mesh, inputs, TFloat,
                                         inputs[:backend]; lERA5=true)

        # Use ERA5 data to initialize the simulation
        initialize_from_era5!(q, era5_fields, mesh, inputs, TFloat)

    else
        @info " ERA5 disabled - using default initialization"
        # Fallback to standard initialization
        initialize_default!(q, mesh, inputs, TFloat)
    end

    return q
end

"""
    initialize_from_era5!(q, era5_fields, mesh, inputs, TFloat)

Initialize simulation state from ERA5 data.
"""
function initialize_from_era5!(q, era5_fields, mesh, inputs, TFloat)

    PhysConst = PhysicalConst{Float64}()

    @info " Initializing from ERA5 data..."
    @info "   Number of mesh points: $(mesh.npoin)"

    # Constants
    R_dry = 287.05    # Gas constant for dry air [J/(kg·K)]
    p0 = 100000.0     # Reference pressure [Pa]
    κ = R_dry / 1004.0  # Poisson constant

    # Initialize each mesh point with ERA5 data
    for ip = 1:mesh.npoin

        x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]

        # Get ERA5 fields at this point (already interpolated to mesh)
        T = era5_fields.temperature[ip]
        u = era5_fields.u_wind[ip]
        v = era5_fields.v_wind[ip]
        w = era5_fields.w_wind[ip]
        qt = era5_fields.specific_humidity[ip]
        P = era5_fields.pressure[ip]

        # Compute potential temperature: θ = T * (p0/P)^κ
        θ = T * (p0 / P)^κ

        # Compute density: ρ = P / (R * T)
        ρ = P / (R_dry * T)

        # Initial cloud liquid water (assume zero, will be adjusted by saturation)
        ql = 0.0

        # Add small random perturbations to trigger turbulence
        # Typical LES practice: add ~0.1 K perturbations
        Random.seed!(ip)
        θ_pert = θ + 0.1 * (rand() - 0.5)
        qt_pert = qt + 1.0e-5 * (rand() - 0.5)

        # Store reference state (unperturbed ERA5 data)
        ρref = ρ
        θref = θ
        qtref = qt
        Pref = P

        # Solution type: perturbation or total
        if inputs[:SOL_VARS_TYPE] == PERT()
            # Perturbation formulation
            q.qn[ip,1] = 0.0                    # ρ' (no density perturbation initially)
            q.qn[ip,2] = 0.0                    # ρu' (no momentum perturbation)
            q.qn[ip,3] = 0.0                    # ρv'
            q.qn[ip,4] = 0.0                    # ρw'
            q.qn[ip,5] = ρref * (θ_pert - θref) # ρθ' (small temperature perturbation)
            q.qn[ip,6] = ρref * (qt_pert - qtref) # ρqt'
            q.qn[ip,7] = 0.0                    # ρql'
            q.qn[ip,end] = P

            # Store background state
            q.qe[ip,1] = ρref
            q.qe[ip,2] = u
            q.qe[ip,3] = v
            q.qe[ip,4] = w
            q.qe[ip,5] = ρref * θref
            q.qe[ip,6] = ρref * qtref
            q.qe[ip,7] = 0.0
            q.qe[ip,end] = Pref

        else
            # Total formulation
            q.qn[ip,1] = ρ
            q.qn[ip,2] = ρ * u
            q.qn[ip,3] = ρ * v
            q.qn[ip,4] = ρ * w
            q.qn[ip,5] = ρ * θ_pert
            q.qn[ip,6] = ρ * qt_pert
            q.qn[ip,7] = ρ * ql
            q.qn[ip,end] = P

            # Store reference state
            q.qe[ip,1] = ρref
            q.qe[ip,2] = ρref * u
            q.qe[ip,3] = ρref * v
            q.qe[ip,4] = ρref * w
            q.qe[ip,5] = ρref * θref
            q.qe[ip,6] = ρref * qtref
            q.qe[ip,7] = 0.0
            q.qe[ip,end] = Pref
        end
    end

    @info " ERA5 initialization complete!"
    @info "   Temperature range: $(minimum(era5_fields.temperature)) - $(maximum(era5_fields.temperature)) K"
    @info "   Wind speed range: $(minimum(sqrt.(era5_fields.u_wind.^2 .+ era5_fields.v_wind.^2))) - $(maximum(sqrt.(era5_fields.u_wind.^2 .+ era5_fields.v_wind.^2))) m/s"
end

"""
    initialize_default!(q, mesh, inputs, TFloat)

Default initialization (fallback if ERA5 is disabled).
"""
function initialize_default!(q, mesh, inputs, TFloat)

    @warn " Using default initialization - idealized atmosphere"

    PhysConst = PhysicalConst{Float64}()

    # Simple isothermal atmosphere
    T0 = 300.0  # K
    P0 = 100000.0  # Pa
    R_dry = 287.05
    g = 9.81
    H = R_dry * T0 / g  # Scale height

    for ip = 1:mesh.npoin
        x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]

        # Exponential atmosphere
        P = P0 * exp(-z / H)
        T = T0
        θ = T * (P0 / P)^(R_dry / 1004.0)
        ρ = P / (R_dry * T)

        # No wind
        u = 0.0
        v = 0.0
        w = 0.0

        # Dry atmosphere
        qt = 0.0
        ql = 0.0

        q.qn[ip,1] = ρ
        q.qn[ip,2] = ρ * u
        q.qn[ip,3] = ρ * v
        q.qn[ip,4] = ρ * w
        q.qn[ip,5] = ρ * θ
        q.qn[ip,6] = ρ * qt
        q.qn[ip,7] = ρ * ql
        q.qn[ip,end] = P

        q.qe[ip,:] .= q.qn[ip,:]
    end
end
