"""
Template definitions for generating Julia code files.
"""

# Template for initialize.jl
const INITIALIZE_TEMPLATE = """
function initialize(SD::NSD_{{spatial_dim}}D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    \"\"\"
    Initialize fields for {{spatial_dim}}D {{problem_name}}
    {{description}}
    \"\"\"
    @info " Initialize fields for {{spatial_dim}}D {{problem_name}} ........................ "

    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    #
    #---------------------------------------------------------------------------------
    qvars    = [{{qvars_list}}]
    qoutvars = [{{qoutvars_list}}]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------
    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()
        if (inputs[:case] === "default")
            # TODO: Define initial conditions here
            # This is a placeholder - modify according to your specific problem
            for iel_g = 1:mesh.nelem
                for j=1:mesh.ngl, i=1:mesh.ngl

                    ip = mesh.connijk[iel_g,i,j]
                    x, y = mesh.x[ip], mesh.y[ip]

                    # Initialize variables
{{#variables}}
                    {{.}}_val = 0.0  # TODO: Set initial value for {{.}}
{{/variables}}

                    if inputs[:SOL_VARS_TYPE] == PERT()
{{#variables_indexed}}
                        q.qn[ip,{{index}}] = {{name}}_val  # {{name}}
{{/variables_indexed}}
                        q.qn[ip,end] = 0.0  # pressure or additional variable
                    else
{{#variables_indexed}}
                        q.qn[ip,{{index}}] = {{name}}_val  # {{name}}
{{/variables_indexed}}
                        q.qn[ip,end] = 0.0  # pressure or additional variable
                    end
                end
            end

        else
            error(" ERROR: {{problem_name}}: initialize.jl:\\n assign value to inputs[:case]")
        end

        if inputs[:CL] == NCL()
            # TODO: Add non-conservative form adjustments if needed
        end

    else
        # GPU initialization
        if (inputs[:SOL_VARS_TYPE] == PERT())
            lpert = true
        else
            lpert = false
        end
        PhysConst = PhysicalConst{TFloat}()

        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, PhysConst, lpert; ndrange = (mesh.npoin))
    end
    @info " Initialize fields for {{spatial_dim}}D {{problem_name}} ........................ DONE "

    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, PhysConst, lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    x_val = x[ip]
    y_val = y[ip]

    # TODO: Define initial conditions for GPU
{{#variables}}
    {{.}}_val = T(0.0)  # TODO: Set initial value for {{.}}
{{/variables}}

    if (lpert)
{{#variables_indexed}}
        qn[ip,{{index}}] = {{name}}_val
{{/variables_indexed}}
        qn[ip,end] = T(0.0)
    else
{{#variables_indexed}}
        qn[ip,{{index}}] = {{name}}_val
{{/variables_indexed}}
        qn[ip,end] = T(0.0)
    end
end
"""

# Template for user_flux.jl
const USER_FLUX_TEMPLATE = """
function user_flux!(F, G, SD::NSD_{{spatial_dim}}D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs={{num_equations}}, ip=1)

    PhysConst = PhysicalConst{Float64}()

    # Extract conservative variables
{{#variables_indexed}}
    {{name}} = q[{{index}}]
{{/variables_indexed}}

{{#has_pressure}}
    # Calculate pressure (TODO: Adjust based on equation of state)
    Press = 0.0  # TODO: Define pressure calculation
{{/has_pressure}}

    # X-direction fluxes
{{#flux_x_indexed}}
    F[{{index}}] = {{flux}}
{{/flux_x_indexed}}

{{#has_y_fluxes}}
    # Y-direction fluxes
{{#flux_y_indexed}}
    G[{{index}}] = {{flux}}
{{/flux_y_indexed}}
{{/has_y_fluxes}}
end

function user_flux!(F, G, SD::NSD_{{spatial_dim}}D,
                    q,
                    qe,
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs={{num_equations}}, ip=1)

    PhysConst = PhysicalConst{Float64}()

    # Extract conservative variables (perturbation form)
{{#variables_indexed}}
    {{name}} = q[{{index}}] + qe[{{index}}]
{{/variables_indexed}}

{{#has_pressure}}
    # Calculate pressure
    Press = 0.0  # TODO: Define pressure calculation
    Press = Press - qe[end]
{{/has_pressure}}

    # X-direction fluxes
{{#flux_x_indexed}}
    F[{{index}}] = {{flux}}
{{/flux_x_indexed}}

{{#has_y_fluxes}}
    # Y-direction fluxes
{{#flux_y_indexed}}
    G[{{index}}] = {{flux}}
{{/flux_y_indexed}}
{{/has_y_fluxes}}
end

function user_flux_gpu(q, qe, PhysConst, lpert)
    T = eltype(q)

    if (lpert)
        # Perturbation form
{{#variables_indexed}}
        {{name}} = q[{{index}}] + qe[{{index}}]
{{/variables_indexed}}

{{#has_pressure}}
        Pressure = T(0.0)  # TODO: Define pressure calculation
{{/has_pressure}}

        return {{flux_x_gpu_pert}}{{#has_y_fluxes}}, {{flux_y_gpu_pert}}{{/has_y_fluxes}}
    else
        # Total form
{{#variables_indexed}}
        {{name}} = q[{{index}}]
{{/variables_indexed}}

{{#has_pressure}}
        Pressure = T(0.0)  # TODO: Define pressure calculation
{{/has_pressure}}

        return {{flux_x_gpu_total}}{{#has_y_fluxes}}, {{flux_y_gpu_total}}{{/has_y_fluxes}}
    end
end
"""

# Template for user_source.jl
const USER_SOURCE_TEMPLATE = """
function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs={{num_equations}}, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    PhysConst = PhysicalConst{Float64}()

    # Extract variables
{{#variables_indexed}}
    {{name}} = q[{{index}}]
{{/variables_indexed}}

    # Source terms
{{#source_indexed}}
    S[{{index}}] = {{term}}
{{/source_indexed}}

end

function user_source!(S,
                      q,
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs={{num_equations}}, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    PhysConst = PhysicalConst{Float64}()

    # Extract variables (perturbation form)
{{#variables_indexed}}
    {{name}} = q[{{index}}]  # Perturbation
{{/variables_indexed}}

    # Source terms
{{#source_indexed}}
    S[{{index}}] = {{term}}
{{/source_indexed}}

end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)

    T = eltype(q)

    # Extract variables
{{#variables_indexed}}
    {{name}} = q[{{index}}]
{{/variables_indexed}}

    # Source terms
    return {{source_gpu}}
end
"""

# Template for user_bc.jl
const USER_BC_TEMPLATE = """
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    # TODO: Define Dirichlet boundary conditions for TOTAL form
    # Example: slip boundary condition for momentum
    # qnl = nx*q[2] + ny*q[3]
    # qbdy[2] = q[2] - qnl*nx
    # qbdy[3] = q[3] - qnl*ny

{{#variables_indexed}}
    # qbdy[{{index}}] = ...  # {{name}}
{{/variables_indexed}}

end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx::AbstractFloat, ny::AbstractFloat, qe, ::PERT)
    # TODO: Define Dirichlet boundary conditions for PERT form

{{#variables_indexed}}
    # qbdy[{{index}}] = ...  # {{name}}
{{/variables_indexed}}

end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    # TODO: Define Neumann boundary conditions
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)
    T = eltype(q)
    # TODO: Define GPU boundary conditions
    if (lpert)
        return {{bc_gpu_vars}}
    else
        return {{bc_gpu_vars}}
    end
end
"""

# Template for user_primitives.jl
const USER_PRIMITIVES_TEMPLATE = """
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    # TODO: Convert conservative to primitive variables
{{#variables_indexed}}
    uprimitive[{{index}}] = u[{{index}}]  # TODO: {{name}}
{{/variables_indexed}}
end

function user_primitives!(u, qe, uprimitive, ::PERT)
    # TODO: Convert conservative to primitive variables (perturbation form)
{{#variables_indexed}}
    uprimitive[{{index}}] = u[{{index}}]  # TODO: {{name}}
{{/variables_indexed}}
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    if (lpert)
        return {{prim_gpu_vars}}
    else
        return {{prim_gpu_vars}}
    end
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)
    # TODO: Define output variables
{{#variables_indexed}}
    uout[{{index}}] = u[{{index}}]  # {{name}}
{{/variables_indexed}}
    uout[end] = u[end]
end
"""

# Template for user_inputs.jl
const USER_INPUTS_TEMPLATE = """
function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :tend                 => 1000.0,
        :ode_solver           => SSPRK54(),
        :Δt                   => 0.5,
        :diagnostics_at_times => [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
        :case                 => "default",
        :lsource              => {{lsource}},
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => false,
        :μ                   => {{viscosity_coefficients}},
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10.msh",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => true,
        :output_dir          => "./output",
        :loutput_pert        => true,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs

end
"""
