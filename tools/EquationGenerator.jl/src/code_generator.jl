"""
Code generator for creating Julia problem directories from equation information.
"""

using Mustache

"""
    generate_code(eq_info::NamedTuple, output_dir::String, category::Union{String,Nothing}=nothing) -> String

Generate a complete problem directory with all Julia files.

# Arguments
- `eq_info::NamedTuple`: Equation information from analyze_equations
- `output_dir::String`: Base output directory (e.g., "problems")
- `category::Union{String,Nothing}`: Category subdirectory (e.g., "CompEuler")

# Returns
- Path to the created problem directory

# Example
```julia
problem_dir = generate_code(eq_info, "../../problems", "CompEuler")
```
"""
function generate_code(
    eq_info::NamedTuple,
    output_dir::String,
    category::Union{String,Nothing}=nothing
)
    # Determine output path
    if isnothing(category)
        category = replace(eq_info.equation_type, " " => "", "-" => "")
        category = uppercasefirst(category)
    end

    problem_name = eq_info.problem_name
    full_path = joinpath(output_dir, category, problem_name)

    # Create directory
    mkpath(full_path)

    # Prepare template context
    context = prepare_context(eq_info)

    # Generate each file
    files = Dict(
        "initialize.jl" => INITIALIZE_TEMPLATE,
        "user_bc.jl" => USER_BC_TEMPLATE,
        "user_flux.jl" => USER_FLUX_TEMPLATE,
        "user_inputs.jl" => USER_INPUTS_TEMPLATE,
        "user_primitives.jl" => USER_PRIMITIVES_TEMPLATE,
        "user_source.jl" => USER_SOURCE_TEMPLATE
    )

    for (filename, template) in files
        generate_file(filename, template, context, full_path)
    end

    # Generate README
    generate_readme(eq_info, full_path)

    return full_path
end

"""
    prepare_context(eq_info::NamedTuple) -> Dict

Prepare template context from equation information.
"""
function prepare_context(eq_info::NamedTuple)
    # Format variables list for Julia
    qvars_list = join(["\"$var\"" for var in eq_info.variables], ", ")

    # Determine output variables
    if !isempty(eq_info.primitive_variables)
        qoutvars = eq_info.primitive_variables
    else
        qoutvars = eq_info.variables
    end
    qoutvars_list = join(["\"$var\"" for var in qoutvars], ", ")

    # Format viscosity coefficients
    num_eqs = eq_info.num_equations
    viscosity = "[" * join(fill("0.0", num_eqs), ", ") * "]"

    # Create indexed variables for templates
    variables_indexed = [
        Dict("index" => i, "name" => var)
        for (i, var) in enumerate(eq_info.variables)
    ]

    flux_x_indexed = [
        Dict("index" => i, "flux" => flux)
        for (i, flux) in enumerate(eq_info.flux_x)
    ]

    flux_y_indexed = [
        Dict("index" => i, "flux" => flux)
        for (i, flux) in enumerate(eq_info.flux_y)
    ]

    source_indexed = [
        Dict("index" => i, "term" => term)
        for (i, term) in enumerate(eq_info.source_terms)
    ]

    # Format GPU flux returns
    flux_x_gpu = join(["T($flux)" for flux in eq_info.flux_x], ", ")
    flux_y_gpu = join(["T($flux)" for flux in eq_info.flux_y], ", ")
    source_gpu = join(["T($term)" for term in eq_info.source_terms], ", ")
    bc_gpu_vars = join(["T(qbdy[$i])" for i in 1:num_eqs], ", ")
    prim_gpu_vars = join(["T(u[$i])" for i in 1:num_eqs], ", ")

    # Check if has source terms
    has_source = any(s != "0" && s != "0.0" for s in eq_info.source_terms)

    context = Dict(
        "problem_name" => eq_info.problem_name,
        "description" => eq_info.description,
        "num_equations" => eq_info.num_equations,
        "spatial_dim" => eq_info.spatial_dimensions,
        "variables" => eq_info.variables,
        "variables_indexed" => variables_indexed,
        "qvars_list" => qvars_list,
        "qoutvars_list" => qoutvars_list,
        "flux_x" => eq_info.flux_x,
        "flux_x_indexed" => flux_x_indexed,
        "flux_y" => eq_info.flux_y,
        "flux_y_indexed" => flux_y_indexed,
        "source_terms" => eq_info.source_terms,
        "source_indexed" => source_indexed,
        "has_pressure" => eq_info.has_pressure,
        "has_gravity" => eq_info.has_gravity,
        "has_source" => has_source,
        "has_y_fluxes" => eq_info.spatial_dimensions >= 2,
        "lsource" => has_source ? "true" : "false",
        "viscosity_coefficients" => viscosity,
        "flux_x_gpu_pert" => flux_x_gpu,
        "flux_x_gpu_total" => flux_x_gpu,
        "flux_y_gpu_pert" => flux_y_gpu,
        "flux_y_gpu_total" => flux_y_gpu,
        "source_gpu" => source_gpu,
        "bc_gpu_vars" => bc_gpu_vars,
        "prim_gpu_vars" => prim_gpu_vars
    )

    return context
end

"""
    generate_file(filename::String, template::String, context::Dict, output_dir::String)

Generate a single Julia file from template.
"""
function generate_file(filename::String, template::String, context::Dict, output_dir::String)
    try
        content = Mustache.render(template, context)
        output_file = joinpath(output_dir, filename)

        open(output_file, "w") do f
            write(f, content)
        end

        println("✓ Generated $filename")
    catch e
        println("✗ Error generating $filename: $e")
        rethrow(e)
    end
end

"""
    generate_readme(eq_info::NamedTuple, output_dir::String)

Generate a README file for the problem.
"""
function generate_readme(eq_info::NamedTuple, output_dir::String)
    content = """
# $(eq_info.problem_name)

$(eq_info.description)

## Equation System

Number of equations: $(eq_info.num_equations)
Spatial dimensions: $(eq_info.spatial_dimensions)

### Variables
$(join(["- $var: $(get(eq_info.variable_descriptions, var, "N/A"))" for var in eq_info.variables], "\n"))

### Flux Terms

#### X-direction (F)
$(join(["$(i). $flux" for (i, flux) in enumerate(eq_info.flux_x)], "\n"))

$(if !isempty(eq_info.flux_y)
"#### Y-direction (G)\n" * join(["$(i). $flux" for (i, flux) in enumerate(eq_info.flux_y)], "\n")
else
""
end)

### Source Terms (S)
$(join(["$(i). $source" for (i, source) in enumerate(eq_info.source_terms)], "\n"))

## Notes

This directory was auto-generated by the Jexpresso Equation Generator tool (Julia version).

### TODO
- [ ] Review and verify flux calculations
- [ ] Implement initial conditions in `initialize.jl`
- [ ] Define boundary conditions in `user_bc.jl`
- [ ] Update `user_inputs.jl` with problem-specific parameters
- [ ] Verify primitive variable transformations in `user_primitives.jl`
- [ ] Test with actual mesh and verify results
"""

    readme_file = joinpath(output_dir, "README.md")
    open(readme_file, "w") do f
        write(f, content)
    end

    println("✓ Generated README.md")
end
