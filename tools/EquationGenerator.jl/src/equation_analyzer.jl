"""
Equation Analyzer using Claude AI to extract structured equation information.
"""

using JSON3

"""
    analyze_equations(pdf_text::String, api_key::Union{String,Nothing}=nothing) -> NamedTuple

Analyze equations from PDF text and extract structured information using Claude AI.

# Arguments
- `pdf_text::String`: Text extracted from PDF
- `api_key::Union{String,Nothing}`: Anthropic API key

# Returns
A NamedTuple containing:
- `problem_name::String`: Suggested name for the problem
- `description::String`: Brief description
- `num_equations::Int`: Number of equations
- `spatial_dimensions::Int`: 1, 2, or 3
- `variables::Vector{String}`: Variable names
- `variable_descriptions::Dict`: Variable descriptions
- `flux_x::Vector{String}`: X-direction flux terms
- `flux_y::Vector{String}`: Y-direction flux terms
- `flux_z::Vector{String}`: Z-direction flux terms
- `source_terms::Vector{String}`: Source terms
- `primitive_variables::Vector{String}`: Primitive variables
- `has_pressure::Bool`: Whether system has pressure
- `has_gravity::Bool`: Whether system has gravity
- `equation_type::String`: Type of equation system

# Example
```julia
pdf_text = parse_pdf("equations.pdf")
eq_info = analyze_equations(pdf_text)
```
"""
function analyze_equations(pdf_text::String, api_key::Union{String,Nothing}=nothing)
    prompt = """You are analyzing a system of partial differential equations from a PDF document.
The equations are in the form:

∂_t(q) + ∂_x(F) + ∂_y(G) = S

where:
- q is the vector of conserved variables
- F is the flux in the x-direction
- G is the flux in the y-direction
- S is the source term vector

Here is the text extracted from the PDF:

$pdf_text

Please analyze this text and extract the following information in JSON format:

{
    "problem_name": "suggested_name_for_problem",
    "description": "brief description of the equation system",
    "num_equations": <number of equations>,
    "spatial_dimensions": <1, 2, or 3>,
    "variables": ["var1", "var2", ...],
    "variable_descriptions": {"var1": "description", ...},
    "flux_x": ["flux_x_1", "flux_x_2", ...],
    "flux_y": ["flux_y_1", "flux_y_2", ...],
    "flux_z": ["flux_z_1", "flux_z_2", ...],
    "source_terms": ["source_1", "source_2", ...],
    "primitive_variables": ["prim_var1", "prim_var2", ...],
    "has_pressure": true/false,
    "has_gravity": true/false,
    "equation_type": "euler/advection-diffusion/shallow-water/other"
}

For example, for the compressible Euler equations with theta and tracers:
- variables: ["ρ", "ρu", "ρv", "ρθ", "c1", "c2"]
- flux_x: ["ρu", "ρu²+p", "ρvu", "ρθu", "c1u", "c2u"]
- flux_y: ["ρv", "ρuv", "ρv²+p", "ρθv", "c1v", "c2v"]
- source_terms: ["0", "0", "-ρg", "0", "0", "0"]

Extract the information and respond with ONLY the JSON object, no additional text."""

    # Call Claude API
    response_text = call_claude(prompt, api_key; max_tokens=4096)

    # Clean up response (remove markdown code blocks if present)
    clean_text = strip(response_text)
    if startswith(clean_text, "```")
        lines = split(clean_text, '\n')
        clean_text = join(filter(l -> !startswith(strip(l), "```"), lines), '\n')
    end

    # Parse JSON
    try
        data = JSON3.read(clean_text)

        # Convert to NamedTuple with proper types
        eq_info = (
            problem_name = String(data.problem_name),
            description = String(get(data, :description, "")),
            num_equations = Int(data.num_equations),
            spatial_dimensions = Int(get(data, :spatial_dimensions, 2)),
            variables = String.(data.variables),
            variable_descriptions = Dict{String,String}(
                String(k) => String(v) for (k,v) in pairs(get(data, :variable_descriptions, Dict()))
            ),
            flux_x = String.(data.flux_x),
            flux_y = String.(get(data, :flux_y, String[])),
            flux_z = String.(get(data, :flux_z, String[])),
            source_terms = String.(data.source_terms),
            primitive_variables = String.(get(data, :primitive_variables, String[])),
            has_pressure = Bool(get(data, :has_pressure, false)),
            has_gravity = Bool(get(data, :has_gravity, false)),
            equation_type = String(get(data, :equation_type, "custom"))
        )

        return eq_info

    catch e
        error("Failed to parse Claude response as JSON: $e\nResponse: $clean_text")
    end
end
