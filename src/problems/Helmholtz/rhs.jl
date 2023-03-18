include("../AbstractProblems.jl")

include("../../kernel/abstractTypes.jl")
include("../../kernel/mesh/mesh.jl")
include("../../kernel/mesh/metric_terms.jl")
include("../../io/print_matrix.jl")
include("./user_flux.jl")
include("./user_source.jl")

function rhs!(du, u, params, time)
    
    T       = params.T
    SD      = params.SD
    QT      = params.QT
    PT      = params.PT
    BCT     = params.BCT
    neqns   = params.neqns
    basis   = params.basis
    mesh    = params.mesh
    metrics = params.metrics
    inputs  = params.inputs
    ω       = params.ω
    M       = params.M
    De      = params.De
    Le      = params.Le

    RHS = build_rhs(SD, QT, PT, BCT, u, neqns, basis.ψ, basis.dψ, ω, mesh, metrics, M, De, Le, time, inputs, T)    
    du .= RHS
    
    return du #This is already DSSed
end

function build_rhs_source(SD::NSD_2D,
                          QT::Inexact,
                          PT::Elliptic,
                          q::Array,
                          mesh::St_mesh,
                          M::AbstractArray, #M is a vector for inexact integration
                          T)

    S = user_source(q, mesh, T)
    
    return M.*S    
end
