using SparseArrays

function DSS_laplace_sparse!(SD::NSD_2D, Lel, ω, mesh, metrics, N, T; llump=false)
    """
    Assembles a sparse matrix for the Laplace operator and returns non-zero count and full matrix size.

    Args:
        mesh: A mesh struct (nelem, ngl, connijk).
        metrics: A struct (dξdx, dydη, dηdy, dxdξ).
        Lel: Element stiffness matrices.
        ω: Gauss-Lobatto quadrature weights.

    Returns:
        A sparse matrix, number of non-zeros, and full matrix size.
    """

    num_nodes = maximum(mesh.connijk)
    rows = Int[]
    cols = Int[]
    data = Float64[]

    for iel in 1:mesh.nelem
        for i in 1:mesh.ngl
            for j in 1:mesh.ngl
                ip = mesh.connijk[iel, i, j]

                # First term (dξdx * Lel * ω * dydη)
                for k in 1:mesh.ngl
                    jp = mesh.connijk[iel, k, j]
                    value = metrics.dξdx[iel, i, k] * Lel[i, k] * ω[j] * metrics.dydη[iel, i, k]
                    push!(rows, ip)
                    push!(cols, jp)
                    push!(data, value)
                end

                # Second term (dηdy * Lel * ω * dxdξ)
                for l in 1:mesh.ngl
                    jp = mesh.connijk[iel, i, l]
                    value = metrics.dηdy[iel, i, l] * Lel[j, l] * ω[i] * metrics.dxdξ[iel, i, l]
                    push!(rows, ip)
                    push!(cols, jp)
                    push!(data, value)
                end
            end
        end
    end

    # Create a sparse matrix in CSR format
    L_sparse = sparse(rows, cols, data, num_nodes, num_nodes)
    num_nonzeros = length(data)
    full_matrix_size = num_nodes * num_nodes
    
    #println("\nNumber of Non-zeros: ", num_nonzeros)
    #println("Full Matrix Size (including zeros): ", full_matrix_size)
    #println("Sparse Matrix Size (non-zeros only): ", sizeof(L_sparse)) # Approximate sparse size

    return L_sparse, num_nonzeros
end
#=

# Example Usage (replace with your actual data):

# Example mesh
mutable struct Mesh
    nelem::Int
    ngl::Int
    connijk::Array{Int, 3}
end

# Example metrics
mutable struct Metrics
    dξdx::Array{Float64, 3}
    dydη::Array{Float64, 3}
    dηdy::Array{Float64, 3}
    dxdξ::Array{Float64, 3}
end

function build_conn!(connijk)
    node_counter = 1
    for iel in 1:nelem
        for j in 1:ngl
            for i in 1:ngl
                connijk[iel, i, j] = node_counter
                node_counter += 1
            end
        end
    end
end

# Example parameters
nelem = 2
ngl = 3
connijk = Array{Int, 3}(undef, nelem, ngl, ngl)

dξdx = rand(nelem, ngl, ngl)
dydη = rand(nelem, ngl, ngl)
dηdy = rand(nelem, ngl, ngl)
dxdξ = rand(nelem, ngl, ngl)

build_conn!(connijk)
metrics = Metrics(dξdx, dydη, dηdy, dxdξ)
Lel = rand(ngl, ngl)
ω = rand(ngl)

mesh = Mesh(nelem, ngl, connijk)

L_sparse, num_nonzeros, full_matrix_size = DSS_laplace_sparse!(mesh, metrics, Lel, ω)

println("Sparse Matrix:")
println(L_sparse)

println("\nNumber of Non-zeros: ", num_nonzeros)
println("Full Matrix Size (including zeros): ", full_matrix_size)
println("Sparse Matrix Size (non-zeros only): ", sizeof(L_sparse)) # Approximate sparse size
=#
