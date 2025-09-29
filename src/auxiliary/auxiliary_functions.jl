function get_memory_usage(position_string::String)

    GC.gc()  # Force garbage collection first for more accurate reading
    stats = Base.gc_num()

    bytes = stats.allocd
    mib =  bytes / (1024^2)
    
    println(RED_FG(string(" ----- Total memory allocated so far ($position_string): $mib MB")))
    
    #return mib  # Total allocated memory in bytes
end


function remove_arrays!(matrix::Matrix{Int64}, rows_to_remove::Matrix{Int64})
    """
        Removes specific rows from a matrix (Matrix{Int64}).

        Args:
            matrix: The input matrix (Matrix{Int64}).
            rows_to_remove: A Matrix{Int64} of rows to remove.

        Returns:
            The modified matrix with the specified rows removed (in-place).
        """
    rows_to_keep = trues(size(matrix, 1))

    for i in 1:size(matrix, 1)
        for j in 1:size(rows_to_remove, 1)
            if matrix[i, :] == rows_to_remove[j, :]
                rows_to_keep[i] = false
                break
            end
        end
    end

    # Create a new matrix with the rows to keep
    new_matrix = matrix[rows_to_keep, :]
    
    return new_matrix
end

function replace_shared_values!(A::Matrix{Int}, B::Matrix{Int})
    # Create a set of all values in B for fast lookup.
    b_values = Set{Int}()
    for row in eachrow(B)
        union!(b_values, row)
    end

    # Iterate through A and replace shared values with -1.
    for i in eachindex(A)
        if A[i] in b_values
            A[i] = -1
        end
    end
    return A
end

function unroll_positive_unique(A::Matrix{<:Int})
    # Preallocate a BitSet for tracking seen values, which is efficient for unique checks.
    seen = BitSet()
    # Preallocate a vector with an estimated size. This avoids excessive reallocations.
    result = Vector{eltype(A)}()
    sizehint!(result, length(A)) # Initial guess, can be refined based on expectations.

    for val in A
        if val > 0 && !(val in seen)
            push!(result, val)
            push!(seen, reinterpret(UInt, val)) # Efficiently store and check floats
        end
    end

    return result
end

# Function to print the matrix in the desired symbolic format
function print_symbolic_matrix(matrix::Matrix{Int})
    rows, cols = size(matrix)  # Get the dimensions of the matrix

    for i in 1:rows
        for j in 1:cols
            # Print "A" followed by the row and column indices.
            # Use @printf for formatted output to ensure consistent spacing.
            @printf "A%d,%d " i j
        end
        println() # Move to the next line after each row is printed.
    end
end


function get_array_pointer(arr)
    return pointer(arr)
end

function norm2(X::SMatrix{3, 3, FT}) where {FT}
    abs2(X[1, 1]) +
    abs2(X[2, 1]) +
    abs2(X[3, 1]) +
    abs2(X[1, 2]) +
    abs2(X[2, 2]) +
    abs2(X[3, 2]) +
    abs2(X[1, 3]) +
    abs2(X[2, 3]) +
    abs2(X[3, 3])
end

# ### [Pricipal Invariants](@id tensor-invariants)
# ```math
# \textit{I}_{1} = \mathrm{tr(X)} \\
# \textit{I}_{2} = (\mathrm{tr(X)}^2 - \mathrm{tr(X^2)}) / 2 \\
# \textit{I}_{3} = \mathrm{det(X)} \\
# ```
"""
    principal_invariants(X)

Calculates principal invariants of a tensor `X`. Returns 3 element tuple containing the invariants.
"""
function principal_invariants(X)
    first = tr(X)
    second = (first^2 - tr(X .^ 2)) / 2
    third = det(X)
    return (first, second, third)
end

# ### [Symmetrize](@id symmetric-tensors)
# ```math
# \frac{\mathrm{X} + \mathrm{X}^{T}}{2} \\
# ```
"""
    symmetrize(X)

Given a (3,3) second rank tensor X, compute `(X + X')/2`, returning a
symmetric `SHermitianCompact` object.
"""
function symmetrize(X)
    SVector(
        X[1, 1],
        (X[2, 1] + X[1, 2]) / 2,
        (X[3, 1] + X[1, 3]) / 2,
        X[2, 2],
        (X[3, 2] + X[2, 3]) / 2,
        X[3, 3],
    )
end


# Given a tensor X, return the tensor dot product
# ```math
# \sum_{i,j} S_{ij}^2
# ```
"""
        norm2(X)
    Given a tensor `X`, computes X:X.
    """
function norm2(X)
        abs2(X[1, 1]) +
        abs2(X[2, 1]) +
        abs2(X[3, 1]) +
        abs2(X[1, 2]) +
        abs2(X[2, 2]) +
        abs2(X[3, 2]) +
        abs2(X[1, 3]) +
        abs2(X[2, 3]) +
        abs2(X[3, 3])
end
function norm2(X::SHermitianCompact{3, FT, 6}) where {FT}
        abs2(X[1, 1]) +
        2 * abs2(X[2, 1]) +
        2 * abs2(X[3, 1]) +
        abs2(X[2, 2]) +
        2 * abs2(X[3, 2]) +
        abs2(X[3, 3])
end
