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
