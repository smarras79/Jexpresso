function check_length(V::Array, v_length::Int64, function_name::String)
    
    if length(V) != v_length
        println(" ERROR in ", function_name, ": the length of V is not as it should be.")
        println("          Check ", function_name)
        error(" ----- ")
    end
    return nothing
end
