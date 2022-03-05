module Mod_operators_fd
export St_fd_derivatives
export ∑

struct St_fd_derivatives{T}

    #1st dertivatives
    ∂f∂x::T
    ∂f∂y::T
    ∂f∂z::T

    #2nd dertivatives
    ∂2f∂x2::T
    ∂2f∂y2::T
    ∂2f∂z2::T
    
end

const ∑ = sum

end #module
