using Noise
"""
    Source for Helmholtz equation:

```math
    ∇²(q(x,y) - k²q(x,y) = W (1)
```
    where W is white noise.
    We solve (1)  on a circular domain with infinite extension.

"""
function user_source(q::Array, mesh::St_mesh, T)

    S = zeros(T, mesh.npoin)

    #
    # S(x,y) .= k²q(x,y) - W 
    #
    k2     = 1.0^2
    
    S .= k2*q .- add_gauss(q[:,1], 0.1)
    
    return  S
    
end
