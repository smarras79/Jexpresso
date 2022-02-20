#
# This file contains all the structs definitions
# S. Marras, Feb 2022
# 
export St_legendre, St_lgl
export f_lgl

struct St_legendre{TFloat}
    
    """
        struct St_legendrea{TFloat<:Real}
            legendre  :: TFloat
            dlegendre :: TFloat
            q         :: TFloat # q  = legendre(p+1)  - legendre(p-1)
            dq        :: TFloat # dq = dlegendre(p+1) - dlegendre(p-1)
        end
    """
    
    legendre  :: TFloat
    dlegendre :: TFloat
    q         :: TFloat # q  = legendre(p+1)  - legendre(p-1)
    dq        :: TFloat # dq = dlegendre(p+1) - dlegendre(p-1)  
end

#=function (l::St_legendre)()
    
    
   return [l.legendre,l.dlegendre]
end
=#


struct St_lgl{TFloat}
    lgl    :: TFloat
    weight :: TFloat
end


function f_lgl(lgl::St_lgl, p::Int8, TFloat::Float64)

"""
2  Evaluate recursion, the Legendre polynomial of order p
  and its Derivatives at coordinate x
 
  L_{p}  --> legendre of order p
  L'_{p} --> dlegendr of order p
  
  and
 
  q  = L_{p+1}  -  L_{p-1}
  q' = L'_{p+1} -  L'_{p-1} 
  
  Algorithm 22+24 of Kopriva's book
  
"""

    size    = p + 1;
    ksi     = zeros(TFloat, size)
    weights = zeros(TFloat, size)
    
    #LGL nodes
    LegendreGaussLobattoNodesAndWeights(lgl, nop);
    
    #LG nodes
    #LegendreGaussNodesAndWeights(lgl, nop);

    return lgl;
end

