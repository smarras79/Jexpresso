module build_lgl

using basis_structs

export lgl, 

function lgl(p::int, ksi::T, w::T)

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
    ksi     = zeros(T, size)
    weights = zeros(T, size)
    
    #LGL nodes
    LegendreGaussLobattoNodesAndWeights(lgl, nop);
    
    #LG nodes
    #LegendreGaussNodesAndWeights(lgl, nop);

    return lgl;
}


end
