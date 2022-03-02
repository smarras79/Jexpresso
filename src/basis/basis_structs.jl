#
# This file contains all the structs definitions
# S. Marras, Feb 2022
# 
export St_legendre
export St_lgl
export build_lgl!
export LegendreGaussLobattoNodesAndWeights

mutable struct St_legendre{TFloat}
    
    """
        struct St_legendrea{TFloat<:Real}
            legendre  :: TFloat
            dlegendre :: TFloat
            q         :: TFloat -> q  = legendre(p+1)  - legendre(p-1)
            dq        :: TFloat -> dq = dlegendre(p+1) - dlegendre(p-1)
        end
    """
    
    legendre  :: TFloat
    dlegendre :: TFloat
    q         :: TFloat
    dq        :: TFloat
end

mutable struct St_lgl{TFloat}
    ξ::Array{TFloat}
    ω::Array{TFloat}
end


function build_lgl!(Legendre::St_legendre, lgl::St_lgl, nop::TInt)

"""
  Evaluate recursion, the Legendre polynomial of order p
  and its Derivatives at coordinate x
 
  L_{p}  --> legendre of order p
  L'_{p} --> dlegendr of order p
  
  and
 
  q  = L_{p+1}  -  L_{p-1}
  q' = L'_{p+1} -  L'_{p-1} 
  
  Algorithm 22+24 of Kopriva's book
  
"""
    size::Int8 = nop+1
    
    lgl.ξ = zeros(Float64, size)
    lgl.ω = zeros(Float64, size)
    
    #LGL nodes
    LegendreGaussLobattoNodesAndWeights(Legendre, lgl, nop);
    
    #LG nodes
    #LegendreGaussNodesAndWeights(lgl, nop);

    return lgl;
end


function LegendreGaussLobattoNodesAndWeights(Legendre::St_legendre, lgl::St_lgl, nop::TInt)
    
    NITER = 100
    TOL = 4*eps()
    
    ξ0 ::Float64=0.0
    ξ1 ::Float64=1.0
    ξP ::Float64=1.0
    ξj ::Float64=0.0
    ξj2::Float64=0.0
    L2 ::Float64=0.0
    
    ω0 ::Float64=0.0
    ω1 ::Float64=0.0
    ωP ::Float64=0.0
    
    Δ  ::Float64=0.0
    
    println( " # Compute LGL nodes ........................")
    
    for j=1:nop+1
	lgl.ξ[j] = 0.0;
	lgl.ω[j] = 1.0;
    end
    
    if (nop == 1)
	ξ0           = -1.0
	ω0           =  1.0
	ξ1           =  1.0
	ω1           =   ω0
        
	lgl.ξ[1]     =   ξ0
	lgl.ξ[nop+1] =   ξ1
	lgl.ω[1]     =   ω0
	lgl.ω[nop+1] =   ω1
    else 
	ξ0            = -1.0
	ω0            =  TFloat(2.0/(nop*(nop + 1)))
	ξP            =  1.0
	ωP            =   ω0
        
	lgl.ξ[1]      =   ξ0
	lgl.ξ[nop+1]  =   ξP
	lgl.ω[1]      =   ω0
	lgl.ω[nop+1]  =   ωP
	
        for jj = 2:floor(Int,(nop + 1)/2) 
	    j = jj - 1
	    ξj = -cos((j + 0.25)*π/nop - 3.0/(8.0*nop*π*(j + 0.25)))
	    lgl.ξ[jj] = ξj;
            
            for k = 0:NITER
	        LegendreAndDerivativeAndQ!(Legendre, nop, ξj)
	        Δ = -Legendre.q/Legendre.dq
	        ξj    =  ξj + Δ
	        
	        if (abs(Δ) <= TOL*abs(ξj))
                    break
                end
	    end
            LegendreAndDerivativeAndQ!(Legendre, nop, ξj)
	    lgl.ξ[jj]      =  ξj
	    lgl.ξ[nop+1-j] = -ξj
	    xj2            =  ξj*ξj
	    L2             = Legendre.legendre*Legendre.legendre
	    lgl.ω[jj]      = 2.0/(nop*(nop + 1.0)*L2)
	    lgl.ω[nop+1-j] = lgl.ω[jj]
                        
        end
    end
    
    if (mod(nop,2) == 0)
	LegendreAndDerivativeAndQ!(Legendre, nop, 0.0);
	lgl.ξ[TInt(nop/2)+1] = 0.0;
	
	L2           = Legendre.legendre*Legendre.legendre;
	lgl.ω[TInt(nop/2)+1] = 2.0/(nop*(nop + 1.0)*L2);
    end

    for j=1:nop+1       
        println( " # ξ, ω =: ", " ", lgl.ξ[j], " " , lgl.ω[j])
    end
    
    println(" # Compute LGL nodes ........................ DONE")
    
end

function LegendreAndDerivativeAndQ!(Legendre::St_legendre, nop::TInt, x::TFloat)
    
    """
         Evaluate by recursion, the Legendre polynomial of order p
         its Derivatives at coordinate x, and q:
 
          L_{p}  --> legendre of order p
          L'_{p} --> dlegendr of order p
          q  = L_{p+1}  -  L_{p-1}
          q' = L'_{p+1} -  L'_{p-1} 
 
          Algorithm 24 of Kopriva's book
 
          Simone Marras, October 2021
     """
    TFloat=Float64
    
    a   ::TFloat=0.0
    b   ::TFloat=0.0
    L   ::TFloat=0.0
    dL  ::TFloat=0.0
    Lp1 ::TFloat=0.0
    dLp1::TFloat=0.0
    Lm1 ::TFloat=0.0
    dLm1::TFloat=0.0
    Lm2 ::TFloat=0.0
    dLm2::TFloat=0.0
    
        
    #st_legendre Legendre;
        
    if (nop == 0)
	L  = 1.0
	dL = 0.0
    elseif (nop == 1)
	L  = x
	dL = 1.0
    else
	Lm2  = 1.0
	Lm1  = x
	dLm2 = 0.0
	dLm1 = 1.0
	
	#Construct Nth Order Legendre Polynomial
	for k=2:nop
            
	    a = TFloat((2.0*k - 1.0)/k)
	    b = TFloat((k - 1.0)/k)
            
	    L  = a*x*Lm1 - b*Lm2
	    dL = dLm2 + (2.0*k - 1.0)*Lm1
	    
	    a = TFloat((2.0*(k+1) - 1.0)/(k+1))
	    b = TFloat((k+1 - 1.0)/(k+1))
	    Lp1  = a*x*L - b*Lm1
	    dLp1 = dLm1 + (2.0*k - 1.0)*L
	    
	    Lm2 = Lm1
	    Lm1 = L

	    dLm2 = dLm1
	    dLm1 = dL
	end
    end
    
    Legendre.legendre  = L
    Legendre.dlegendre = dL
    
    Legendre.q  =  Lp1 -  Lm2 #WARNING: FIX THIS FOR nop=0 and nop=1
    Legendre.dq = dLp1 - dLm2 #WARNING: FIX THIS FOR nop=0 and nop=1
end



