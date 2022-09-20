#
# This file contains all the structs definitions
# S. Marras, Feb 2022
# 
export St_Lagrange
export St_Legendre
export St_lgl
export build_lgl!
export LegendreGaussLobattoNodesAndWeights

mutable struct St_Legendre{TFloat}
    
    """
    struct St_Legendrea{TFloat<:Real}
          legendre  :: TFloat
          dlegendre :: TFloat
           q        :: TFloat -> q  = legendre(p+1)  - legendre(p-1)
           dq       :: TFloat -> dq = dlegendre(p+1) - dlegendre(p-1)
    end
    """
    
    legendre  :: TFloat
    dlegendre :: TFloat
    q         :: TFloat
    dq        :: TFloat
end

abstract type AbstractIntegrationPointAndWeights end
abstract type AbstractInterpolationBasis end
abstract type AbstractSpaceDimensions end

struct LagrangeBasis <: AbstractInterpolationBasis end
struct NSD_1D <: AbstractSpaceDimensions end
struct NSD_2D <: AbstractSpaceDimensions end
struct NSD_3D <: AbstractSpaceDimensions end

mutable struct St_Chebyshev{TFloat} <:AbstractIntegrationPointAndWeights
    chebyshev ::TFloat
end

mutable struct St_lgl{TFloat} <:AbstractIntegrationPointAndWeights
    ξ::Array{TFloat}
    ω::Array{TFloat}
end

mutable struct St_lg{TFloat} <:AbstractIntegrationPointAndWeights
    ξ::Array{TFloat}
    ω::Array{TFloat}
end

mutable struct St_cg{TFloat} <:AbstractIntegrationPointAndWeights
    ξ::Array{TFloat}
    ω::Array{TFloat}
end

mutable struct St_cgl{TFloat} <:AbstractIntegrationPointAndWeights
    ξ::Array{TFloat}
    ω::Array{TFloat}
end

mutable struct St_Lagrange{TFloat} <:AbstractInterpolationBasis
    ψ::Matrix{TFloat}
    dψ::Matrix{TFloat}
end

function build_Integration_points!(::AbstractIntegrationPointAndWeights,nop::TInt) end
 
function build_Integration_points!(cg::St_cg, nop::TInt)
  build_cg!(cg,nop)
end

function build_Integration_points!(cgl::St_cgl,nop::TInt)
  build_cgl!(cgl,nop)
end

function build_Integration_points!(lg::St_lg,nop::TInt)
  Legendre = St_Legendre{TFloat}(0.0,0.0,0.0,0.0)
  build_lg!(Legendre,lg,nop)
end

function build_Integration_points!(lgl::St_lgl,nop::TInt)
  Legendre = St_Legendre{TFloat}(0.0,0.0,0.0,0.0)
  build_lgl!(Legendre,lgl,nop)
end

function build_Interpolation_basis!(TP::LagrangeBasis, SD::NSD_1D, T::Type{Float64}, ξ, ξq)

    Nξ = size(ξ,1)  - 1
    Qξ = size(ξq,1) - 1

    N  = (Nξ + 1)
    Q  = (Qξ + 1)
    
    basis = St_Lagrange{T}(zeros(N,Q), zeros(N,Q))
    (basis.ψ, basis.dψ) = LagrangeInterpolatingPolynomials_classic(ξ, ξq, T)
    
    return basis
end


function build_Interpolation_basis!(TP::LagrangeBasis, SD::NSD_2D, T::Type{Float64}, ξ, ξq)

    Nξ = size(ξ,1)  - 1
    Qξ = size(ξq,1) - 1
    
    Nη = Nξ
    Qη = Qξ

    N  = (Nξ + 1)*(Nη + 1)
    Q  = (Qξ + 1)*(Qη + 1)
    
    basis = St_Lagrange{T}(zeros(N,Q), zeros(N,Q))    
    (ψ, dψ) = LagrangeInterpolatingPolynomials_classic(ξ, ξq, T)
    
    for i = 1:Nξ+1
        ψ[i] .= basis.ψ[i,:]
        for j = 1:Nη+1
            l = i + 1 + j*(Nξ + 1)
            basis.ψ[l] = ψ[i]*ψ[j]
        end
    end
    
    return basis
end


function build_Interpolation_basis!(TP::LagrangeBasis, SD::NSD_3D, T::Type{Float64}, ξ, ξq)

    Nξ = size(ξ,1)  - 1
    Qξ = size(ξq,1) - 1
    
    Nη = Nξ #NOTICE THAT THESE SHOULD Be size(η,1) - 1 when we add different orders in different directions.
    Qη = Qξ # " " " 
    
    Nζ = Nξ #NOTICE THAT THESE SHOULD Be size(ζ,1) - 1 when we add different orders in different directions.
    Qζ = Qξ # " " " 

    N  = (Nξ + 1)*(Nη + 1)*(Nζ + 1)
    Q  = (Qξ + 1)*(Qη + 1)*(Qζ + 1)
    
    basis = St_Lagrange{T}(zeros(N,Q), zeros(N,Q))
    (ψ, dψ) = LagrangeInterpolatingPolynomials_classic(ξ, ξq, T)

    for i = 1:Nξ+1
        for j = 1:Nη+1
            for k = 1:Nζ+1
                l = i + 1 + j*(Nξ + 1) + k*(Nξ + 1)*(Nη + 1)
                basis.ψ[l] = ψ[i]*ψ[j]*ψ[k]
            end
        end
    end
    
    return basis
end



function build_lgl!(Legendre::St_Legendre, lgl::St_lgl, nop::TInt)

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
    LegendreGaussLobattoNodesAndWeights!(Legendre, lgl, nop);
    
    return lgl;
end


function build_cg!(cg::St_cg, nop::TInt)
    size::Int8=nop+1
    cg.ξ = zeros(Float64, size)
    cg.ω = zeros(Float64, size)

    #CG nodes
    ChebyshevGaussNodesAndWeights!(cg,nop)
    return cg
end

function build_cgl!(cgl::St_cgl, nop::TInt)
    size::Int8=nop+1
    cgl.ξ = zeros(Float64, size)
    cgl.ω = zeros(Float64, size)

    #CGL nodes
    ChebyshevGaussLobattoNodesAndWeights!(cgl,nop)
    return cgl
end

function build_lg!(Legendre::St_Legendre,lg::St_lg,nop)
    size::Int8=nop+1
    lg.ξ = zeros(Float64, size)
    lg.ω = zeros(Float64, size)

    #LG nodes
    LegendreGaussNodesAndWeights(Legendre,lg,nop)
    return lg
end
function ChebyshevGaussNodesAndWeights!(cg::St_lg, nop::TInt)
    """
         Compute the Nodes for the Chebyshev Gauss Quadrature
         using Algorithm 26 of Kopriva's book 
      """
    for j=0:nop
        cg.ξ[j+1]=-cospi((2*j+1)/(2*nop+2))
        cg.ω[j+1]=π/(nop+1)
    end
end

function ChebyshevGaussLobattoNodesAndWeights!(cgl::St_cgl,nop::TInt)
    for j=0:nop
        cgl.ξ[j+1]=-cospi(j/nop)
        cgl.ω[j+1]=π/nop
    end
    cgl.ω[1]=cgl.ω[1]/2
    cgl.ω[nop+1]=cgl.ω[nop+1]/2
end


function LegendreGaussNodesAndWeights!(Legendre::St_Legendre, lg::St_lg, nop::TInt)
    """
          Compute the Nodes and Weights for the Legendre Gauss Quadrature
          using Algorithm 23 of Kopriva's book valid for nop ≤  200
       """
    NITER = 100
    TOL = 4*eps(Float64)
    Δ::Float64=0.0
    if (nop == 0)
        lg.ξ[1]=0
        lg.ω[1]=2
    elseif (nop == 1)
        lg.ξ[1] = -sqrt(1/3)
        lg.ω[1] = 1
        lg.ξ[2] = - lg.ξ[1]
        lg.ω[2] = lg.ω[1]
    else
        for j=0:floor(Int,(nop+1)/2)-1
            lg.ξ[j+1] = - cospi((2*j+1)/(2*nop+2))
            for k=0:NITER
                LegendreAndDerivativeAndQ!(Legendre, nop+1, lg.ξ[j+1])
                Δ = -Legendre.legendre/Legendre.dlegendre
                lg.ξ[j+1]=lg.ξ[j+1]+Δ
                if (abs(Δ)≤TOL)
                    break
                end
            end
            LegendreAndDerivativeAndQ!(Legendre, nop+1, lg.ξ[j+1])
            lg.ξ[nop-j+1] = - lg.ξ[j+1]
            lg.ω[j+1]=2/(1-lg.ξ[j+1]^2)/(Legendre.dlegendre)^2
            lg.ω[nop-j+1]=lg.ω[j+1]
        end
    end
    if (mod(nop,2) == 0)
        LegendreAndDerivativeAndQ!(Legendre, nop+1, 0.0)
        lg.ξ[Tint(nop/2)] = 0
        lg.ω[Tint(nop/2)] = 2/Legendre.dlegendre^2
    end
end

function LegendreGaussLobattoNodesAndWeights!(Legendre::St_Legendre, lgl::St_lgl, nop::TInt)
    
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

function LegendreAndDerivativeAndQ!(Legendre::St_Legendre, nop::TInt, x::TFloat)
    
    """
             Evaluate by recursion, the Legendre polynomial of order p
             its Derivatives at coordinate x, and q:
     
              L_{p}  --> legendre of order p
              L'_{p} --> dlegendr of order p
              q  = L_{p+1}  -  L_{p-1}
              q' = L'_{p+1} -  L'_{p-1} 
     
              Algorithm 24 of Kopriva's book
     
              Simone Marras, October 2021

              
              Note that this algorithm looses precision at high enough nop
              if we intendt to use particularly large nop we should examine Yakimiw 1996 (Yassine) 

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
    
    if (nop == 0) #Order 0 case
	Legendre.legendre  = 1.0
	Legendre.dlegendre = 0.0
        Legendre.q = x
        Legendre.dq = 1.0
    elseif (nop == 1)
	Legendre.legendre  = x
	Legendre.dlegendre = 1.0
        Legendre.q = 0.5*(3*x^2-2)-1
        Legendre.dq = 3*x
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
        Legendre.legendre = L
        Legendre.dlegendre = dL
        Legendre.q = Lp1 - Lm2
        Legendre.dq = dLp1 - dLm2
    end
end

function ChebyshevPolynomial!(Chebyshev::St_Chebyshev,nop::TInt,x::TFloat,Ks::TInt)
    """
          Evaluate by recursion the Chebyshev Polynomial
          and switch to direct evaluation when nop > Ks
          uses Algorithm 21 of Kopriva's book, 
          NOTE Ks depends should never exceed 70 regardless of machine architecture
          Yassine Tissaoui March 2022
       """
    T2::TFloat=0.0
    T1::TFloat=0.0
    T::TFloat=0.0
    if (nop == 0) #case for order 0
        Chebyshev.chebyshev = 1
    elseif (nop ==1) #case for order 1
        Chebyshev.chebyshev = x
    elseif (nop ≤ Ks)
        T2=1
        T1=x
        for j=2:nop
            T=2*x*T1 - T2
            T2=T1
            T1=T
        end
        Chebyshev.chebyshev = T
    else
        T=cos(nop*acos(x))
        Chebyshev.chebyshev = T
    end
end

"""
    LagrangeInterpolatingPolynomials_classic(ξ, ξq, N, Q, TFloat)
    ξ::set of N interpolation points (e.g. LGL points)
    ξq::point to interpolate to (e.g. quadrature points of points within the element)
    
    Algorithm 3.1 + 3.2 from Giraldo's book

    from https://github.com/fxgiraldo/Element-based-Galerkin-Methods/blob/master/Projects/Project_01_1D_Interpolation/For_Instructors/julia/lagrange_basis.jl

"""
function LagrangeInterpolatingPolynomials_classic(ξ, ξq, TFloat)

    N = size(ξ,1) - 1
    Q = size(ξq,1) - 1
    
    #Initialize arrays
    L    = zeros(TFloat, N+1, Q+1)
    dLdx = zeros(TFloat, N+1, Q+1)
 
    for l=1:Q+1
        xl = ξq[l]

        #Construct Basis
        for i=1:N+1
            
            xi        = ξ[i]
            L[i,l]    = 1.0
            dLdx[i,l] = 0.0
            for j=1:N+1
                xj = ξ[j]

                #L
                if (j != i)
                    L[i,l] = L[i,l]*(xl - xj)/(xi - xj)
                end
                
                ddL=1
                if (j != i)
                    for k=1:N+1
                        xk = ξ[k]
                        
                        #dL/dx
                        if (k !=i && k !=j)
                            ddL = ddL*(xl - xk)/(xi - xk)
                        end
                    end
                    dLdx[i, l] = dLdx[i, l] + ddL/(xi - xj)
                end
            end
        end
    end

    return (L, dLdx)
end
