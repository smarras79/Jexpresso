#
# This file contains all the structs definitions
# S. Marras, Feb 2022
#
using  Quadmath
using  Polynomials
using  AMRVW 
export St_Lagrange
export St_Legendre
export St_lgl
export St_Laguerre
export St_gr
export build_lgl!
export LegendreGaussLobattoNodesAndWeights
export build_gr!

export St_lg
export St_cg


abstract type AbstractIntegrationPointAndWeights end
abstract type AbstractInterpolationBasis end
abstract type AbstractSpaceDimensions end

mutable struct St_Legendre{TFloat}
    
    """
    struct St_Legendre{TFloat<:Real}
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

mutable struct St_Laguerre{Float128}
   
   """
   struct St_Laguerre{TFloat<:Real}
        Laguerre  :: TFloat
        dLaguerre :: TFloat
   end
   """
   Laguerre  :: Polynomial{Float128}
   dLaguerre :: Polynomial{Float128}
   d2Laguerre :: Polynomial{Float128}
   d3Laguerre :: Polynomial{Float128}
end
 
abstract type AbstractIntegrationPointAndWeights end
abstract type AbstractInterpolationBasis end
abstract type AbstractSpaceDimensions end

struct LagrangeBasis <: AbstractInterpolationBasis end
struct ScaledLaguerreBasis <: AbstractInterpolationBasis end

mutable struct St_Chebyshev{TFloat} <:AbstractIntegrationPointAndWeights
    chebyshev ::TFloat
end

Base.@kwdef mutable struct St_lgl{TFloat,backend} <:AbstractIntegrationPointAndWeights
    ξ = KernelAbstractions.zeros(backend, TFloat, 0)
    ω = KernelAbstractions.zeros(backend, TFloat, 0)
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

Base.@kwdef mutable struct St_gr{TFloat, backend} <:AbstractIntegrationPointAndWeights
    ξ = KernelAbstractions.zeros(backend, TFloat, 0)
    ω = KernelAbstractions.zeros(backend, TFloat, 0)
end


Base.@kwdef mutable struct St_Lagrange{TFloat, backend} <:AbstractInterpolationBasis
    ψ  = KernelAbstractions.zeros(backend, TFloat, 0, 0)#mesh.ngl, mesh.ngl)
    dψ = KernelAbstractions.zeros(backend, TFloat, 0, 0)#mesh.ngl, mesh.ngl)    
end

mutable struct St_ScaledLaguerre{TFloat} <:AbstractInterpolationBasis
    ψ::Matrix{TFloat}
    dψ::Matrix{TFloat}
end

function basis_structs_ξ_ω!(::AbstractPointsType, nop, backend) end

function basis_structs_ξ_ω!(ξωtype::LG, nop, backend)
    
    lg = St_lg{TFloat}(KernelAbstractions.zeros(backend, TFloat, nop+1),
                       KernelAbstractions.zeros(backend, TFloat, nop+1))
    
    build_Integration_points!(lg, nop, backend)

    return lg
end

function basis_structs_ξ_ω!(ξωtype::LGL, nop, backend)
    #
    # Note: `nop` means `nq` when building the quadrature points in sem_setup.jl
    #
    lgl = St_lgl{TFloat, backend}(KernelAbstractions.zeros(backend, TFloat, nop+1),
                         KernelAbstractions.zeros(backend, TFloat, nop+1))
    
    build_Integration_points!(lgl, nop, backend)

    return lgl
end

function basis_structs_ξ_ω!(ξωtype::LGR, nop, beta, backend)

    lgr = St_gr{TFloat, backend}(KernelAbstractions.zeros(backend, TFloat, nop+1),
                         KernelAbstractions.zeros(backend, TFloat, nop+1))

    build_Integration_points!(lgr, nop, beta, backend)

    return lgr
end

function basis_structs_ξ_ω!(ξωtype::CG, nop, backend)
    
    #
    # Note: `nop` means `nq` when building the quadrature points in sem_setup.jl
    #
    cg = St_cg{TFloat}(KernelAbstractions.zeros(backend, TFloat, nop+1),
                       KernelAbstractions.zeros(backend, TFloat, nop+1))
    
    build_Integration_points!(cg, nop, backend)

    return cg
end


function basis_structs_ξ_ω!(ξωtype::CGL, nop, backend)
    
    cgl = St_cgl{TFloat}(KernelAbstractions.zeros(backend, TFloat, nop+1),
                         KernelAbstractions.zeros(backend, TFloat, nop+1))
    
    build_Integration_points!(cgl, nop, backend)

    return cgl
end


function build_Integration_points!(::AbstractIntegrationPointAndWeights,nop) end
 
function build_Integration_points!(cg::St_cg, nop, backend)
  build_cg!(cg,nop, backend)
end

function build_Integration_points!(cgl::St_cgl,nop, backend)
  build_cgl!(cgl,nop, backend)
end

function build_Integration_points!(lg::St_lg,nop, backend)
  Legendre = St_Legendre{TFloat}(0.0,0.0,0.0,0.0)
  build_lg!(Legendre,lg,nop, backend)
end

function build_Integration_points!(lgl::St_lgl,nop, backend)
  Legendre = St_Legendre{TFloat}(0.0,0.0,0.0,0.0)
  build_lgl!(Legendre,lgl,nop, backend)
end

function build_Integration_points!(lgr::St_gr,nop,beta, backend)
    Laguerre = St_Laguerre(Polynomial(Float128(2.0)),Polynomial(Float128(2.0)),Polynomial(Float128(2.0)),Polynomial(Float128(2.0)))
    build_gr!(Laguerre,lgr,nop,beta, backend)
end

function build_Interpolation_basis!(TP::LagrangeBasis, ξ, ξq, T, backend)

    Nξ = size(ξ,1)  - 1
    Qξ = size(ξq,1) - 1

    N  = (Nξ + 1)
    Q  = (Qξ + 1)
    
    basis = St_Lagrange{T, backend}(KernelAbstractions.zeros(backend, TFloat, N,Q), KernelAbstractions.zeros(backend, TFloat, N,Q))
    (basis.ψ, basis.dψ) = LagrangeInterpolatingPolynomials_classic(ξ, ξq, T, backend)
    
    return basis
end

function build_Interpolation_basis!(TP::ScaledLaguerreBasis, ξ, ξq, beta, T, backend)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0 @info "built laguerre basis ..... " end
    
    Nξ = size(ξ,1)  - 1
    Qξ = size(ξq,1) - 1

    N  = (Nξ + 1)
    Q  = (Qξ + 1)
    basis = St_Lagrange{T, backend}(KernelAbstractions.zeros(backend, TFloat, N,Q), KernelAbstractions.zeros(backend, TFloat, N,Q))
    (basis.ψ, basis.dψ) = LagrangeLaguerreBasis(ξ, ξq, beta,T, backend)
    if rank == 0 @info "built laguerre basis ..... DONE" end
    return basis
end


function build_lgl!(Legendre::St_Legendre, lgl::St_lgl, nop, backend)

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
    size::Int64 = nop+1
    
    lgl.ξ = KernelAbstractions.zeros(backend, TFloat, size)
    lgl.ω = KernelAbstractions.zeros(backend, TFloat, size)
    
    #LGL nodes
    LegendreGaussLobattoNodesAndWeights!(Legendre, lgl, nop, backend);
    
    return lgl;
end


function build_cg!(cg::St_cg, nop, backend)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0 println(" # Compute Chebyshev-Gauss nodes ........................ ") end
        
    size::Int8=nop+1
    cg.ξ = KernelAbstractions.zeros(backend, TFloat, size)
    cg.ω = KernelAbstractions.zeros(backend, TFloat, size)

    #CG nodes
    ChebyshevGaussNodesAndWeights!(cg,nop, backend)
    
    for j=1:size
        println( " # ξ cheby, ω =: ", " ", cg.ξ[j], " " , cg.ω[j])
    end
    
    if rank == 0 println(" # Compute Chebyshev-Gauss nodes ........................ END") end

    return cg
end

function build_cgl!(cgl::St_cgl, nop, backend)
    
    if rank == 0 println(" # Compute Chebyshev-Gauss-Lobatto nodes ........................ ") end
    
    size::Int8=nop+1
    cgl.ξ = KernelAbstractions.zeros(backend, TFloat, size)
    cgl.ω = KernelAbstractions.zeros(backend, TFloat, size)

    #CGL nodes
    ChebyshevGaussLobattoNodesAndWeights!(cgl,nop, backend)
    
    for j=1:size        
        println( " # ξ cheby, ω =: ", " ", cgl.ξ[j], " " , cgl.ω[j])
    end
    
    if rank == 0 println(" # Compute Chebyshev-Gauss-Lobatto nodes ........................ DONE") end
    
    return cgl
end

function build_lg!(Legendre::St_Legendre,lg::St_lg,nop, backend)
    size::Int8=nop+1
    lg.ξ = KernelAbstractions.zeros(backend, TFloat, size)
    lg.ω = KernelAbstractions.zeros(backend, TFloat, size)

    #LG nodes
    LegendreGaussNodesAndWeights!(Legendre,lg,nop, backend)
    return lg
end

function build_gr!(Laguerre::St_Laguerre,gr::St_gr,nop,beta, backend)
    size::Int8=nop+1
    gr.ξ = KernelAbstractions.zeros(backend, TFloat, Int64(size))
    gr.ω = KernelAbstractions.zeros(backend, TFloat, Int64(size))

    #LG nodes
    GaussRadauLaguerreNodesAndWeights!(Laguerre,gr,nop,beta, backend)
    return gr
end
function ChebyshevGaussNodesAndWeights!(cg::St_cg, nop, backend)
    """
         Compute the Nodes for the Chebyshev-Gauss Quadrature
         using Algorithm 26 of Kopriva's book 
      """
    for j=0:nop
        cg.ξ[j+1]=-cospi((2*j+1)/(2*nop+2))
        cg.ω[j+1]=π/(nop+1)
    end
end

function ChebyshevGaussLobattoNodesAndWeights!(cgl::St_cgl,nop, backend)
    """
         Compute the Nodes for the Chebyshev-Gauss-Lobatto Quadrature
         using Algorithm 26 of Kopriva's book 
      """
    for j=0:nop
        cgl.ξ[j+1]=-cospi(j/nop)
        cgl.ω[j+1]=π/nop
    end
    cgl.ω[1]=cgl.ω[1]/2
    cgl.ω[nop+1]=cgl.ω[nop+1]/2
end


function LegendreGaussNodesAndWeights!(Legendre::St_Legendre, lg::St_lg, nop, backend)
    """
          Compute the Nodes and Weights for the Legendre-Gauss Quadrature
          using Algorithm 23 of Kopriva's book valid for nop ≤  200
    """

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    
    println_rank( " # Compute LG nodes ........................"; msg_rank = rank)
    NITER = 100
    TOL = 4*eps(TFloat)
    Δ::TFloat=0.0
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
        lg.ξ[TInt(nop/2)+1] = 0
        lg.ω[TInt(nop/2)+1] = 2/Legendre.dlegendre^2
    end

    
    for j=1:nop+1       
        println_rank( " # ξ, ω =: ", " ", lg.ξ[j], " " , lg.ω[j]; msg_rank = rank)
    end
    
    println_rank(" # Compute LG nodes ........................ DONE"; msg_rank = rank)
    
end

function LegendreGaussLobattoNodesAndWeights!(Legendre::St_Legendre, lgl::St_lgl, nop, backend)
     """
          Compute the Nodes and Weights for the Legendre-Gauss-Lobatto Quadrature
     """
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    
    NITER = 100
    TOL = 4*eps()
    
    ξ0 ::TFloat=0.0
    ξ1 ::TFloat=1.0
    ξP ::TFloat=1.0
    ξj ::TFloat=0.0
    ξj2::TFloat=0.0
    L2 ::TFloat=0.0
    
    ω0 ::TFloat=0.0
    ω1 ::TFloat=0.0
    ωP ::TFloat=0.0
    
    Δ  ::TFloat=0.0
    ξ = zeros(TFloat,nop+1)
    ω = zeros(TFloat,nop+1)
    println_rank( " # Compute LGL nodes ........................"; msg_rank = rank)
    
    for j=1:nop+1
	ξ[j] = 0.0;
	ω[j] = 1.0;
    end
    
    if (nop == 1)
	ξ0           = -1.0
	ω0           =  1.0
	ξ1           =  1.0
	ω1           =   ω0
        
	ξ[1]     =   ξ0
	ξ[nop+1] =   ξ1
	ω[1]     =   ω0
	ω[nop+1] =   ω1
    else 
	ξ0            = -1.0
	ω0            =  TFloat(2.0/(nop*(nop + 1)))
	ξP            =  1.0
	ωP            =   ω0
        
	ξ[1]      =   ξ0
	ξ[nop+1]  =   ξP
	ω[1]      =   ω0
	ω[nop+1]  =   ωP
	
        for jj = 2:floor(Int,(nop + 1)/2) 
	    j = jj - 1
	    ξj = -cos((j + 0.25)*π/nop - 3.0/(8.0*nop*π*(j + 0.25)))
	    ξ[jj] = ξj;
            
            for k = 0:NITER
	        LegendreAndDerivativeAndQ!(Legendre, nop, ξj)
	        Δ = -Legendre.q/Legendre.dq
	        ξj    =  ξj + Δ
	        
	        if (abs(Δ) <= TOL*abs(ξj))
                    break
                end
	    end
            LegendreAndDerivativeAndQ!(Legendre, nop, ξj)
	    ξ[jj]      =  ξj
	    ξ[nop+1-j] = -ξj
	    xj2            =  ξj*ξj
	    L2             = Legendre.legendre*Legendre.legendre
	    ω[jj]      = 2.0/(nop*(nop + 1.0)*L2)
	    ω[nop+1-j] = ω[jj]
            
        end
    end
    
    if (mod(nop,2) == 0)
	LegendreAndDerivativeAndQ!(Legendre, nop, 0.0);
	ξ[TInt(nop/2)+1] = 0.0;
	
	L2           = Legendre.legendre*Legendre.legendre;
	ω[TInt(nop/2)+1] = 2.0/(nop*(nop + 1.0)*L2);
    end

    if (backend == CPU())
        lgl.ξ .= ξ
        lgl.ω .= ω
    else
        KernelAbstractions.copyto!(backend,lgl.ξ,ξ)
        KernelAbstractions.copyto!(backend,lgl.ω,ω)
    end
    for j=1:nop+1       
        println_rank( " # ξ, ω =: ", " ", ξ[j], " " , ω[j]; msg_rank = rank)
    end
    
    println_rank(" # Compute LGL nodes ........................ DONE"; msg_rank = rank)
    
end

function LegendreAndDerivativeAndQ!(Legendre::St_Legendre, nop, x)
    
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
    #TFloat=TFloat
    
    a   ::TFloat=0.0
    b   ::TFloat=0.0
    ϕ   ::TFloat=0.0
    dϕ  ::TFloat=0.0
    ϕp1 ::TFloat=0.0
    dϕp1::TFloat=0.0
    ϕm1 ::TFloat=0.0
    dϕm1::TFloat=0.0
    ϕm2 ::TFloat=0.0
    dϕm2::TFloat=0.0
    
    
    #st_legendre Legendre;  
    if (nop == 0) #Order 0 case
        ϕ    = 1.0
        dϕ   = 0.0
        q    = x
        dq   = 1.0
    elseif (nop == 1)
        ϕ    = x
        dϕ   = 1.0
        q    = 0.5*(3*x^2 - 2) - 1
        dq   = 3*x
    else
        ϕm2  = 1.0
	ϕm1  = x
	dϕm2 = 0.0
	dϕm1 = 1.0

        #Construct Nth Order Legendre Polynomial
	for k=2:nop
             
            ϕ    = x*ϕm1*(2.0*k - 1.0)/k - ϕm2*(k - 1.0)/k
	    dϕ   = dϕm2 + (2.0*k - 1.0)*ϕm1
	    
	    ϕp1  = x*ϕ*(2.0*k + 1.0)/(k+1) - ϕm1*k/(k+1)
	    dϕp1 = dϕm1 + (2.0*k - 1.0)*ϕ
            
            q  =  ϕp1 -  ϕm1
            dq = dϕp1 - dϕm1
            
	    ϕm2 = ϕm1
	    ϕm1 = ϕ

	    dϕm2 = dϕm1
	    dϕm1 = dϕ
            
	end
    end
    
    Legendre.legendre  =  ϕ
    Legendre.dlegendre = dϕ
    Legendre.q         =  q
    Legendre.dq        = dq
    
end

function ChebyshevPolynomial!(Chebyshev::St_Chebyshev,nop,x::TFloat,Ks)
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

function LagrangeInterpolatingPolynomials_classic(ξ, ξq, TFloat, backend)
"""
    LagrangeInterpolatingPolynomials_classic(ξ, ξq, N, Q, TFloat)
    ξ::set of N interpolation points (e.g. LGL points)
    ξq::point to interpolate to (e.g. quadrature points of points within the element)
    
    Algorithm 3.1 + 3.2 from Giraldo's book

    from https://github.com/fxgiraldo/Element-based-Galerkin-Methods/blob/master/Projects/Project_01_1D_Interpolation/For_Instructors/julia/lagrange_basis.jl

"""
    N = size(ξ,1) - 1
    Q = size(ξq,1) - 1
    
    #Initialize arrays
    L_1 = zeros(TFloat, N+1, Q+1)
    dLdx_1 = zeros(TFloat, N+1, Q+1)
    L    = KernelAbstractions.zeros(backend, TFloat, N+1, Q+1)
    dLdx = KernelAbstractions.zeros(backend, TFloat, N+1, Q+1)
    ξq_1 = zeros(TFloat,Q+1)
    ξ_1 = zeros(TFloat,N+1)
    KernelAbstractions.copyto!(CPU(),ξq_1,ξq)
    KernelAbstractions.copyto!(CPU(),ξ_1,ξ)

    for l=1:Q+1
        xl = ξq_1[l]

        #Construct Basis
        for i=1:N+1
            
            xi        = ξ_1[i]
            L_1[i,l]    = 1.0
            dLdx_1[i,l] = 0.0
            for j=1:N+1
                xj = ξ_1[j]

                #L
                if (j != i)
                    L_1[i,l] = L_1[i,l]*(xl - xj)/(xi - xj)
                end
                
                ddL=1
                if (j != i)
                    for k=1:N+1
                        xk = ξ_1[k]
                        
                        #dL/dx
                        if (k !=i && k !=j)
                            ddL = ddL*(xl - xk)/(xi - xk)
                        end
                    end
                    dLdx_1[i, l] = dLdx_1[i, l] + ddL/(xi - xj)
                end
            end
        end
    end
    if (backend == CPU())
        L .= L_1
        dLdx .= dLdx_1
    else
        KernelAbstractions.copyto!(backend,L,L_1)
        KernelAbstractions.copyto!(backend,dLdx,dLdx_1)
    end

    return (L, dLdx)
end

function ScaledLaguerreAndDerivative!(nop,SL::St_Laguerre,beta, backend)
  Laguerre = zeros(Float128,nop+1)
  if (nop == 0)
     Laguerre[1] = 1.0
  elseif (nop == 1)
     Laguerre[1] = 1.0
     Laguerre[2] = -1.0
  else
    Lkm2 = zeros(Float128,nop+1,1);
    Lkm2[nop+1] = 1;
    Lkm1 = zeros(Float128,nop+1,1);
    Lkm1[nop] = -1;
    Lkm1[nop+1] = 1;

    for k=2:nop

        Laguerre = zeros(Float128,nop+1);

        for e=nop-k+1:nop
            Laguerre[e] = (2*k-1)*Lkm1[e] - beta*Lkm1[e+1] + (1-k)*Lkm2[e];
        end

        Laguerre[nop+1] = (2*k-1)*Lkm1[nop+1] + (1-k)*Lkm2[nop+1];
        Laguerre = Laguerre/k;

        Lkm2 .= Lkm1;
        Lkm1 .= Laguerre;
    end
    for k=1:nop+1
      Laguerre[k] = Lkm1[nop+2-k]
    end
  end
  SL.Laguerre = Polynomial(Laguerre)
  SL.dLaguerre = Polynomials.derivative(SL.Laguerre)
  SL.d2Laguerre = Polynomials.derivative(SL.dLaguerre)
  SL.d3Laguerre = Polynomials.derivative(SL.d2Laguerre)
end

function LaguerreAndDerivative!(nop,SL::St_Laguerre, backend)
  Laguerre = zeros(Float128, nop+1)
  if (nop == 0)
     Laguerre[1] = 1.0
  elseif (nop == 1)
     Laguerre[1] = 1.0
     Laguerre[2] = -1.0
  else
    Lkm2 = zeros(nop+1,1);
    Lkm2[nop+1] = 1;
    Lkm1 = zeros(nop+1,1);
    Lkm1[nop] = -1;
    Lkm1[nop+1] = 1;

    for k=2:nop
        
        Laguerre = zeros(Float128, nop+1);

        for e=nop-k+1:nop
            Laguerre[e] = (2*k-1)*Lkm1[e] - Lkm1[e+1] + (1-k)*Lkm2[e];
        end
        
        Laguerre[nop+1] = (2*k-1)*Lkm1[nop+1] + (1-k)*Lkm2[nop+1];
        Laguerre = Laguerre/k;
        
        Lkm2 .= Lkm1;
        Lkm1 .= Laguerre;
    end
    for k=1:nop+1
      Laguerre[k] = Lkm1[nop+2-k]
    end
  end
  SL.Laguerre = Polynomial(Laguerre)
  SL.dLaguerre = Polynomials.derivative(SL.Laguerre)
end

function GaussRadauLaguerreNodesAndWeights!(Laguerre::St_Laguerre, gr::St_gr, nop,beta, backend)
    Pp1 = nop+1
    n = zeros(Float128,nop+1)
    bn = zeros(Float128,nop)
    an = zeros(Float128,nop+1)
    filler = zeros(Float128,nop+1)
    for i = 0:nop
       n[i+1] = i
       an[i+1] = (2 * n[i+1] + 1)/beta
    end
    for i = 1:nop
       bn[i] = i/beta
    end
    an[nop+1] = nop/beta
    J = zeros(Float128,nop+1,nop+1)
    J .= diagm(an) .+ Bidiagonal(filler,bn,:U) .+ Bidiagonal(filler,bn,:L)
    xi = eigen(J)
    ξ = Float128.(xi.values)
    ngr = length(gr.ξ)
    thresh = 1e-13
    if (backend != CPU())
    #if (backend == MetalBackend())
        thresh = 1e-5
    end
        
    x0 = 0.0
    x1 = 0.0
    for k=1:ngr
      x0 = ξ[k]
      diff1 = 1.0
      stuck = 1.0
      while(diff1 > thresh)
          ScaledLaguerreAndDerivative!(nop+1,Laguerre,beta,backend)
          L1 = Laguerre.Laguerre
          L3 = Laguerre.dLaguerre
          L4 = Laguerre.d2Laguerre
          L5 = Laguerre.d3Laguerre
          ScaledLaguerreAndDerivative!(nop,Laguerre,beta,backend)
          L2 = Laguerre.Laguerre
          #x1 = x0 - L1(x0)/L3(x0)#+ (L1(x0) - L2(x0))/L2(x0)
          #x1 = x0 -  L3(x0)/L4(x0) 
          # Laguerre root finder
          x1 = x0 - L3(x0)/L4(x0) - L3(x0)*L5(x0)/(2*(L4(x0)^2)) 
          diff1 = abs(x1 -x0)
          x0 = x1
          repeat = 1.0
          #=for k1=1:k
            if AlmostEqual(x1,gr.ξ[k1])
              repeat = 0.0
            end
          end
          if (repeat == 0.0)
             stuck +=1
             x0 = gr.ξ[k]*stuck
             diff = 1.0
          end=#
          
      end
      ξ[k] = x1
    end
    ξ[1] = 0
    #ScaledLaguerreAndDerivative!(nop+1,Laguerre,beta)
    #gr.ξ[2:ngr] = AMRVW.roots((coeffs(Laguerre.dLaguerre)))
    #for k=1:ngr
    #  @info Laguerre.dLaguerre(gr.ξ[k]), gr.ξ[k]
    #end
    ω = zeros(Float128,ngr) 
    ScaledLaguerreAndDerivative!(nop,Laguerre,beta,backend)
    Lkx = zeros(nop+1,1)
    for i=1:nop+1
      Lkx[i] = scaled_laguerre(ξ[i],nop,beta,backend)
      #Lkx[i] = Laguerre.Laguerre(gr.ξ[i])
      ω[i] = 1/(beta*Pp1*Lkx[i]^2)
      #gr.ω[i] = exp(gr.ξ[i]*beta)/(beta*Pp1*Lkx[i]^2)
    
    end
    if (backend == CPU())
        gr.ξ .= ξ
        gr.ω .= ω
    else
        T = eltype(gr.ξ)
        KernelAbstractions.copyto!(backend, gr.ξ, T.(ξ))
        KernelAbstractions.copyto!(backend, gr.ω, T.(ω))
    end
    #@info gr.ξ
    #@info gr.ω
    #gr.ω[1] = 1-sum(gr.ω[2:nop+1])
    #@info gr.ω
    #if(scale)
      #gr.ω .= exp.(gr.ξ*beta).*gr.ω 
    #end
end

function LagrangeLaguerreBasis(ξ, ξq, beta, TFloat, backend)
    nbasis = size(ξq,1)
    N = nbasis -1
    Np1 = N+1
    
    psi_1 = ones(TFloat, nbasis, nbasis)
    dpsi_1 = zeros(TFloat, nbasis, nbasis)
    ξq_1 = zeros(TFloat, nbasis)
    ξ_1 = zeros(TFloat, nbasis)
    psi = KernelAbstractions.ones(backend,TFloat, nbasis,nbasis)
    dpsi = KernelAbstractions.zeros(backend, TFloat, nbasis,nbasis)
    KernelAbstractions.copyto!(CPU(), ξq_1, ξq)
    KernelAbstractions.copyto!(CPU(), ξ_1, ξ)
    
    dpsi_1[1,1]=-beta*((N +1)./2.0)
    for i = 1:nbasis
        xi = ξq_1[i]
        for j = 1:nbasis
            xj = ξ_1[j]
            if(i != j)
                psi_1[i,j] = 0.0
                dpsi_1[j,i] = scaled_laguerre(xi,Np1,beta,backend)/(scaled_laguerre(xj,Np1,beta,backend)*(xi -xj));
            end
        end
    end
    if (backend == CPU())
        psi .= psi_1
        dpsi .= dpsi_1
    else
        KernelAbstractions.copyto!(backend,psi,psi_1)
        KernelAbstractions.copyto!(backend,dpsi,dpsi_1)
    end

    return (psi,dpsi)


end


function scaled_laguerre(x,n,beta,backend)
    Laguerre = St_Laguerre(Polynomial(TFloat(2.0)),Polynomial(TFloat(2.0)),Polynomial(TFloat(2.0)),Polynomial(TFloat(2.0)))
    ScaledLaguerreAndDerivative!(n,Laguerre,beta,CPU())
    #Lkx = Laguerre.Laguerre(x)
    Lkx = Real(Laguerre.Laguerre(x))
    y = exp(-(beta*x)/2)*Lkx#exp(-x)*Lkx
    return y
end
