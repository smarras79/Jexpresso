include("basis_structs.jl")
include("Kopriva_functions.jl")
mutable struct Nodal2DStoragecgl
  N::Int64
  M::Int64
  ξ::St_cgl
  η::St_cgl
  Dξ::Array{Float64}
  Dη::Array{Float64}
  D2ξ::Array{Float64}
  D2η::Array{Float64}
end

function buildNodal2DStorage(N,M,N2D::Nodal2DStoragecgl)
  Chebyshev = St_chebyshev{TFloat}(0.0)
  cgl      = St_cgl{TFloat}(zeros(TFloat, N+1),
                       zeros(TFloat, N+1))
  ChebyshevGaussLobattoNodesAndWeights(Chebyshev,cgl,N)
  N2D.N = N
  N2D.M = M
  N2D.ξ = cgl
  N2D.Dξ = PolynomialDerivativeMatrix(cgl.ξ)
  N2D.D2ξ = mthOrderPolynomialDerivativeMatrix(2,cgl.ξ)
  Chebyshev = St_chebyshev{TFloat}(0.0)
  cgl      = St_cgl{TFloat}(zeros(TFloat, M+1),
                       zeros(TFloat, M+1))
  ChebyshevGaussLobattoNodesAndWeights(Chebyshev,cgl,M)
  N2D.η = cgl
  N2D.Dη = PolynomialDerivativeMatrix(cgl.ξ)
  N2D.D2η = mthOrderPolynomialDerivativeMatrix(2,cgl.ξ)
end

mutable struct Nodal2DStoragelgl
  N::Int64
  M::Int64
  ξ::St_lgl
  η::St_lgl
  Dξ::Array{Float64}
  Dη::Array{Float64}
  D2ξ::Array{Float64}
  D2η::Array{Float64}
 end
 
 function buildNodal2DStorage(N,M,N2D::Nodal2DStoragelgl)
   Legendre = St_legendre{TFloat}(0.0,0.0,0.0,0.0)
   lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                      zeros(TFloat, N+1))
   LegendreGaussLobattoNodesAndWeights(Legendre,lgl,N)
   N2D.N = N
   N2D.M = M
   N2D.ξ = lgl
   N2D.Dξ = PolynomialDerivativeMatrix(lgl.ξ)
   N2D.D2ξ = mthOrderPolynomialDerivativeMatrix(2,lgl.ξ)
   Legendre = St_legendre{TFloat}(0.0,0.0,0.0,0.0)
   lgl      = St_lgl{TFloat}(zeros(TFloat, M+1),
                       zeros(TFloat, M+1))
   LegendreGaussLobattoNodesAndWeights(Legendre,lgl,M)
   N2D.η = lgl
   N2D.Dη = PolynomialDerivativeMatrix(lgl.ξ)
   N2D.D2η = mthOrderPolynomialDerivativeMatrix(2,lgl.ξ)
end
