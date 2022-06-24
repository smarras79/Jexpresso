TInt=Int64
TFloat=Float64
include("../src/Infrastructure/2D_3D_structures.jl")
N=40
M=40
Nit=10
Tol=0.001
P1 = LGL2D()
P2 = CGL2D()
T1 = Collocation()
T2 = Nodal_Galerkin()
dim=2
dims=[N M]
NP = build_Nodal_Potential(dim,dims,P1,T2)
for i=2:N
  for j=2:M
    NP.s[i,j] = -8*π^2*cospi(2*NP.ND.ξ.ξ[i])*sinpi(2*NP.ND.η.ξ[j])
    NP.Φ[i,j] = cospi(2*NP.ND.ξ.ξ[i])*sinpi(2*NP.ND.η.ξ[j])
  end
end

NP.Φ[:,1].= 0.0
NP.Φ[N+1,:].= sinpi.(2 .*NP.ND.η.ξ[:])
NP.Φ[:,M+1].= 0.0
NP.Φ[1,:].= sinpi.(2 .*NP.ND.η.ξ[:])
NP.mask = [1 1 1 1]
#rhs = CollocationRHSComputation(NP)
#y=LUSolve(A1,p,rhs)
rhs = NodalGalerkinRHS(NP)
A = NodalGalerkinMatrix(NP)
A2,p1,info = LinearAlgebra.LAPACK.getrf!(A)
LinearAlgebra.LAPACK.getrs!('N',A2,p1,rhs)
@info rhs
@info NP.Φ
rhs = NodalGalerkinRHS(NP)
A = NodalGalerkinMatrix(NP)
@info rhs
@info A\rhs
