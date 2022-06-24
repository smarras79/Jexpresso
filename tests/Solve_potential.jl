TInt=Int64
TFloat=Float64
include("../src/Infrastructure/2D_3D_structures.jl")
N=40
M=40
Nit=100
Tol=0.1
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
NP.Φ[:,1] .= 0.0 #NP.ND.ξ.ξ[:].^3 .- 1
NP.Φ[N+1,:].= sinpi.(2 .*NP.ND.η.ξ[:])
NP.Φ[:,M+1].= 0.0 #NP.ND.ξ.ξ[:].^3 .+ 1
NP.Φ[1,:].= sinpi.(2 .*NP.ND.η.ξ[:])
@info NP.Φ
NP.mask = [1 1 1 1]
FDP = initializeFDPrecondtioner2D(N,M)
buildPreconditioner!(FDP,NP.ND.ξ.ξ,NP.ND.η.ξ)
#CollocationPotentialDriver2D!(NP,FDP,Nit,Tol)
NP = PreconditionedConjugateGradientSolve(Nit,Tol,NP);
@info NP.Φ

