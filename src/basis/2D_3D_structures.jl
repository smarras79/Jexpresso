using LinearAlgebra
include("basis_structs.jl")
include("Kopriva_functions.jl")
abstract type NodalStorage end
abstract type Abstract_Integration_Points end
mutable struct Nodal2DStorage <:NodalStorage
  N::Int64
  M::Int64
  ξ::AbstractIntegrationPointAndWeights
  η::AbstractIntegrationPointAndWeights
  Dξ::Array{Float64}
  Dη::Array{Float64}
  D2ξ::Array{Float64}
  D2η::Array{Float64}
end


mutable struct Nodal3DStorage <: NodalStorage
  N::Int64
  M::Int64
  L::Int64
  ξ::AbstractIntegrationPointAndWeights
  η::AbstractIntegrationPointAndWeights
  ζ::AbstractIntegrationPointAndWeights
  Dξ::Array{Float64}
  Dη::Array{Float64}
  Dζ::Array{Float64}
  D2ξ::Array{Float64}
  D2η::Array{Float64}
  D2ζ::Array{Float64}
end


struct LGL2D <: Abstract_Integration_Points end
struct CGL2D <: Abstract_Integration_Points end
struct LGL3D <: Abstract_Integration_Points end
struct CGL3D <: Abstract_Integration_Points end
abstract type FDPreconditioner end

mutable struct NodalPotential
  PT::Abstract_Integration_Points
  dim::Int64
  dims::Array{Int64}
  ND::NodalStorage
  Φ::Array{Float64}
  s::Array{Float64}
  mask::Array{Int64}
end


function build_nodal_2DStorage_cgl(N,M)
  cgl      = St_cgl{TFloat}(zeros(TFloat, N+1),
                       zeros(TFloat, N+1))
  build_Integration_points!(cgl,N)
  ξ = cgl
  Dξ = PolynomialDerivativeMatrix(cgl.ξ)
  D2ξ = mthOrderPolynomialDerivativeMatrix(2,cgl.ξ)
  
  cgl      = St_cgl{TFloat}(zeros(TFloat, M+1),
                       zeros(TFloat, M+1))
  build_Integration_points!(cgl,M)
  η = cgl
  Dη = PolynomialDerivativeMatrix(cgl.ξ)
  D2η = mthOrderPolynomialDerivativeMatrix(2,cgl.ξ)
  ND = Nodal2DStorage(N,M,ξ,η,Dξ,Dη,D2ξ,D2η)
  return ND
end

function build_nodal_2DStorage_lgl(N,M)
   lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                      zeros(TFloat, N+1))
   build_Integration_points!(lgl,N)
   ξ = lgl
   Dξ = PolynomialDerivativeMatrix(lgl.ξ)
   D2ξ = mthOrderPolynomialDerivativeMatrix(2,lgl.ξ)
   lgl      = St_lgl{TFloat}(zeros(TFloat, M+1),
                       zeros(TFloat, M+1))
   build_Integration_points!(lgl,M)
   η = lgl
   Dη = PolynomialDerivativeMatrix(lgl.ξ)
   D2η = mthOrderPolynomialDerivativeMatrix(2,lgl.ξ)
   ND = Nodal2DStorage(N,M,ξ,η,Dξ,Dη,D2ξ,D2η)
   return ND
end

function build_nodal_3DStorage_cgl(N,M,L)
  cgl      = St_cgl{TFloat}(zeros(TFloat, N+1),
                       zeros(TFloat, N+1))
  build_Integration_points!(cgl,N)
  ξ=cgl
  Dξ=PolynomialDerivativeMatrix(cgl,ξ)
  D2ξ=mthOrderPolynomialDerivativeMatrix(2,cgl.ξ)
  cgl      = St_cgl{TFloat}(zeros(TFloat, M+1),
                            zeros(TFloat, M+1))
  build_Integration_points!(cgl,M)
  η=cgl
  Dη=PolynomialDerivativeMatrix(cgl.ξ)
  D2η=mthOrderPolynomialDerivativeMatrix(2,cgl.ξ)
  cgl      = St_cgl{TFloat}(zeros(TFloat, L+1),
                            zeros(TFloat, L+1))
  ζ=cgl
  Dζ=PolynomialDerivativeMatrix(cgl.ξ)
  D2ζ=mthOrderPolynomialDerivativeMatrix(2,cgl.ξ)
  ND=Nodal3DStorage(N,M,L,ξ,η,ζ,Dξ,Dη,Dζ,D2ξ,D2η,D2ζ)
  return ND
end

function build_nodal_3DStorage_lgl(N,M,L)
  lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                            zeros(TFloat, N+1))
  build_Integration_points!(lgl,N)
  ξ=lgl
  Dξ=PolynomialDerivativeMatrix(lgl.ξ)
  D2ξ=mthOrderPolynomialDerivativeMatrix(2,lgl.ξ)
  lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                            zeros(TFloat, N+1))
  build_Integration_points!(lgl,N)
  η=lgl
  Dη=PolynomialDerivativeMatrix(lgl.ξ)
  D2η=mthOrderPolynomialDerivativeMatrix(2,lgl.ξ)
  lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                            zeros(TFloat, N+1))
  build_Integration_points!(lgl,N)
  ζ=lgl
  Dζ=PolynomialDerivativeMatrix(lgl.ξ)
  D2ζ=mthOrderPolynomialDerivativeMatrix(2,lgl.ξ)
  ND=Nodal3DStorage(N,M,L,ξ,η,ζ,Dξ,Dη,Dζ,D2ξ,D2η,D2ζ)
  return ND
end 


function build_nodal_Storage(dims,PT::LGL2D)
  ND = build_nodal_2DStorage_lgl(dims[1],dims[2])
  return ND
end

function build_nodal_Storage(dims,PT::CGL2D)
  ND = build_nodal_2DStorage_cgl(dims[1],dims[2])
  return ND
end

function build_nodal_Storage(dims,PT::LGL3D)
  ND=  build_nodal_3DStorage_lgl(dims[1],dims[2],dims[3])
  return ND
end

function build_nodal_Storage(dims,PT::CGL3D)
  ND = build_nodal_3DStorage_cgl(dims[1],dims[2],dims[3])
  return ND
end

function build_Nodal_Potential(dim,dims,PT::Abstract_Integration_Points)
  if (dim == 2)
    Φ=zeros(Float64,dims[1]+1,dims[2]+1)
    s=zeros(Float64,dims[1]+1,dims[2]+1)
    mask=zeros(Float64,4)
  else
    Φ=zeros(Float64,dims[1]+1,dims[2]+1,dims[3]+1)
    mask =zeros(Float64,8)
    s=zeros(Float64,dims[1]+1,dims[2]+1,dims[3]+1)
  end
  ND = build_nodal_Storage(dims,PT)
  NP = NodalPotential(PT,dim,dims,ND,Φ,s,mask)
  return NP
end

function  ComputeLaplacian(U,NP::NodalPotential)
  dim=NP.dim
  N=NP.ND.N
  M=NP.ND.M
  if (dim ==2)
    Uxx=zeros(Float64,N+1,M+1)
    for j=1:M+1
      Uxx[:,j] = MxVDerivative(NP.ND.D2ξ,U[:,j])
    end
    Uyy = zeros(Float64,N+1,M+1)
    for i=1:N+1
      Uyy[i,:] = MxVDerivative(NP.ND.D2η,U[i,:])
    end
    @info Uxx
    @info Uyy
    LapU = zeros(Float64,N+1,M+1)
    for j=1:M+1
      for i=1:N+1
        LapU[i,j]=Uxx[i,j]+Uyy[i,j]
      end
    end

  else
    L=NP.ND.L
    Uxx=zeros(Float64,N+1,M+1,L+1)
    Uyy=zeros(Float64,N+1,M+1,L+1)
    Uzz=zeros(Float64,N+1,M+1,L+1)
    for j=1:M+1
      for l=1:L+1
        Uxx[:,j,l] = MxVDerivative(NP.ND.D2ξ,U[:,j,l])
      end
    end
    for i=1:N+1
      for l=1:L+1
        Uyy[i,:,l] = MxVDerivative(NP.ND.D2η,U[i,:,l])
      end
    end
    for i=1:N+1
      for j=1:M+1
        Uzz[i,j,:] = MxVDerivative(NP.ND.D2ζ,U[i,j,:])
      end
    end
    LapU=zeros(Float64,N+1,M+1,L+1)
    LapU .= Uxx .+ Uyy .+ Uzz
  end
  return LapU
end

function MaskSidesFaces(U,NP::NodalPotential)
  mask=NP.mask
  dim=NP.dim
  N=NP.ND.N
  M=NP.ND.M
  if (dim ==2)
    if (mask[1]==1)
      U[:,1] .=0.0
    end
    if (mask[2]==1)
      U[N+1,:] .=0.0
    end
    if (mask[3]==1)
      U[:,M+1] .=0.0
    end
    if (mask[4]==1)
      U[1,:] .= 0.0
    end
  else
    L=NP.ND.L
    if (mask[1]==1)
      U[:,1,1] .= 0.0
    end
    if (mask[2]==1)
      U[N+1,:,L+1] .=0.0
    end
    if (mask[3]==1)
      U[:,M+1,L+1] .=0.0
    end
    if (mask[4]==1)
      U[1,:,1] .=0.0
    end
    if (mask[5]==1)
      U[1,1,:] .= 0.0
    end
    if (mask[6]==1)
      U[N+1,M+1,:] .=0.0
    end
  end
  return U
end

function MatrixAction(U,NP::NodalPotential)
  Action = ComputeLaplacian(U,NP)
  Action = MaskSidesFaces(U,NP)
  return Action
end

function CollocationRHSComputation(NP::NodalPotential)
  N=NP.ND.N
  M=NP.ND.M
  dim=NP.dim
  if (dim==2)
    L=N*M
    rhs=zeros(Float64,L+1)
    for j=1:M
      for i=1:N
        n=i+(j-1)*M
        rhs[n] = NP.s[i,j] - NP.ND.D2ξ[i,1] * NP.Φ[1,j] - NP.ND.D2ξ[i,N+1] * NP.Φ[N+1,j] - NP.ND.D2η[j,1] * NP.Φ[i,1] - NP.ND.D2η[j,M+1] * NP.Φ[i,M+1]
      end
    end
    return rhs
  else
    Ll=NP.ND.L
    L = N*M*Ll
    rhs=zeros(Float64,L+1)
    for l=1:Ll
      for j=1:M
        for i=1:N
          n=i+M*(j-1+Ll*(l-1))
          rhs[n] = NP.s[i,j] - NP.ND.D2ξ[i,1] * NP.Φ[1,j,1] - NP.ND.D2ξ[i,N+1] * NP.Φ[N+1,j,Ll+1] - NP.ND.D2η[j,1] * NP.Φ[i,1,1] - NP.ND.D2η[j,M+1] * NP.Φ[i,M+1,L+1],
          - NP.ND.D2ζ[l,1]*NP.Φ[1,1,l] - NP.ND.D2ζ[l,Ll+1] * NP.Φ[1,1,Ll+1]
        end
      end
    end
    return rhs
  end
end

function LaplaceCollocationMatrix(NP::NodalPotential)
  N=NP.ND.N
  M=NP.ND.M
  dim=NP.dim
  if (dim==2)
    L=N*M
    A=zeros(Float64,L+1,L+1)
    for j=1:M
      for i=1:N
        n=i+(j-1)*M
        for k=1:N
          m=k+(j-1)*M
          A[n,m] = NP.ND.D2ξ[i,k]
        end
        for k=1:M
          m=i+(k-1)*M
          A[n,m] = A[n,m] + NP.ND.D2η[j,k]
        end
      end
    end
  else
    Ll=NP.ND.L
    L=N*M*Ll
    A=zeros(Float64,L+1,L+1)
    for l=1:Ll
      for j=1:M
        for i=1:N
          n=i+M*(j-1+Ll*(l-1))
          for k=1:N
            m=k+M*(j-1+Ll*(l-1))
            A[n,m] = NP.ND.D2ξ[i,k]
          end
          for k=1:M
            n=i+M*(k-1+Ll*(l-1))
            A[n,m] = A[n,m] +NP.ND.D2η[j,k]
          end
          for k=1:Ll
            n=i+M*(j-1+Ll*(k-1))
            A[n,m] = A[n,m] + NP.ND.D2ζ[l,k]
          end
        end
      end
    end
  end
  return A
end

function Residual(NP::NodalPotential)
  N=NP.ND.N
  M=NP.ND.M
  dim=NP.dim
  if (dim==2)
    L=(N+1)*(M+1)
    r=ComputeLaplacian(NP.Φ,NP)
    r=LinearAlgebra.BLAS.scal(L,-1.0,r,1)
    LinearAlgebra.axpy!(1,NP.s,r)
    r=MaskSidesFaces(r,NP)
  else
    Ll=NP.ND.L
    L=(N+1)*(M+1)*(Ll+1)
    r=ComputeLaplacian(NP.Φ,NP)
    r=BLAS.scal(L,-1.0,r,1)
    LinearAlgebra.axpby!(1,NP.s,1,r)
    r=MaskSidesFaces(r,NP)
  end
  return r
end

mutable struct FDPreconditioner2D <:FDPreconditioner
  N::Int64
  M::Int64
  a::Array{Float64}
  Δx::Array{Float64}
  Δy::Array{Float64}
  A::Array{Float64}
  B::Array{Float64}
  C::Array{Float64}
  E::Array{Float64}
  F::Array{Float64}
end

mutable struct FDPreconditioner3D <:FDPreconditioner
  N::Int64
  M::Int64
  L::Int64
  a::Array{Float64}
  Δx::Array{Float64}
  Δy::Array{Float64}
  Δz::Array{Float64}
end

function buildPreconditioner!(FDP::FDPreconditioner) end

function initializeFDPrecondtioner2D(N,M)
  a=zeros(Float64,N-1,M-1)
  Δx=zeros(Float64,N)
  Δy=zeros(Float64,M)
  A=zeros(Float64,N-1,M-1)
  B=zeros(Float64,N-1,M-1)
  C=zeros(Float64,N-1,M-1)
  E=zeros(Float64,N-1,M-1)
  F=zeros(Float64,N-1,M-1)
  FDP = FDPreconditioner2D(N,M,a,Δx,Δy,A,B,C,E,F)
  return FDP
end

function buildPreconditioner!(FDP::FDPreconditioner2D,x,y)
  N=FDP.N
  M=FDP.M
  for i=1:N
    FDP.Δx[i]=x[i+1]-x[i]
  end
  for j=1:M
    FDP.Δy[j] = y[j+1]-y[j]
  end
  for i=1:N-1
    for j=1:M-1
      FDP.A[i,j]=-2*(1/(FDP.Δx[i]*FDP.Δx[i+1])+1/(FDP.Δy[j]*FDP.Δy[j+1]))
      FDP.B[i,j]=2/(FDP.Δx[i]*(FDP.Δx[i]+FDP.Δx[i+1]))
      FDP.C[i,j]=2/(FDP.Δy[j]*(FDP.Δy[j]+FDP.Δy[j+1]))
      FDP.E[i,j]=2/(FDP.Δx[i+1]*(FDP.Δx[i]+FDP.Δx[i+1]))
      FDP.F[i,j]=2/(FDP.Δy[j+1]*(FDP.Δy[j]+FDP.Δy[j+1]))
    end
  end
  FDP.a[1,1]=FDP.A[1,1]
  for i=2:N-1
    FDP.a[i,1] = FDP.A[i,1] - FDP.B[i,1]*FDP.E[i-1,1]/FDP.a[i-1,1] - FDP.B[i,1]*FDP.F[i-1,1]/FDP.a[i-1,1]
  end
  for j=2:M-1
    FDP.a[1,j]=FDP.A[1,j] - FDP.C[1,j]*FDP.F[1,j-1]/FDP.a[1,j-1] - FDP.C[1,j]*FDP.E[1,j-1]/FDP.a[1,j-1]
    for i=2:N-1
      FDP.a[i,j]=FDP.A[i,j] - FDP.B[i,j]*FDP.E[i-1,j]/FDP.a[i-1,j]-FDP.C[i,j]*FDP.F[i,j-1]/FDP.a[i,j-1]-FDP.B[i,j]*FDP.F[i-1,j]/FDP.a[i-1,j]-FDP.C[i,j]*FDP.E[i,j-1]/FDP.a[i,j-1]
    end
  end
end

#BUILD this for 3D once working on square

function Solve(FDP::FDPreconditioner2D,R)
  N=FDP.N
  M=FDP.M
  w=zeros(Float64,N-1,M-1)
  w[1,1] = R[2,2]/FDP.a[1,1]
  u=zeros(Float64,N+1,M+1)
  for i=2:N-1
    w[i,1] = (R[i,1] -FDP.B[i,1]*w[i-1,1])/FDP.a[i,1]
  end
  for j=2:M-1
    w[1,j] = (R[1,j]-FDP.C[1,j]*w[1,j-1])/FDP.a[1,j]
    for i=2:N-1
      w[i,j]=(R[i,j] -FDP.B[i,j]*w[i-1,j] - FDP.C[i,j]*w[i,j-1])/FDP.a[i,j]
    end
  end
  u[N-1,M-1] = w[N-1,M-1]
  for i=N-2:-1:1
    u[i,M-1] = w[i,M-1] - FDP.E[i,M-1]*u[i+1,M-1]/FDP.a[i,M-1]
  end
  for j=M-2:-1:1
    u[N-1,j]=w[N-1,j]-FDP.F[N-1,j]*u[N-1,j+1]/FDP.a[N-1,j]
    for i=N-2:1:-1
      u[i,j] = w[i,j] - FDP.E[i,j]*u[i+1,j]/FDP.a[i,j] - FDP.F[i,j]*u[i,j+1]/FDP.a[i,j]
    end
  end
  return u
end

#TODO build this for 3D once working 

function BiCGStabSolve(Nit::Int64,Tol,NP::NodalPotential,FDP::FDPreconditioner2D)
  N=NP.ND.N
  M=NP.ND.M
  L=(N+1)*(M+1)
  ρ=1.0
  α=1.0
  ω=1.0
  β=1.0
  r= Residual(NP)
  r̄=zeros(Float64,N+1,M+1)
  r̄=LinearAlgebra.BLAS.blascopy!(L,r,1,r̄,1)
  v=zeros(Float64,N+1,M+1)
  s=zeros(Float64,N+1,M+1)
  p=zeros(Float64,N+1,M+1)
  ρ=BLAS.dot(L,r̄,1,r,1)
  ρ̂=ρ
  for k=1:Nit
    @info BLAS.nrm2(L,r,1)
    #β=ρ*α/(ρ̂*ω)
    p=LinearAlgebra.axpy!(-ω,v,p)
    p=LinearAlgebra.BLAS.scal(L,β,p,1)
    p=LinearAlgebra.axpy!(1,r,p)
    y=Solve(FDP,p)
    v=MatrixAction(y,NP)
    α=ρ/BLAS.dot(L,r̄,1,v,1)
    LinearAlgebra.BLAS.blascopy!(L,r,1,s,1)
    LinearAlgebra.axpy!(-α,v,s)
    z=Solve(FDP,s)
    t=MatrixAction(z,NP)
    ω=BLAS.dot(L,t,1,s,1)/BLAS.dot(L,t,1,t,1)
    NP.Φ=LinearAlgebra.axpy!(α,y,NP.Φ)
    NP.Φ=LinearAlgebra.axpy!(ω,z,NP.Φ)
    r=LinearAlgebra.BLAS.blascopy!(L,s,1,r,1)
    r=LinearAlgebra.axpy!(-ω,t,r)
    ρ̂=ρ
    ρ=BLAS.dot(L,r̄,1,r,1)
    β=ρ*α/(ρ̂*ω)
    if (BLAS.nrm2(L,r,1)<Tol)
      return NP
    end
  end
  @info "did not converge", Nit,Tol,BLAS.nrm2(L,r,1)
end

#TODO see if the same would apply fine to 3D

function CollocationPotentialDriver2D!(NP,FDP,Nit,Tol)
  N=NP.ND.N
  M=NP.ND.M
  NP=BiCGStabSolve(Nit,Tol,NP,FDP)
end
