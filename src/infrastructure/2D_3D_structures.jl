using LinearAlgebra
include("../basis/basis_structs.jl")
include("Kopriva_functions.jl")
abstract type NodalStorage end
abstract type Abstract_Integration_Points end
abstract type Abstract_Method_Type end
abstract type Nodal_Solver end


mutable struct Nodal1DStorage <:NodalStorage
    N::Int64
    ξ::AbstractIntegrationPointAndWeights
    Dξ::Array{Float64}
    D2ξ::Array{Float64}
end

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


struct LGL_1D <: Abstract_Integration_Points end
struct CGL_1D <: Abstract_Integration_Points end
struct LGL_2D <: Abstract_Integration_Points end
struct CGL_2D <: Abstract_Integration_Points end
struct LGL_3D <: Abstract_Integration_Points end
struct CGL_3D <: Abstract_Integration_Points end
struct Collocation <: Abstract_Method_Type end
struct NodalGalerkin <: Abstract_Method_Type end

abstract type FDPreconditioner end

mutable struct NodalPotential <: Nodal_Solver
    PT::Abstract_Integration_Points
    dim::Int64
    dims::Array{Int64}
    ND::NodalStorage
    Φ::Array{Float64}
    s::Array{Float64}
    mask::Array{Int64}
    T::Abstract_Method_Type
end

mutable struct NodalAdvDiff <: Nodal_Solver
    u::Float64
    v::Float64
    ν::Float64
    PT::Abstract_Integration_Points
    dim::Int64
    dims::Array{Int64}
    ND::NodalStorage
    Φ::Array{Float64}
    transport::Array{Float64}
    RHS::Array{Float64}
    mask::Array{Float64}
    p::Array{Int64}
    T::Abstract_Method_Type
end

#1D
function build_nodal_1DStorage_cgl(N,T::Collocation)
    cgl      = St_cgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    build_Integration_points!(cgl,N)
    ξ = cgl
    Dξ = PolynomialDerivativeMatrix(cgl.ξ)
    D2ξ = mthOrderPolynomialDerivativeMatrix(2,cgl.ξ)
    
    ND = Nodal1DStorage(N,ξ,Dξ,D2ξ,ψ,dψ)
    return ND
end

function build_nodal_1DStorage_cgl(N,T::NodalGalerkin)
    cgl      = St_cgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    build_Integration_points!(cgl,N)
    ξ = cgl
    Dξ = PolynomialDerivativeMatrix(cgl.ξ)
    D2ξ = CGDerivativeMatrix(N)
    
    ND = Nodal1DStorage(N,ξ,Dξ,D2ξ)
    return ND
end

function build_nodal_1DStorage_lgl(N,T::Collocation)
    lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    build_Integration_points!(lgl,N)
    ξ = lgl
    Dξ = PolynomialDerivativeMatrix(lgl.ξ)
    D2ξ = mthOrderPolynomialDerivativeMatrix(2,lgl.ξ)
    
    build_Integration_points!(lgl,M)
    ND = Nodal1DStorage(N,ξ,Dξ,D2ξ)
    
    return ND
end
function build_nodal_1DStorage_lgl(N,T::NodalGalerkin)
    lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    build_Integration_points!(lgl,N)
    ξ = lgl
    
    Dξ = PolynomialDerivativeMatrix(lgl.ξ)
    D2ξ = CGDerivativeMatrix(N)
    
    ND = Nodal1DStorage(N,ξ,Dξ,D2ξ)
    return ND
end


#2D
function build_nodal_2DStorage_cgl(N,M,T::Collocation)
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

function build_nodal_2DStorage_cgl(N,M,T::NodalGalerkin)
    cgl      = St_cgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    build_Integration_points!(cgl,N)
    ξ = cgl
    Dξ = PolynomialDerivativeMatrix(cgl.ξ)
    D2ξ = CGDerivativeMatrix(N)
    
    cgl      = St_cgl{TFloat}(zeros(TFloat, M+1),
                              zeros(TFloat, M+1))
    build_Integration_points!(M)
    η = cgl
    Dη = PolynomialDerivativeMatrix(cgl.ξ)
    D2η = CGDerivativeMatrix(M)
    ND = Nodal2DStorage(N,M,ξ,η,Dξ,Dη,D2ξ,D2η)
    return ND
end

function build_nodal_2DStorage_lgl(N,M,T::Collocation)
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

function build_nodal_2DStorage_lgl(N,M,T::NodalGalerkin)
    lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    build_Integration_points!(lgl,N)
    ξ = lgl
    Dξ = PolynomialDerivativeMatrix(lgl.ξ)
    D2ξ = CGDerivativeMatrix(N)
    lgl      = St_lgl{TFloat}(zeros(TFloat, M+1),
                              zeros(TFloat, M+1))
    build_Integration_points!(lgl,M)
    η = lgl
    Dη = PolynomialDerivativeMatrix(lgl.ξ)
    D2η = CGDerivativeMatrix(M)
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

function build_nodal_3DStorage_lgl(N,M,L,T::Collocation)
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

function build_nodal_3DStorage_lgl(N,M,L,T::NodalGalerkin)
    lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    build_Integration_points!(lgl,N)
    ξ=lgl
    Dξ=PolynomialDerivativeMatrix(lgl.ξ)
    D2ξ=mthOrderPolynomialDerivativeMatrix(N)
    lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    build_Integration_points!(lgl,N)
    η=lgl
    Dη=PolynomialDerivativeMatrix(lgl.ξ)
    D2η=CGDerivativeMatrix(M)
    lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    build_Integration_points!(lgl,N)
    ζ=lgl
    Dζ=PolynomialDerivativeMatrix(lgl.ξ)
    D2ζ=CGDerivativeMatrix(L)
    ND=Nodal3DStorage(N,M,L,ξ,η,ζ,Dξ,Dη,Dζ,D2ξ,D2η,D2ζ)
    return ND
end 


#1D LGL
function build_nodal_Storage(dims,PT::LGL_1D,T)
    ND = build_nodal_1DStorage_lgl(dims[1],T)
    return ND
end

#1D CGL
function build_nodal_Storage(dims,PT::CGL_1D,T)
    ND = build_nodal_1DStorage_cgl(dims[1],T)
    return ND
end

#2DLGL
function build_nodal_Storage(dims,PT::LGL_2D,T)
    ND = build_nodal_2DStorage_lgl(dims[1],dims[2],T)
    return ND
end

#2D CGL
function build_nodal_Storage(dims,PT::CGL_2D,T)
    ND = build_nodal_2DStorage_cgl(dims[1],dims[2],T)
    return ND
end

#3D LGL
function build_nodal_Storage(dims,PT::LGL_3D,T)
    ND=  build_nodal_3DStorage_lgl(dims[1],dims[2],dims[3],T)
    return ND
end

#3D CGL
function build_nodal_Storage(dims,PT::CGL_3D,T)
    ND = build_nodal_3DStorage_cgl(dims[1],dims[2],dims[3])
    return ND
end

function build_Nodal_Potential(dim,dims,PT::Abstract_Integration_Points,T::Abstract_Method_Type)
    if (dim == 2)
        Φ=zeros(Float64,dims[1]+1,dims[2]+1)
        s=zeros(Float64,dims[1]+1,dims[2]+1)
        mask=zeros(Float64,4)
    else
        Φ=zeros(Float64,dims[1]+1,dims[2]+1,dims[3]+1)
        mask =zeros(Float64,8)
        s=zeros(Float64,dims[1]+1,dims[2]+1,dims[3]+1)
    end
    ND = build_nodal_Storage(dims,PT,T)
    NP = NodalPotential(PT,dim,dims,ND,Φ,s,mask,T)
    return NP
end

function build_Nodal_AdvDiff(u,v,ν,dim,dims,PT::Abstract_Integration_Points,T::Abstract_Method_Type)
    if (dim == 2)
        Φ=zeros(Float64,dims[1]+1,dims[2]+1,3)
        transport=zeros(Float64,dims[1]+1,dims[2]+1,3)
        RHS = zeros(Float64,dims[1]+1,dims[2]+1)
        mask = zeros(Float64,4)
        p = [1 2 3]
    else
        Φ = zeros(Float64, dims[1]+1, dims[2]+1, dims[3]+1,3)
        transport = zeros(Float64, dims[1]+1, dims[2]+1, dims[3]+1,3)
        RHS = zeros(Float64, dims[1]+1, dims[2]+1, dims[3]+1)
        mask = zeros(Float64,8)
        p = [1 2 3]
    end
    ND = build_nodal_Storage(dims,PT,T)
    NAD = NodalAdvDiff(u,v,ν,PT,dim,dims,ND,Φ,transport,RHS,mask,p,T)
    return NAD
end

function transport!(NAD::NodalAdvDiff,T::Collocation,k)
    N=NAD.ND.N
    M=NAD.ND.M
    Φx = zeros(Float64,N+1,N+1)
    Φy = zeros(Float64,M+1,M+1)
    for j=1:M+1
        Φx[:,j] = MxVDerivative(NAD.ND.Dξ,NAD.Φ[:,j])
    end
    for i=1:N+1
        Φy[i,:] = MxVDerivative(NAD.ND.Dη,NAD.Φ[i,:])
    end
    for j=1:M+1
        for i=1:N+1
            NAD.transport[i,j,k] = NAD.u * Φx[i,j] + NAD.v * Φy[i,j]
        end
    end
    NAD.transport = MaskSidesFaces(NAD.transport, NAD)
end

function transport!(NAD::NodalAdvDiff,T::NodalGalerkin,k)
    N=NAD.ND.N
    M=NAD.ND.M
    Φx = zeros(Float64,N+1,N+1)
    Φy = zeros(Float64,M+1,M+1)
    for j=1:M+1
        Φx[:,j] = MxVDerivative(NAD.ND.Dξ,NAD.Φ[:,j])
    end
    for i=1:N+1
        Φy[i,:] = MxVDerivative(NAD.ND.Dη,NAD.Φ[i,:])
    end
    for j=1:M+1
        for i=1:N+1
            NAD.transport[i,j,k] = NAD.u * Φx[i,j] + NAD.v * Φy[i,j]
            NAD.transport[i,j,k] = NAD.ND.ξ.ω[i] * NAD.ND.η.ω[j] * NAD.transport[i,j,k]
        end
    end
    NAD.transport = MaskSidesFaces(NAD.transport, NAD)
end

function AdvDiffExplicitRHS!(NAD::NodalAdvDiff,T::Collocation,Δt)
    n=NAD.p[1]
    nm1=NAD.p[2]
    nm2=NAD.p[3]
    M=NAD.ND.M
    N=NAD.ND.N
    for j=1:M+1
        for i =1:N+1
            NAD.RHS[i,j] = (1/11) * (18*NAD.Φ[i,j,n] - 9* NAD.Φ[i,j,nm1] + 2*NAD.Φ[i,j,nm2])
            NAD.RHS[i,j] = NAD.RHS[i,j] - (6*Δt/11) * (3*NAD.transport[i,j,n]-3*transport[i,j,nm1] + NAD.transport[i,j,nm2])
        end
    end
    NAD.RHS = MaskSidesFaces(NAD.RHS,NAD)
end

function AdvDiffExplicitRHS!(NAD::NodalAdvDiff,T::NodalGalerkin,Δt)
    n=NAD.p[1]
    nm1=NAD.p[2]
    nm2=NAD.p[3]
    M=NAD.ND.M
    N=NAD.ND.N
    for j=1:M+1
        for i =1:N+1
            NAD.RHS[i,j] = (1/11) * (18*NAD.Φ[i,j,n] - 9* NAD.Φ[i,j,nm1] + 2*NAD.Φ[i,j,nm2])
            NAD.RHS[i,j] = NAD.ND.ξ.ω[i] * NAD.ND.η.ω[j] * NAD.RHS[i,j]
            NAD.RHS[i,j] = NAD.RHS[i,j] - (6*Δt/11) * (3*NAD.transport[i,j,n]-3*transport[i,j,nm1] + NAD.transport[i,j,nm2])
        end
    end
    NAD.RHS = MaskSidesFaces(NAD.RHS,NAD)
end


function  ComputeLaplacian(U,NP::Nodal_Solver,::Collocation)
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

function MaskSidesFaces(U,NP::Nodal_Solver)
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
    Action = ComputeLaplacian(U,NP,NP.T)
    Action = MaskSidesFaces(Action,NP)
    return Action
end

function CollocationRHSComputation(NP::NodalPotential)
    N=NP.ND.N
    M=NP.ND.M
    dim=NP.dim
    if (dim==2)
        L=(N-1)*(M-1)
        rhs=zeros(Float64,L)
        for j=1:M-1
            for i=1:N-1
                n=j+(i-1)*(N-1)
                rhs[n] = NP.s[i+1,j+1] - NP.ND.D2ξ[i+1,1] * NP.Φ[1,j+1] - NP.ND.D2ξ[i+1,N+1] * NP.Φ[N+1,j+1] - NP.ND.D2η[1,j+1] * NP.Φ[i+1,1] - NP.ND.D2η[M+1,j+1] * NP.Φ[i+1,M+1]
            end
        end
        return rhs
    else
        Ll=NP.ND.L
        L = N*M*Ll
        rhs=zeros(Float64,L)
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

function NodalGalerkinRHS(NP::NodalPotential)
    N=NP.ND.N
    M=NP.ND.M
    dim=NP.dim
    L=(N-1)*(M-1)
    wx=NP.ND.ξ.ω
    wy=NP.ND.η.ω
    Gx=NP.ND.D2ξ
    Gy=NP.ND.D2η
    Φ=NP.Φ
    rhs=zeros(Float64,L)
    for j=1:M-1
        for i=1:N-1
            n=j+(i-1)*(N-1)
            rhs[n]=NP.s[i+1,j+1]*wx[i+1]*wy[j+1] - wy[j+1]*Gx[i+1,1]*Φ[1,j+1]-wy[j+1]*Gx[i+1,N+1]*Φ[N+1,j+1] - wx[i+1]*Gy[j+1,1]*Φ[i+1,1]-wx[i+1]*Gy[j+1,M+1]*Φ[i+1,M+1]
        end
    end
    return rhs
end

function LaplaceCollocationMatrix(NP::NodalPotential)
    N=NP.ND.N
    M=NP.ND.M
    dim=NP.dim
    if (dim==2)
        L=(N-1)*(M-1)
        A=zeros(Float64,L,L)
        for j=1:M-1
            for i=1:N-1
                n=j+(i-1)*(M-1)
                for k=1:N-1
                    if (k != i)
                        m=j+(k-1)*(M-1)
                        A[n,m] = NP.ND.D2ξ[i+1,k+1]
                    end
                end
                for k=1:M-1
                    if (k!= j)
                        m=k+(i-1)*(M-1)
                        A[n,m] = NP.ND.D2η[j+1,k+1]
                    end
                end
                A[n,n] = NP.ND.D2ξ[i+1,i+1] + NP.ND.D2η[j+1,j+1]
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

function NodalGalerkinMatrix(NP::NodalPotential)
    N=NP.ND.N
    M=NP.ND.M
    Gx=NP.ND.D2ξ
    Gy=NP.ND.D2η
    wx=NP.ND.ξ.ω
    wy=NP.ND.η.ω
    L=(N-1)*(M-1)
    A=zeros(Float64,L,L)
    for j=1:M-1
        for i=1:N-1
            n=j+(i-1)*(M-1)
            for k=1:N-1
                if (k != i)
                    m=j+(k-1)*(M-1)
                    A[n,m] = wy[j+1]*Gx[i+1,k+1]
                end
            end
            for k=1:M-1
                if (k !=j)
                    m=k+(i-1)*(M-1)
                    A[n,m] = wx[i+1]*Gy[j+1,k+1]
                end
            end
            A[n,n] = wy[j+1]*Gx[i+1,i+1] + wx[i+1]*Gy[j+1,j+1]
        end
    end
    return A
end

function Residual(NP::NodalPotential,::Collocation)
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
        w[i,1] = (R[i+1,1] -FDP.B[i,1]*w[i-1,1])/FDP.a[i,1]
    end
    for j=2:M-1
        w[1,j] = (R[2,j+1]-FDP.C[1,j]*w[1,j-1])/FDP.a[1,j]
        for i=2:N-1
            w[i,j]=(R[i+1,j+1] -FDP.B[i,j]*w[i-1,j] - FDP.C[i,j]*w[i,j-1])/FDP.a[i,j]
        end
    end
    u[N,M] = w[N-1,M-1]
    for i=N-1:-1:2
        u[i,M] = w[i-1,M-1] - FDP.E[i-1,M-1]*u[i+1,M]/FDP.a[i-1,M-1]
    end
    for j=M-1:-1:2
        u[N,j]=w[N-1,j-1]-FDP.F[N-1,j-1]*u[N,j+1]/FDP.a[N-1,j-1]
        for i=N-1:-1:2
            u[i,j] = w[i-1,j-1] - FDP.E[i-1,j-1]*u[i+1,j]/FDP.a[i-1,j-1] - FDP.F[i-1,j-1]*u[i,j+1]/FDP.a[i-1,j-1]
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
    r= Residual(NP)
    r̄=zeros(Float64,N+1,M+1)
    r̄=LinearAlgebra.BLAS.blascopy!(L,r,1,r̄,1)
    v=zeros(Float64,N+1,M+1)
    s=zeros(Float64,N+1,M+1)
    p=zeros(Float64,N+1,M+1)
    for k=1:Nit
        @info "norm",BLAS.nrm2(L,r,1)
        ρ̂=ρ
        ρ=BLAS.dot(L,r̄,1,r,1)
        β=ρ*α/(ρ̂*ω)
        r̄=LinearAlgebra.BLAS.blascopy!(L,r,1,r̄,1)
        @info "ρ, ρ̂,β",ρ, ρ̂,β
        p=LinearAlgebra.axpy!(-ω,v,p)
        @info "p=-ωv+p",p
        p=LinearAlgebra.BLAS.scal(L,β,p,1)
        @info "p=β(-ωv+p)",p
        p=LinearAlgebra.axpy!(1,r,p)
        @info "p=r+β(-ωv+p)",p
        y=Solve(FDP,p)
        @info "y, Ay=p", Solve(FDP,p)
        v=MatrixAction(y,NP)
        @info "v=Lap(y)",v
        α=ρ/BLAS.dot(L,v,1,r̄,1)
        @info "α",α
        LinearAlgebra.BLAS.blascopy!(L,r,1,s,1)
        @info "s=r",r
        @info "s=r",s
        LinearAlgebra.axpy!(-α,v,s)
        @info "s=-αv+s",s
        @info "norm s", BLAS.nrm2(L,s,1)
        z=Solve(FDP,s)
        @info "z, Az=s",z
        t=MatrixAction(z,NP)
        @info "t=Lap(z)",t
        ω=BLAS.dot(L,t,1,s,1)/BLAS.dot(L,t,1,t,1)
        @info "ω",ω
        NP.Φ=LinearAlgebra.axpy!(α,y,NP.Φ)
        @info "Φ=Φ+αp",NP.Φ
        NP.Φ=LinearAlgebra.axpy!(ω,z,NP.Φ)
        @info "Φ=Φ+αp+ωz", NP.Φ
        r=LinearAlgebra.BLAS.blascopy!(L,s,1,r,1)
        @info "r=s",r
        r=LinearAlgebra.axpy!(-ω,t,r)
        @info "r=r -ωt",r
        @info "Res", Residual(NP)
        if (BLAS.nrm2(L,r,1)<Tol)
            return NP
        end
    end
    @info "did not converge", Nit,Tol,BLAS.nrm2(L,r,1)
end

function CGSolve(Nit,Tol,NP,FDP)
    N=NP.ND.N
    M=NP.ND.M
    L=(N+1)*(M+1)
    ρ=1.0
    α=1.0
    ω=1.0
    r= Residual(NP)
    p=zeros(Float64,N+1,M+1)
    p=LinearAlgebra.BLAS.blascopy!(L,r,1,p,1)
    for j=1:Nit
        α=BLAS.dot(L,r,1,r,1)/BLAS.dot(L,MatrixAction(Solve(FDP,p),NP),1,p,1)
        LinearAlgebra.axpy!(α,p,NP.Φ)
        r1=LinearAlgebra.axpy!(-α,MatrixAction(Solve(FDP,p),NP),r)
        if (BLAS.nrm2(L,r,1)<Tol)
            return NP
        end
        β=BLAS.dot(L,r1,1,r1,1)/BLAS.dot(L,r,1,r,1)
        LinearAlgebra.axpby!(1,r1,β,p)
        r=LinearAlgebra.BLAS.blascopy!(L,r1,1,r,1)
    end
end
#TODO see if the same would apply fine to 3D

function CollocationPotentialDriver2D!(NP,FDP,Nit,Tol)
    N=NP.ND.N
    M=NP.ND.M
    #NP=BiCGStabSolve(Nit,Tol,NP,FDP)
    NP=CGSolve(Nit,Tol,NP,FDP)
end

function FactorizeLU(A)
    N= size(A,1)
    p=zeros(Int64,N)
    for k=1:N
        p[k] = k
        for i=k+1:N
            if abs(A[i,k]) > abs(A[p[k],k]) 
                p[k] = i
            end
        end
        if (p[k] != k)
            for j=1:N
                t=A[k,j]
                A[k,j] = A[p[k],j]
                A[p[k],j] = t
            end
        end
    end
    for k=1:N
        for j=k:N
            s=0.0
            for n=1:N
                s=s+A[k,n]*A[n,j]
            end
            A[k,j]=A[k,j]-s
        end
        for i=k+1:N
            s=0.0
            for n=1:k-1
                s=s+A[i,n]*A[n,k]
            end
            A[i,k] = (A[i,k]-s)/A[k,k]
        end
    end
    return A,p
end

function LUSolve(A,p,rhs)
    N=size(p,1)
    w=zeros(Float64,N)
    y=zeros(Float64,N)
    y[:].=rhs[:]
    for i=1:N
        if (p[i] != i)
            for m=1:NRHS
                t = y[i,m]
                y[i,m] = y[p[i],m]
                y[p[i],m] = t
            end
        end
    end
    w[1] = y[1]
    for i=2:N
        s=0.0
        for j=1:i-1
            s=s+A[i,j]*w[j]
        end
        w[i]=y[i]-s
    end
    y[N] = w[N]/A[N,N]
    for i=N-1:-1:1
        s=0.0
        for j=i+1:N
            s=s+A[i,j]*y[j]
        end
        y[i] = (w[i]-s)/A[i,j]
    end
    return y
end

function ComputeLaplacian(U,NP::Nodal_Solver,T::NodalGalerkin)
    N=NP.ND.N
    M=NP.ND.M
    Uxx = zeros(Float64,N+1,M+1)
    Uyy = zeros(Float64,N+1,M+1)
    Lap = zeros(Float64,N+1,M+1)
    for j=1:M+1
        Uxx[:,j] = MxVDerivative(NP.ND.D2ξ,U[:,j])
        for i =1:N+1
            Uxx[i,j] = NP.ND.η.ω[j] * Uxx[i,j]
        end
    end
    for i=1:N+1
        Uyy[i,:] = MxVDerivative(NP.ND.D2η,U[i,:])
    end
    for j=1:M+1
        for i=1:N+1
            Lap[i,j] = -Uxx[i,j] - NP.ND.ξ.ω[i] * Uyy[i,j]
        end
    end
    return Lap
end

function Residual(NP::NodalPotential,T::NodalGalerkin)
    N=NP.ND.N
    M=NP.ND.M
    L=(N+1)*(M+1)
    r=ComputeLaplacian(NP.Φ,NP,NP.T)
    r=LinearAlgebra.BLAS.scal(L,-1.0,r,1)
    LinearAlgebra.axpy!(1,NP.s,r)
    r=MaskSidesFaces(r,NP)
end

function Ψ_Ξ(k,l,η)
    result = -(-1)^k * (1-l-(-1)^l * η)
    return result
end

function Ψ_Η(k,l,ξ)
    result = -(-1)^l * (1-k-(-1)^k * ξ)
    return result
end

function LocalStiffnessMatrix(Δx,Δy)
    Ŝ=zeros(Float64,4,4)
    for m=0:1
        for n=0:1
            q=n+2*m+1
            for l=0:1
                for k=0:1
                    p=k+2*l+1
                    t=0.0
                    for s=0:1
                        for r=0:1
                            t = t+Ψ_Ξ(k,l,s)*Ψ_Ξ(n,m,s)/Δx^2 + Ψ_Η(k,l,r)*Ψ_Η(n,m,r)/Δy^2 
                        end
                    end
                    Ŝ[p,q] = Δx * Δy * t / 4
                end
            end
        end
    end
    return Ŝ
end

function ApproximateFEMStencil(i,j,x,y)
    pŜ=zeros(Float64,4,4,4)
    C=zeros(Float64,3,3)
    for m=0:1
        for n=0:1
            p=n+m*2+1 
            Δx = x[i-n+2]-x[i-n+1]
            Δy = y[j-m+2]-y[j-m+1]
            pŜ[:,:,p] = LocalStiffnessMatrix(Δx,Δy)
        end
    end
    for m=0:1
        for n=0:1
            p=n+2*m+1
            for k=-n:-n+1
                for l=-m:-m+1
                    C[k+2,l+2] = C[k+2,l+2] + pŜ[k+n+1,l+m+1,p]
                end
            end
        end
    end
    return C
end

function SSORSweep(r,ω,NP::NodalPotential)
    N=NP.ND.N
    M=NP.ND.M
    z=zeros(Float64,N+1,M+1)
    for j=1:M-1
        for i=1:N-1
            s=0.0
            C=ApproximateFEMStencil(i,j,NP.ND.ξ.ξ,NP.ND.η.ξ)
            for k=-1:1
                for l=-1:1
                    s=s+C[k+2,l+2]*z[i+k+1,j+l+1]
                end
            end

            z[i+1,j+1] = z[i+1,j+1] + ω*(r[i+1,j+1]-s)/C[2,2]
        end
    end
    for j=M-1:-1:1
        for i=N-1:-1:1
            s= 0.0
            C=ApproximateFEMStencil(i,j,NP.ND.ξ.ξ,NP.ND.η.ξ)
            for k=-1:1
                for l=-1:1
                    s= s+C[k+2,l+2] * z[i+k+1,j+l+1]
                end
            end
            z[i+1,j+1] = z[i+1,j+1] + ω*(r[i+1,j+1]-s)/C[2,2]
        end
    end
    return z
end

function PreconditionedConjugateGradientSolve(Nit,Tol,NP::NodalPotential)
    N=NP.ND.N
    M=NP.ND.M
    L=(N+1)*(M+1)
    v=zeros(Float64,N+1,M+1)
    r = Residual(NP,NP.T)
    z=SSORSweep(r,1.43,NP)
    LinearAlgebra.BLAS.blascopy!(L,z,1,v,1)
    c = BLAS.dot(L,r,1,z,1)
    for k=1:Nit
        @info BLAS.nrm2(L,r,1)
        z=MatrixAction(v,NP)
        ω=c/BLAS.dot(L,v,1,z,1)
        LinearAlgebra.axpy!(ω,v,NP.Φ)
        LinearAlgebra.axpy!(-ω,z,r)
        if (BLAS.nrm2(L,r,1) <= Tol)
            @info "converged"
            return NP
        end
        z = SSORSweep(r,1+k/Nit,NP)
        d = BLAS.dot(L,r,1,z,1)
        v = LinearAlgebra.BLAS.scal(L,d/c,v,1)
        LinearAlgebra.axpy!(1,z,v)
        c = d
    end
    @info "did not converge"
    return NP
end

function MatrixAction(NAD::NodalAdvDiff,Δt,U,T::Collocation)
    Lap = ComputeLaplacian(U,NAD,T)
    M=NAD.ND.M
    N=NAD.ND.N
    action = zeros(Float64,N+1,M+1)
    for j=1:M+1
        for i=1:N+1
            action[i,j] = U[i,j] - (6*Δt/11) * Lap[i,j]
        end
    end
    action = MaskSidesFaces(Action,NAD)
    return Action
end

function MatrixAction(NAD::NodalAdvDiff,Δt,U,T::NodalGalerkin)
    Lap = ComputeLaplacian(U,NAD,T)
    M=NAD.ND.M
    N=NAD.ND.N
    action = zeros(Float64,N+1,M+1)
    for j=1:M+1
        for i=1:N+1
            action[i,j] = NAD.ND.ξ.ω[i] * NAD.ND.η.ω[j] * U[i,j] - (6*Δt/11) * Lap[i,j]
        end
    end
    action = MaskSidesFaces(Action,NAD)
    return Action
end

function Residual(NAD::NodalAdvDiff,Δt,U)
    action = MatrixAction(NAD,Δt,U,NAD.T)
    M=NAD.ND.M
    N=NAD.ND.N
    r = zeros(Float64,N+1,M+1)
    for j=1:M+1
        for i=1:N+1
            r[i,j] = NAD.RHS[i,j] - action[i,j]
        end
    end
    r[i,j] = MaskSidesFaces(r,NAD)
    return r
end

function MultistepIntegration(NAD::NodalAdvDiff, FD::FDPreconditioner, Δt, T::Collocation)
    NAD.RHS = AdvDiffExplicitRHS(NAD,T,Δt)
    tmp = NAD.p[3]
    NAD.p[3] = NAD.p[2]
    NAD.p[2] = NAD.p[1]
    NAD.p[1] = tmp
    NAD.Φ = SetBoundaryConditions(t+Δt, NAD)
    NAD = BiCGSSTABSolve(Nit,Tol,NAD, FD)
    NAD.transport[:,:,NAD.p[1]] = Transport(NAD,NAD.p[1])
    return NAD
end
