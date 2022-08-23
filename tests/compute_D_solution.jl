TInt=Int64
TFloat=Float64
include("../src/Infrastructure/Kopriva_functions.jl")
include("../src/basis/basis_structs.jl")
a=zeros(Float64,3)
b=zeros(Float64,3)
g=zeros(Float64,3)
a[1] = 0.0
b[1] = 0.0
g[1] = 1/3
a[2] = -5/9
b[2] = 1/3
g[2] = 15/16
a[3] = -153/128
b[3] = 3/4
g[3] = 8/15
Nout = 50
N_2 = floor(Int64,Nout/2)
T = 0.02
NT= 100.0
ν=0.2
xer = zeros(Float64,Nout+1)
Φr = zeros(Float64,Nout+1)
N=15
Φ = zeros(Float64,N+1)
Legendre = St_legendre{TFloat}(0.0, 0.0, 0.0, 0.0)
lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                       zeros(TFloat, N+1))
LegendreGaussLobattoNodesAndWeights!(Legendre,lgl,N)
x=lgl.ξ
w=lgl.ω
err = 0.0
for j=1:N+1
  Φ[j] = sinpi(x[j]+1)
end
DG = NodalDiscontinuousGalerkin(N,zeros(Float64,N+1,N+1),zeros(Float64,N+1),zeros(Float64,N+1),zeros(Float64,N+1),zeros(Float64,N+1))
buildNodalDiscontinuousGalerkin!(N,DG)
DG.Φ=Φ
ad=2
Φ = LegendreCollocationIntegrator(N,NT,Nout,T,Φ,a,b,g,0.0,0.0,ad)
#Φ = CGDriver(N,NT,Nout,T,Φ,a,b,g,0.0,0.0)
#Φ=DGDriver!(N,NT,Nout,T,DG,a,b,g,1.0,0.0)
for i=1:Nout+1
  xer[i] = -1+2*(i-1)/Nout
  Φr[i] = sinpi(xer[i]+1)*exp(-π^2*T)
  #Φr[i] = sinpi(xer[i]+1 - T)
  #Φr[i] = sinpi(xer[i]+1 - T)*exp(-π^2*T)
end
for i=1:Nout+1
  global err = err +(2*π/Nout)*(Φr[i]-Φ[i])^2
  # global err = err +(2*π/Nout)*(Φr[i]-DG.Φ[i])^2
end

@info err

