TInt=Int64
TFloat=Float64
include("Kopriva_functions.jl")
include("TimeIntegrators.jl")
rk3 = RK3_Integrator{TFloat}(zeros(TFloat,3),zeros(TFloat,3),zeros(TFloat,3))
buildRK3Integrator!(rk3)
N=16
Nout = 100
xer = zeros(Float64,Nout+1)
Φr = zeros(Complex, Nout+1)
x= zeros(Float64, N+1)
Φ=zeros(Float64, N+1)
ν=0.2
T=1.0
for i=0:N
  x[i+1] = 2*π*i/N
  Φ[i+1] = 3/(5-4*cos(x[i+1]))
end
for i=0:Nout
  xer[i+1]=2*π*i/Nout
  for k = -256:256
    Φr[i+1] = Φr[i+1] + 2.0^(-abs(k))*exp(im*k*(xer[i+1]-T)-ν*k^2*T)
  end
end
NT=200
Φ = FourierCollocationDriver(N,NT,T,Φ,ν,rk3.a,rk3.b,rk3.g)
Φint = zeros(Float64,Nout+1)
w=BarycentricWeights(x)
for i=0:Nout
  Φint[i+1] = LagrangeInterpolation(xer[i+1],x,Φ,w)
end
er = 0.0
for i=0:Nout
  global er = er + (2*pi/Nout)*(Φr[i+1]-Φint[i+1])
end
@info log(er)
