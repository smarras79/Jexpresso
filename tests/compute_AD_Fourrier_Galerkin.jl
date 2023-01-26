include("../src/kernel/infrastructure/Kopriva_functions.jl")
include("../src/kernel/basis/basis_structs.jl")
include("../src/io/plotting/jeplots.jl")

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
N_2  = floor(Int64,Nout/2)
T    = 2.0
NT   = 400.0
ν    = 0.2
xer  = zeros(Float64,Nout+1)
Φr   = zeros(Complex,Nout+1)
#Φ    = zeros(Complex,Nout+1)
er   = zeros(Complex,Nout+1)

for i=1:Nout+1
    xer[i] = 2*π*(i-1)/Nout
    for j=-128:128
        Φr[i] = Φr[i] + 2.0^(-abs(j))*exp(im*j*(xer[i]-T) - ν*j^2*T)
    end
end
for N=6:6 #Nout
    N2=floor(Int64,N/2)
    Φ̂ = zeros(Complex,N+1)
    for k=1:N+1
        iter = k-1-N2
        Φ̂[k] = 2.0^(-abs(iter))
    end
    Φ = FourierGalerkinDriver(Φ̂,N,NT,T,Nout,ν,a,b,g)

  #=  plot_error(range(1, Nout, length=length(er)), Φ,
           "Error",
           [""],
           yscale=:log10)
=#
  
    for i=1:Nout+1
        er[N-5] = er[N-5] +(2*pi/Nout)*(Φr[i]-Φ[i])^2
    end
end
#=
@info length(er)
@info er
plot_error(range(1, Nout, length=length(er)), Φ,
           "Error",
           [""],
           yscale=:log10)

=#
#@info er

