""" DiscreteFourrierCoefficients(f)
    f:: array of length N
    computes the Discrete fourrier coefficients of f
    f̃ using Algorithm 1 of Kopriva
"""
function DiscreteFourierCoefficients(f)
  N=Int32(size(f,1))
  for k=-(N/2):N/2
    s=0
    for j = 0:N-1
      s=s+f[j+1]*exp(-2*pi*Complex(0,1)*j*k/N)
    end 
    f̃[k+N/2+1] =s/N
  end
  return f̃
end
""" FourrierInterpolantFromModes(f̃,x)
     f̃:: array of length N of discrete fourrier coefficients
     computes the fourrier interpolant of f at the point x using the coefficients f̃
     using Algorithm 2 of Kopriva
 """
function FourierInterpolantFromModes(f̃,x)
  i=Complex(0,1)
  N=size(f̂)
  s=(f̃[1] * exp(-i * N * x / 2) + f̃[N] * exp(i * N * x * 2)) / 2
  for k = -(N/2):(N/2)-1
    s=s + f̃ * exp(i * k * x)
  end
  return Real(s)
end
"""
   AlmostEqual(a,b)
   uses algorithm 139 of Kopriva's book to
   determine if two floating point numbers a and b are nearly equal or not
"""
function AlmostEqual(a,b)
  ϵ=eps(typeof(a))
  if (a == 0) || (b == 0)
    if (abs(a-b) ≤ 2*ϵ)
      return true
    else
      return false
    end
  else
    if (abs(a-b)≤ϵ*abs(a)) && (abs(a-b)≤\epsilon*b)
      return true
    else
      return false
    end
  end
end
""" FourrierInterpolantFromNodes(f,x,x_j)
     f:: array of length N
     x::real
     computes the fourrier interpolant of f at the point x using the coefficients the value of f at the nodes x_j
     using Algorithm 3 of Kopriva
 """
function FourierInterpolantFromNodes(f,x,x_j)
  N=Int32(size(f,1))
  for j=0:N-1
    if AlmostEqual(x,x_j[j+1])
      return f[j+1]
    end
  end
  s=0
  for j=0:N-1
    t=(x-x_j[j+1])/2
    s=s+f[j+1]*sin(N*t)*cot(t)/N
  end
  Intf = s
  return Intf
end
"""
  LegendreDerivativeCoefficients(f̂)
  f̂:: coefficients due to polynomial truncation of f
  determines the coefficients of the derivative of a truncated Legendre polynomial series
  using Algorithm 4 of Kopriva's book
"""
function LegendreDerivativeCoefficients(f̂)
  N=Int32(size(f,1))-1
  f̅_1[N+1]=0.0
  f̅_1[N] = (2*N-1)*f̂[N+1]
  for k=N-2:0:-1
    f̅_1[k+1]=(2*k+1)*(f̂[k+2]+f̅_1[k+3]/(2*k+5))
  end
  return f̅_1
end
"""
   ChebyshevDerivativeCoefficients(f̂)
   f̂:: coefficients to due to polynomial truncation of f
   determines the coefficients of the derivative of a truncated Chebyshev polynomial series
   using Algorithm 5 of Kopriva's book
"""
function ChebyshevDerivativeCoefficients(f̂)
  N=Int32(size(f,1))-1
  f̅_1[N+1]=0.0
  f̅_1[N] = 2*N*f̂[N+1]
  for k=N-2:0:-1
    f̅_1[k+1] = 2*(k+1)*f̂[k+2]+f̅_1[k+3]
  end
  f̅_1=f̂[2]+f̅_1[3]/2
  return f̅_1
end
"""
   DFT(f,s)
   f::array of size N
   s::Int32
   computes the discrete forward or backwards fourrier transform of f
   depending on the sign of s using Algorithm 6 of Kopriva's book
"""
function  DFT(f,s)
  N=size(f)
  if (s>0)
    s=1
  else
    s=-1
  end
  for k=0:N-1
    F[k+1]=0
    for j=0:N-1
      F[k+1]=F[k+1]+f[j+1]*exp(-2*s*π*im*j*k/N)
    end
  end
  return F
end
"""
   InitializeFFT(N,s)
   N:Int32
   s:Int32
   Initializes the exp(-2*s*π*i*j/N) terms necessary for an FFT
   using Algorithm 7 of Kopriva's book
"""
function InitializeFFT(N,s)
  if (s>0)
    s=1
  else
    s=-1
  end
  w=exp(-2*s*π*im/N)
  for j=0:N-1
    w_j[j+1]=w^j
  end
  return w_j
end
"""
   Radix2FFT(f,w)
   f::Array of size N
   w::Array of size N
   Computes Temperton's Radix 2 Self Sorting Complex FFT of the array f
   using the complex trigonometric factors computed by InitializeFFT
   uses Algorithm 8 of Kopriva's book
"""
function Radix2FFT(f,w)
  N=size(f)
  Ndiv2 = Int32(N/2)
  m = Int32(log2(N))
  for l=1:Int32((m+1)/2)
    noPtsAtLevel = 2^(l-1)
    a=0
    b=Ndiv2/noPtsAtLevel
    p=0
    for k=0:b-1
      W=w[p+1]
      for i=k:N-1:Int32(N/noPtsAtLevel)
        z=W*(f[a+i+1]-f[b+i+1])
        f[a+i+1]=f[a+i+1]+f[b+i+1]
        f[b+i+1]=z
      end
      p=p+noPtsAtLevel
    end
  end
  for l=(m+3)/2:m
    noPtsAtLevel=2^(l-1)
    a=0
    b=Ndiv2/noPtsAtLevel
    c=noPtsAtLevel
    d=b+noPtsAtLevel
    p=0
    for k=0:b-1
      W=w[p+1]
      for j=k:noPtsAtLevel-1:Int32(N/noPtsAtLevel)
        for i=j:N-1:2*noPtsAtLevel
          z=W*(f[a+i+1]-f[b+i+1])
          f[a+i+1]=f[a+i+1]+f[b+i+1]
          f[b+i+1]=f[c+i+1]+f[d+i+1]
          f[d+i+1]=W*(f[c+i+1]-f[d+i+1])
          f[c+i+1]=z
        end
      end
      p=p+noPtsAtLevel
    end
  end
  return f
end
"""
  FFTOfTwoRealVectors_forward(x,y,w)
  x::real array of size N
  y::real array of size N
  w::trignometric coefficients
  computes the forward FFT of two real vectors x and y using Algorithm 9 of Kopriva
"""
function  FFTOfTwoRealVectors_forward(x,y,w)
  N=size(x)
  for j =0:N-1
    Z[j+1]=x[j+1],y[j+1]im
  end
  Z=Radix2FFT(Z,w)
  X[1]=Z[1].re
  Y[1]=Z[1].im
  for k=1:N-1
    X[k+1]=0.5*(Z[k+1]+conj(Z[N-k-1])/N
    Y[k+1]=(-im/2)*(Z[k+1]+conj(Z[N-k+1]))/N
  end
  return X,Y
end
"""
    FFTOfTwoRealVectors_backward(x,y,w)
    x::real array of size N
    y::real array of size N
    w::trignometric coefficients
    computes the backward FFT of two real vectors x and y using Algorithm 10 of Kopriva
"""
function  FFTOfTwoRealVectors_backward(X,Y,w)
  N=size(x)
  for j =0:N-1
  Z[j+1]=X[j+1],Y[j+1]im
  end
  Z=Radix2FFT(Z,w)
  for k=0:N-1
    x[k+1]=Z[k+1].re
    y[k+1]=Z[k+1].im
  end
  return x,y
end
"""
   FFFTEO(f,w)
   f::real vector of size N
   wf:trigonometric coefficients for a forward transform
   computes the forward FFT of a real vector f using even-odd decomposition
   using algorithm 11 of Kopriva's book
"""
function FFFTEO(f,wf)
  N=size(f)
  for j=0:Int32(N/2)-1
    Z[j+1]=f[2*j+1]+f[2*j+2]im
    w[j+1]=wf[2*j+1]
  end
  Z=Radix2FFT(Z,w)
  F[1]=(Z[1].re+Z[1].im)/N
  F[Int32(N/2)+1]=(Z[1].re-Z[1].im)/N
  for k=1:Int32(N/2)-1
    F[k+1]=(1/(2*N))*((Z[k+1]+conj(Z[Int32(N/2)-k+1]))-wf[k+1]*(Z[k+1]-conj(Z[Int32(N/2)-k+1]))im)
  end
  for k=1:Int32(N/2)-1
    F[N-k+1]=conj(F[k+1])
  end
  return F
end
"""
   BFFTEO(F,wb)
   F:real vecotr of size N
   wb:trigonometric coefficients for a backwards transform
   computes the real backward FFT of of the Vector F using even-odd decomposition
   using algorithm 12 of Kopriva's book
"""
function BFFTEO(F,wb)
  N=size(F)
  M=Int32(N/2)
  for k=0:M-1
    E=F[k+1]+conj(F[M-k+1])
    O=wb[k+1]*(F[k+1]-conj(F[M-k+1]))
    Z[k+1]=E+O
    w[k+1]=wb[2*k+1]
  end
  Z=Radix2FFT(Z,w)
  for j=0:M-1
    f[2*j+1]=Z[j+1].re
    f[2*j+2]=Z[j+1].im
  end
  return f
end
"""
   Forward2DFFT(f,w1,w2)
   f::real matrix of size NxM
   w1,w2::forward trigonometric coefficients corresponding to N and M respectively
   computes the forward FFT of a real 2D array with an even number of points in each direction
   using Algorithm 13 of Kopriva's book
"""
function Forward2DFFT(f,w1,w2)
  N=size(f,1)
  M=size(f,2)
  for k=0:M-2:2
    F̅[:,k+1],F̅[:,k+2]=FFTOfTwoRealVectors_forward(f[:,k+1],f[:,k+2],w1)
  end
  for n=0:N-1
    F[n+1,:]=Radix2FFT(F̅[n+1,:],w2)
  end
  for m=0:M-1
    for n=0:N-1
      F[n+1,m+1]=F[n+1,m+1]/M
    end
  end
  return F
end
"""
   Backward2DFFT(F,w1,w2)
   F::complex matrix of size NxM
   w1,w2::backward trigonometric coefficients corresponding to N and M respectively
   computes the real backward FFT of a 2D array with an even number of points in each direction
   using Algorithm 14 of Kopriva's book
"""
function Backward2DFFT(F,w1,w2)
  for n=0:N-1
    F̄[n+1,:]=Radix2FFT(F[n+1,:],w2)
  end
  for k=0:M-2:2
    f[:,j+1],f[:,j+2]=FFTOfTwoRealVectors_backward(F[:,k+1],F[:,k+2],w1)
  end
  return f
end
"""
   ForwardRealFFT(x,w)
   x::real vector of size N
   w::forward trigonometric coefficients corresponding to N
   computes the forward Real transform of x using algorithm 15 of Kopriva's book
"""
function ForwardRealFFT(x,w)
  X=FFFTEO(x,w)
  for k=0:N/2
    a[k+1]=2*X[k+1].re
    b[k+1]=-2*X[k+1].im
  end
  b[1]=0
  b[Int32(N/2)+1]=0
  return a,b
end
"""
   ForwardRealFFT_efficient(x,w)
   x::real vector of size N
   w::forward trigonometric coefficients corresponding to N
   computes the forward Real transform of x using algorithm 15 of Kopriva's book
   but excludes the computation of unneccesary coefficients by FFFTEO
"""
function ForwardRealFFT_efficient(x,w)
  N=size(X)
  for j=0:Int32(N/2)-1
    Z[j+1]=x[2*j+1]+f[2*j+2]im
    w_j[j+1]=w[2*j+1]
  end
  Z=Radix2FFT(Z,w_j)
  X[1]=(Z[1].re+Z[1].im)/N
  X[Int32(N/2)+1]=(Z[1].re-Z[1].im)/N
  for k=1:Int32(N/2)-1
    X[k+1]=(1/(2*N))*((Z[k+1]+conj(Z[Int32(N/2)-k+1])-w[k+1]*(Z[k+1]-Z[Int32(N/2)-k+1])im))
  end
  for k=0:N/2
    a[k+1]=2*X[k+1].re
    b[k+1]=-2*X[k+1].im
  end
  b[1]=0
  b[Int32(N/2)+1]=0
  return a,b
end
"""
   BackwardRealFFT(a,b)
   a,b::real vectors of size N/2
   computes the backward real transform of two vectors a,b of equal size N/2
   using algorithm 16 of Kopriva's book
"""
function BackwardRealFFT(a,b)
  N2=size(a)
  w=InitializeFFT(N2*2,1)
  X[1]=a[1]/2
  X[N2]=a[N2]/2
  for k=1:N2-2
    X[k+1]=(a[k+1]+b[k+1]im)/2
  end
  for k=1:N2-2
    X[N2*2-k+1]=conj(X[k+1])
  end
  x=BFFTEO(X,w)
  return x
end
"""
   FourierDerviativeByFFT(f,wf,wb)
   f::real vector of size N
   wf::forward trigonometric coefficients
   wb::backward trigonometric coefficients
   computes the derivative of F via FFT using algorithm 17 of Kopriva's book
"""
   
function  FourierDerivativeByFFT(f,wf,wb)
  F=FFFTEO(f,wf)
  N=size(f)
  for k=0:Int32(N/2)-1
    F[k+1]=k*im*F[k+1]
  end
  F[1]=0
  for k=Int32(N/2)+1:N-1
    F[k+1]=im*(k-N)*F[k+1]
  end
  Df=BFFTEO(F,wb)
  return DF
end

"""
   FourierDerivativeMatrix(N)
   N::Integer
   computes the Fourier Derivative Matrix for an N point discretization
   using algorithm 18 of Kopriva's book
   TODO!!!!!!
   Implement sorted and ordered version if this is important for later
"""
function FourierDerivativeMatrix(N)
  for i=0:N-1
    D[i+1,i+1]=0
    for j=0:N-1
      if (j≂̸i)
        D[i+1,j+1]=0.5*(-1)^(i+j)*cot((i-j)*π/N)
        D[i+1,i+1]=D[i+1,i+1]-D[i+1,j+1]
      end
    end
  end
  return D
end
"""
   MxVDerivative(D,f)
   D::NxN Matrix
   V:vector of size N
   Computes the multiplication of V by D using algorithm 19 of Kopriva's book
   TODO!!!
   either use use a BLAS library or implement our own optimized version (BLAS much more likely) This is for demonstration purposes and will likely be abandoned
"""
function MxVDerivative(D,f)
  N=size(f)
  for i=1:N
    t=0
    for j=1:N
      t=t+D[i,j]*f[j]
    end
    Idf[i]=t
  end
  return Idf
end
"""
NOTE this ends the section on fourrier transform implementations,
TODO implement 3D routines as necessary, Test vs FFTW
Determine what can used from existing libraries FFTW, NUFFT etc.
"""

 





