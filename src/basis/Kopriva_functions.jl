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
"""
   InitializeFCosT(N)
   N::Integer
   Determines the N cosine and sine coefficients necessary for a fast
   cosine transform of an array of size N+1
"""
function InitialzeCosT(N::TInt)
  C=zeros(Float64,N+1)
  S=zeros(Float64,N+1)
  for j=0:N
    C[j+1]=cospi(j/N)
    S[j+1]=sinpi(j/N)
  end
  return C,S
end
"""
   FastCosineTransform(f,w,C,S,s)
   f::vector of size N+1
   w::trigonometric coefficients for FFT
   C::cosine coefficients
   S::Sine coefficients
   s::forward or backward determiner
   computes the fast cosine transform of a vector f using 
   algorithm 28 of Kopriva's book
"""
function FastCosineTransform(f,w,C,S,s::TInt)
  N=size(f)-1
  for j=0:N-1
    e[j+1]=0.5*(f[j+1]+f[N-j+1])-S[j+1]*(f[j+1]-f[N-j+1])
  end
  a̅,b̅ = ForwardRealFFT(e,w)
  for k=0:floor(TInt,N/2)
    a[2*k+1]= a̅[k+1]
  end
  a[2] = f[1]-f[N+1]
  for j=1:N-1
    a[2] = a[2] + 2 * C[j+1]*f[j+1]
  end
  a[2] = a[2]/N
  for k=1:floor(TInt,N/2)-1
    a[2*k+2] = b̅[k+1]+a[2*k]
  end
  if (s<0)
    for k=0:N
      a[k+1]=N*a[k+1]/2
    end
  end
  return a
end
"""
   FastChebyshevTransform(f,w,C,S,s)
   f::vector of size N+1
   w::trigonometric coefficients for FFT
   C::cosine coefficients
   S::Since coefficients
   s::forward or backward determiner
   computes the fast Chebyshev transform of a vector f using
   algorithm 29 of Kopriva's book
"""
function FastChebyshevTransform(f,w,C,S,s)
  N=size(f)-1
  for j=0:N
    g[j+1]=f[j+1]
  end
  if (s<0)
    g[1] = 2*g[1]
    g[N+1] = 2*g[N+1]
  end
  a=FastCosineTransform(g,w,C,S,s)
  if (s>0)
    a[1] = a[1]/2
    a[N+1] = a[N+1]/2
  end
  return a
end
"""
   BarycentricWeights(x)
   x:set points x[j]=x_j
   computes the Barycentric weights for a set of interpolation points x_j
   using algorithm 30 of Kopriva's book
"""
function  BarycentricWeights(x)
  N=size(x)-1
  for j=1:N+1
    w[j]=1
  end
  for j=2:N+1
    for k=1:j-1
      w[k]=w[k]*(x[k]-x[j])
      w[j]=w[j]*(x[j]-x[k])
    end
  end
  for j=1:N+1
    w[j]=1/w[j]
  end
  return w
end
"""
   LagnrangeInterpolation(x,x_j,f,w)
   x::point to interpolate to
   x_j::set of interpolation points
   f::values of given function to be interpolated at the points x_j
   w:: barycentric weights associated with x
   compute the value of the lagrange interpolation of f using the points x_j and weight w
   to onto the point x
   using Algorithm 31 of Kopriva's book
"""
function LagrangeInterpolation(x,x_j,f,w)
  num::Float64=0.0
  den::Float64=0.0
  N=size(x)-1
  for j=1:N+1
    if (AlmostEqual(x,x_j[j]))
      return f[j]
    end
    t=w[j]/(x-x_j[j])
    num = num + t*f[j]
    den = den + t
  end
  return num/den
end
"""
   PolynomialInterpolationMatrix(x,w,ξ)
   x::set of interpolation points
   w::barcyentric weights associated with x
   ξ::set of points to interpolate to
   computes the interpolation matrix needed to interpolate between the points x and ξ
   using algorithm 32 of Kopriva's book
"""
function PolynomialInterpolationMatrix(x,w,ξ)
  M=size(ξ)
  N=size(x)
  T = zeros(Float64,M,N)
  for k=1:M
    rowHasMatch = false
    for j=1:N
      T[k,j] = 0.0
      if (AlmostEquals(ξ[k],x[j]))
        rowHasMatch = true
        T[k,j] = 1.0
      end
    end
    if (rowHasMatch == false)
      s=0
      for j=1:N
        t=w[j]/(ξ[k]-x[j])
        T[k,j] = t
        s=s+t
      end
      for j =0:N
        T[k,j] = T[k,j]/s
      end
    end
  end
  return T
end
"""
   IterpolateToNewPoints(T,f)
   T::Interpolation matrix
   f::values to use for interpolation
   Interpolates f to new points using the interpolation points T
   uses Algorithm 33 of Kopriva's book
   NOTE use BLAS for this if possible
"""
function  InterpolateToNewPoints(T,f)
  M=size(T,1)
  N=size(T,2)
  for i=1:M
    t=0.0
    for j=1:N
      t=t+T[i,j]*f[j]
    end
    finterp[i]=t
  end
  return finterp
end
"""
   LagrangeInterpolatingPolynomials(x,x_j,w)
   x::point to interpolate to
   x_j::set of N interpolation points
   w::barycentric weights associated with x_j
   compute the value of N lagrangian interpolating polynomials associated 
   with x_j and w at x
   using algorithm 34 of Kopriva's book
"""
function LagrangeInterpolatingPolynomials(x,x_j,w)
  N=size(x_j)
  xMatchesNode = false
  for j=1:N
    l[j]=0.0
    if (AlmostEqual(x,x_j[j]))
      l[j]=1.0
      xMatchesNode = true
    end
  end
  if (xMatchesNode)
    return l
  end
  s=0.0
  for j=1:N
    t=w[j]/(x-x[j])
    l[j] = t
    s=s+t
  end
  for j=1:N
    l[j]=l[j]/s
  end
  return l
end
"""
   2DCoarseToFineInterpolation(x,y,f,ξ,η)
   x::x coordinates on the coarse grid
   y::y coordinates on the coarse grid
   f:: values of f(x,y)
   ξ::x coordinates on the fine grid
   η::y coordinates on the fine grid
   interpolates from coarse to a fine grid in 2D using algorithm 35 of Kopriva's book
"""
function 2DCoarseToFineInterpolation(x,y,f,ξ,η)
  N_o = size(x)
  M_o = size(y)
  N_n = size(ξ)
  M_n = size(ξ)
  w_x = BarycentricWeights(x)
  Tx = PolynomialInterpolationMatrix(x,w_x,ξ)
  for j=1:M_o
    F̅[:,j]=InterpolateToNewPoints(Tx,f[:,j])
  end
  w_y = BarycentricWeights(y)
  Ty = PolynomialInterpolationMatrix(y,w_y,η)
  for n=1:N_n
    F[n,:] = InterpolateToNewPoints(Ty,F̅[n,:])
  end
  return F
end
"""
   3DCoarseToFineInterpolation(x,y,z,f,ξ,η,ζ)
   x::x coordinates on the coarse grid
   y::y coordinates on the coarse grid
   z::z coordinates on the coarse grid
   f:: values f(x,y,z)
   ξ::x coordinates on the fine grid
   η::y coordinates on the fine grid
   ζ::z coordinates on the fine grid
   interpolates from a coarse to a fine grid in 3D extending algorithm 35 of Kopriva's book
"""
function 3DCoarseToFineInterpolation(x,y,z,f,ξ,η,ζ)
  N_o = size(x)
  M_o = size(y)
  P_o = size(z)
  N_n = size(ξ)
  M_n = size(η)
  P_n = size(ζ)
  w_x = BarycentricWeights(x)
  Tx = PolynomialInterpolationMatrix(x,w_x,ξ)
  for i=1:M_o
    for j=1:P_o
      F̅[:,i,j] = InterpolateToNewPoints(T,f[:,i,j])
    end
  end
  w_y = BarycentricWeights(y)
  Ty = PolynomialInterpolationMatrix(y,w_y,η)
  for i=1:N_n
    for j=1:P_o
      F̃[i,:,j]=InterpolateToNewPoints(Ty,F̅[i,:,j])
    end
  end
  w_z = BarycentricWeights(z)
  Tz = PolynomialInterpolationMatrix(z,w_z,ζ)
  for i=1:N_n
    for j=1:M_n
      F[i,j,:] = InterpolateToNewPoints(Tz,F̃[i,j,:])
    end
  end
  return F
end
"""
   LagrangeInterpolantDerivative(x,x_j,f,w)
   x::evaluation point
   x_j::interpolation points
   f::functions values at x_j
   w::barycentric weights associated with x_j
   evaluates the derivative of the lagrange interpolant of f at x using the point x_j of weight w
   using algorithm 36 of Kopriva's book
"""
function LagrangeInterpolantDerivative(x,x_j,f,w)
  N=size(x)
  atNode = false
  num::Float64=0.0
  for j=1:N
    if (AlmostEqual(x,x_j[j]))
      atNode = true
      p=f[j]
      den = -w[j]
      i=j
    end
  end
  if (atNode)
    for j=1:N
      if (j≂̸i)
        num = num+w[j]*(p-f[j])/(x-x_j[j])
      end
    end
  else
    den=0
    p=LagrangeInterpolation(x,N,x_j,f,w)
    for j=1:N
      t=w[j]/(x-x_j[j])
      num=num+t*(p-f[j])/(x-x_j[j])
      den = den +t
    end
  end
  return num/den
end
"""
   PolynomialDerivativeMatrix(x)
   x::Interpolation points
   computes the Polynomial derivative matrix for the set of point x
   using algorithm 37 of Kopriva's book
   !!!!TODO implement sorting to reduce round-off errors
"""
function PolynomialDerivativeMatrix(x)
  N=size(x)
  w=BarycentricWeights(x)
  for i=1:N
    D[i,i] = 0.0
    for j=1:N
      if (j≂̸i)
        D[i,j] = (w[j]/w[i])/(x[i]-x[j])
        D[i,i] = D[i,i] - D[i,j]
      end
    end
  end
  return D
end
"""
   CGLDerivativeMatrix(N)
   x::CGL points
   returns the Derivative matrix for the N Chebyshev Gauss-Lobatto points
"""
function CGLDerivativeMatrix(x,N)
  for i=1:N+1
    if (i>1 && i<N+1)
      D[i,i] = -0.5*x[i]/(sinpi(i/N))^2
    else
      D[i,i] = (2*N^2+1)/6
    end
    for j =1:N+1
      if(j≂̸i)
        if (j==1 || j == N+1)
          c̃j = 2
        else
          c̃j = 1
        end
        if (i==1 || i ==N+1)
          c̃i = 2
        else
          c̃i = 1
        end
        D[i,j] = -0.5 * (c̃i/c̃j)*(-1)^(i-1+j-1)/(sinpi((i-1+j-1)/(2*N))*sinpi((i-j)/(2*N)))
      end
    end
  end
  return D
end
"""
   mthOrderPolynomialDerivativeMatrix(m,x)
   m::order of the derivative
   x::interpolation points
   computes the mth Polynomial Derivative Matrix using the set of interpolation points x
   using algorithm 38 of Kopriva's book
"""
function MthOrderPolynomialDerivativeMatrix(m,x)
  w = BarycentricWeights(x)
  D = PolynomialDerivativeMatrix(x)
  N=size(x)
  if (m == 1)
    return D
  end
  for k=2:m
    for i=1:N
      D[i,i] = 0
      for j=1:N
        if (j≂̸i)
          D[i,i] = (k/(x[i]-x[j]))*((w[i]/w[j])*D[i,i] - D[i,j])
          D[i,i] = D[i,i] - D[i,j]
        end
      end
    end
  end
  return D
end
"""
   EOMatrixDerivative(D,f)
   D::Derivative matrix
   f:: function array
   computes the first derivative of f using Even Odd decomposition for speedup
   using algorithm 39 of Kopriva's book
"""
function EOMatrixDerivative(D,f)
  N=size(f)

  M::Int32=floor(Int32,(N)/2)
  for j=1:M+1
    e[j] = (f[j] + f[N-j])/2
    o[j] = (f[j] - f[N-j])/2
  end
  for i=1:M
    De[i] = 0.0
    Do[i] = 0.0
    for j=1:M
      De[i] = De[i] + (D[i,j] + D[i,N-j])*e[j]
      Do[i] = De[i] + (D[i,j] - D[i,N-j])*o[j]
    end
  end
  if (mod(N,2)>0)
    for i=1:M
      De[i] = De[i] + D[i,M+1]*e[M+1]
    end
    De[M+1] = 0.0
    for j=1:M
      Do[M+1] = Do[M+1] + (D[M+1,j] - D[M+1,N-j])*o[j]
    end
  end
  for j=1:M
    Df[j] = De[j] + Do[j]
    Df[N-j] = -De[j] + Do[j]
  end
  if (mod(N,2)>0)
    Df[M+1] = De[M+1] + Do[M+1]
  end
  return Df
end
"""
   FastChebyshevDerivative(f)
   f:: function array
   computes the Fast Chebysheve Derivative of a function f
   Note that this is not efficient until N>60
   at least not without building a tuned FFT
"""
function FastChebyshevDerivative(f)
  N=size(f)-1
  w=InitialzeFFT(N,1)
  C,S = InitializeCosT(N)
  f̃ = FastChebyshevTransform(f,w,C,S,1)
  f̃1 = ChebyshevDerivativeCoefficients(f̃)
  Df = FastChebyshevTransform(f̃1,w,C,S,-1)
  return Df
end





