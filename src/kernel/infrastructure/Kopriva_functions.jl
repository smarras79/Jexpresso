function DiscreteFourierCoefficients(f)
    """ DiscreteFourrierCoefficients(f)
            f:: array of length N
            computes the Discrete fourrier coefficients of f
            f̃ using Algorithm 1 of Kopriva
    """
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

function FourierInterpolantFromModes(f̃,x)
    
    """ FourrierInterpolantFromModes(f̃,x)
         f̃:: array of length N of discrete fourrier coefficients
         computes the fourrier interpolant of f at the point x using the coefficients f̃
         using Algorithm 2 of Kopriva
    """
    i=Complex(0,1)
    N=size(f̂)
    s=(f̃[1] * exp(-i * N * x / 2) + f̃[N] * exp(i * N * x * 2)) / 2
    for k = -(N/2):(N/2)-1
        s=s + f̃ * exp(i * k * x)
    end
    return Real(s)
end

function AlmostEqual(a,b)
    
    """
       AlmostEqual(a,b)
       uses algorithm 139 of Kopriva's book to
       determine if two floating point numbers a and b are nearly equal or not
    """
    #ϵ=eps(typeof(a))
    ϵ = 0.000001
    if (a == 0) || (b == 0) || (a <=ϵ) || (b <= ϵ) 
        if (abs(a-b) ≤ 2*ϵ)
            return true
        else
            return false
        end
    else
        if (abs(a-b)≤ϵ*abs(a)) && (abs(a-b)≤ϵ*abs(b))
            return true
        else
            return false
        end
    end
end

function FourierInterpolantFromNodes(f,x,x_j)
    
    """ FourrierInterpolantFromNodes(f,x,x_j)
         f:: array of length N
         x::real
         computes the fourrier interpolant of f at the point x using the coefficients the value of f at the nodes x_j
         using Algorithm 3 of Kopriva
     """
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

function LegendreDerivativeCoefficients(f̂)
    
    """
      LegendreDerivativeCoefficients(f̂)
      f̂:: coefficients due to polynomial truncation of f
      determines the coefficients of the derivative of a truncated Legendre polynomial series
      using Algorithm 4 of Kopriva's book
    """
    N=Int32(size(f,1))-1
    f̅_1[N+1]=0.0
    f̅_1[N] = (2*N-1)*f̂[N+1]
    for k=N-2:0:-1
        f̅_1[k+1]=(2*k+1)*(f̂[k+2]+f̅_1[k+3]/(2*k+5))
    end
    return f̅_1
end

function ChebyshevDerivativeCoefficients(f̂)
    
    """
       ChebyshevDerivativeCoefficients(f̂)
       f̂:: coefficients to due to polynomial truncation of f
       determines the coefficients of the derivative of a truncated Chebyshev polynomial series
       using Algorithm 5 of Kopriva's book
    """
    N=Int32(size(f,1))-1
    f̅_1[N+1]=0.0
    f̅_1[N] = 2*N*f̂[N+1]
    for k=N-2:0:-1
        f̅_1[k+1] = 2*(k+1)*f̂[k+2]+f̅_1[k+3]
    end
    f̅_1=f̂[2]+f̅_1[3]/2
    return f̅_1
end


function  DFT(f,s)
    """
       DFT(f,s)
       f::array of size N
       s::Int32
       computes the discrete forward or backwards fourrier transform of f
       depending on the sign of s using Algorithm 6 of Kopriva's book

       Convention used throughout the FFT block (kept consistent with
       Jexpresso's spectral Poisson solver, see solvers/fft_laplace.jl):
         s > 0  (forward) :  F_k = Σ_j f_j exp(-2π i j k / N)      (UNnormalized)
         s < 0  (backward):  f_j = Σ_k F_k exp(+2π i j k / N)      (UNnormalized)
       so a forward followed by a backward transform returns N·f. The 1/N
       (or 1/(N·M) in 2D) normalization is applied by the caller / by the
       Forward2DFFT wrapper. This O(N²) routine is the reference against which
       the radix-2 FFT (Radix2FFT) is verified.
    """
    N  = length(f)
    ss = s > 0 ? 1 : -1
    F  = zeros(ComplexF64, N)
    for k = 0:N-1
        acc = zero(ComplexF64)
        for j = 0:N-1
            acc += f[j+1]*exp(-2*ss*π*im*j*k/N)
        end
        F[k+1] = acc
    end
    return F
end

function InitializeFFT(N,s)
    """
       InitializeFFT(N,s)
       N:Int32
       s:Int32
       Initializes the exp(-2*s*π*i*j/N) terms necessary for an FFT
       using Algorithm 7 of Kopriva's book.

       Returns the length-N vector of twiddle factors
         w_j = exp(-2 s π i j / N),  j = 0 … N-1,
       with s collapsed to ±1 (forward / backward). These are consumed by
       Radix2FFT (the stage twiddle for a block of size m is w[k·N/m + 1]).
    """
    ss = s > 0 ? 1 : -1
    w0 = exp(-2*ss*π*im/N)
    w_j = Vector{ComplexF64}(undef, N)
    for j = 0:N-1
        w_j[j+1] = w0^j
    end
    return w_j
end

function Radix2FFT(f,w)
    """
       Radix2FFT(f,w)
       f::Array of size N (N a power of 2)
       w::length-N twiddle vector from InitializeFFT(N,s), w_j = exp(-2 s π i j/N)
       Computes the Radix-2 complex FFT of f, i.e.
         F_k = Σ_j f_j exp(-2 s π i j k / N).
       (Algorithm 8 of Kopriva's book — re-expressed here as the standard
       decimation-in-time iterative radix-2 transform, which is numerically
       robust and matches the DFT reference. The sign s is carried in by the
       precomputed twiddles `w`, so the same routine serves forward (s=+1) and
       backward (s=-1) transforms. No normalization is applied here.)
    """
    N = length(f)
    (N > 0 && (N & (N-1)) == 0) || error("Radix2FFT: N=$N must be a power of 2")
    X = Vector{ComplexF64}(undef, N)
    @inbounds for i = 1:N
        X[i] = f[i]
    end

    # ── bit-reversal permutation (decimation in time) ───────────────────────
    j = 0
    @inbounds for i = 1:N-1
        bit = N >> 1
        while (j & bit) != 0
            j ⊻= bit
            bit >>= 1
        end
        j ⊻= bit
        if i < j
            X[i+1], X[j+1] = X[j+1], X[i+1]
        end
    end

    # ── Danielson–Lanczos butterflies; stage twiddle W_m^k = w[k·N/m + 1] ───
    m = 2
    @inbounds while m <= N
        half = m >> 1
        step = N ÷ m
        for k0 = 0:m:N-1
            for k = 0:half-1
                tw = w[k*step + 1]
                a  = X[k0 + k + 1]
                b  = tw * X[k0 + k + half + 1]
                X[k0 + k + 1]        = a + b
                X[k0 + k + half + 1] = a - b
            end
        end
        m <<= 1
    end
    return X
end

function  FFTOfTwoRealVectors_forward(x,y,w)
    """
      FFTOfTwoRealVectors_forward(x,y,w)
      x::real array of size N
      y::real array of size N
      w::trignometric coefficients
      computes the forward FFT of two real vectors x and y using Algorithm 9 of Kopriva
    """
    N=size(x)
    for j =0:N-1
        Z[j+1]=x[j+1],y[j+1]im
    end
    Z=Radix2FFT(Z,w)
    X[1]=Z[1].re
    Y[1]=Z[1].im
    for k=1:N-1
        X[k+1]=0.5*(Z[k+1]+conj(Z[N-k-1]))/N
        Y[k+1]=(-im/2)*(Z[k+1]+conj(Z[N-k+1]))/N
    end
    return X,Y
end

function  FFTOfTwoRealVectors_backward(X,Y,w)
    """
        FFTOfTwoRealVectors_backward(x,y,w)
        x::real array of size N
        y::real array of size N
        w::trignometric coefficients
        computes the backward FFT of two real vectors x and y using Algorithm 10 of Kopriva
    """
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

function FFFTEO(f,wf)
    """
       FFFTEO(f,w)
       f::real vector of size N
       wf:trigonometric coefficients for a forward transform
       computes the forward FFT of a real vector f using even-odd decomposition
       using algorithm 11 of Kopriva's book
    """
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


function BFFTEO(F,wb)
    """
       BFFTEO(F,wb)
       F:real vecotr of size N
       wb:trigonometric coefficients for a backwards transform
       computes the real backward FFT of of the Vector F using even-odd decomposition
       using algorithm 12 of Kopriva's book
    """
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

function Forward2DFFT(f,w1,w2)
    """
       Forward2DFFT(f,w1,w2)
       f::real or complex matrix of size NxM
       w1,w2::FORWARD twiddles (InitializeFFT(N,1), InitializeFFT(M,1))
       Computes the NORMALIZED forward 2D FFT (the Fourier coefficients)
         F_{kx,ky} = (1/(N·M)) Σ_{nx,ny} f_{nx,ny} exp(-2π i nx kx/N) exp(-2π i ny ky/M),
       so that  f_{nx,ny} = Σ_{kx,ky} F_{kx,ky} exp(+2π i nx kx/N) exp(+2π i ny ky/M)
       is recovered by Backward2DFFT (Algorithm 13 of Kopriva's book, re-expressed
       as a separable complex transform: 1-D Radix2FFT along each column, then each
       row). The separable form is order-independent and avoids the real-pair
       packing, keeping it consistent with Jexpresso's spectral solver.
    """
    N = size(f,1)
    M = size(f,2)
    F = Matrix{ComplexF64}(undef, N, M)
    @inbounds for m = 1:M
        F[:,m] = Radix2FFT(view(f,:,m), w1)     # transform along dim 1 (columns)
    end
    @inbounds for n = 1:N
        F[n,:] = Radix2FFT(view(F,n,:), w2)     # transform along dim 2 (rows)
    end
    F ./= (N*M)
    return F
end

function Backward2DFFT(F,w1,w2)
    """
       Backward2DFFT(F,w1,w2)
       F::complex matrix of size NxM
       w1,w2::BACKWARD twiddles (InitializeFFT(N,-1), InitializeFFT(M,-1))
       Computes the (UNnormalized) inverse 2D FFT
         f_{nx,ny} = Σ_{kx,ky} F_{kx,ky} exp(+2π i nx kx/N) exp(+2π i ny ky/M),
       the exact inverse of Forward2DFFT (Algorithm 14 of Kopriva's book,
       re-expressed as a separable complex transform). Returned array is complex;
       take real() when the physical field is known to be real.
    """
    N = size(F,1)
    M = size(F,2)
    f = Matrix{ComplexF64}(undef, N, M)
    @inbounds for m = 1:M
        f[:,m] = Radix2FFT(view(F,:,m), w1)
    end
    @inbounds for n = 1:N
        f[n,:] = Radix2FFT(view(f,n,:), w2)
    end
    return f
end


function ForwardRealFFT(x,w)
    """
       ForwardRealFFT(x,w)
       x::real vector of size N
       w::forward trigonometric coefficients corresponding to N
       computes the forward Real transform of x using algorithm 15 of Kopriva's book
    """
    X=FFFTEO(x,w)
    for k=0:N/2
        a[k+1]=2*X[k+1].re
        b[k+1]=-2*X[k+1].im
    end
    b[1]=0
    b[Int32(N/2)+1]=0
    return a,b
end


function ForwardRealFFT_efficient(x,w)
    """
       ForwardRealFFT_efficient(x,w)
       x::real vector of size N
       w::forward trigonometric coefficients corresponding to N
       computes the forward Real transform of x using algorithm 15 of Kopriva's book
       but excludes the computation of unneccesary coefficients by FFFTEO
    """
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


function BackwardRealFFT(a,b)
    """
       BackwardRealFFT(a,b)
       a,b::real vectors of size N/2
       computes the backward real transform of two vectors a,b of equal size N/2
       using algorithm 16 of Kopriva's book
    """
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

function  FourierDerivativeByFFT(f,wf,wb)
    """
       FourierDerviativeByFFT(f,wf,wb)
       f::real vector of size N
       wf::forward trigonometric coefficients
       wb::backward trigonometric coefficients
       computes the derivative of F via FFT using algorithm 17 of Kopriva's book
    """

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

function FourierDerivativeMatrix(N)
    """
       FourierDerivativeMatrix(N)
       N::Integer
       computes the Fourier Derivative Matrix for an N point discretization
       using algorithm 18 of Kopriva's book
       TODO!!!!!!
       Implement sorted and ordered version if this is important for later
    """
    D=zeros(Float64,N,N)
    for i=0:N-1
        D[i+1,i+1]=0
        for j=0:N-1
            if (j != i)
                D[i+1,j+1]=0.5*(-1)^(i+j)*cot((i-j)*π/N)
                D[i+1,i+1] = D[i+1,i+1] - D[i+1,j+1]
            end
        end
    end
    @info D
    return D
end


function MxVDerivative(D,f)
    """
       MxVDerivative(D,f)
       D::NxN Matrix
       V:vector of size N
       Computes the multiplication of V by D using algorithm 19 of Kopriva's book
       TODO!!!
       either use use a BLAS library or implement our own optimized version (BLAS much more likely) This is for demonstration purposes and will likely be abandoned
    """
    N=size(f,1)
    Idf = zeros(Float64,N)
    for i=1:N
        t=0.0
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

function InitializeCosT(N)
    """ 
    InitializeFCosT(N)
       N::Integer
       Determines the N cosine and sine coefficients necessary for a fast
       cosine transform of an array of size N+1
    """
    C=zeros(Float64,N+1)
    S=zeros(Float64,N+1)
    for j=0:N
        C[j+1]=cospi(j/N)
        S[j+1]=sinpi(j/N)
    end
    return C,S
end

function FastCosineTransform(f,w,C,S,s)
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

function FastChebyshevTransform(f,w,C,S,s)
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


function  BarycentricWeights(x)
    """
       BarycentricWeights(x)
       x:set points x[j]=x_j
       computes the Barycentric weights for a set of interpolation points x_j
       using algorithm 30 of Kopriva's book
    """
    N=size(x,1)-1
    w=zeros(Float64,N+1)
    mini = minimum(x)
    maxi = maximum(x)
    scale = (maxi - mini)/2
    for j=1:N+1
        w[j]=1
    end
    for j=2:N+1
        for k=1:j-1
            w[k] = w[k]*(x[k]-x[j])/scale
            w[j] = w[j]*(x[j]-x[k])/scale
        end
    end
    for j=1:N+1
        w[j]=1/w[j]
    end
    return w
end

function  BarycentricWeights!(x,w)
    """
       BarycentricWeights(x)
       x:set points x[j]=x_j
       computes the Barycentric weights for a set of interpolation points x_j
       using algorithm 30 of Kopriva's book
    """
    N=size(x,1)-1
    mini = minimum(x)
    maxi = maximum(x)
    scale = (maxi - mini)/2
    for j=1:N+1
        w[j]=1
    end
    for j=2:N+1
        for k=1:j-1
            w[k] = w[k]*(x[k]-x[j])/scale
            w[j] = w[j]*(x[j]-x[k])/scale
        end
    end
    for j=1:N+1
        w[j]=1/w[j]
    end
    return w
end


function BarycentricWeights_gpu!(x,N,ω)
    j = @index(Global, Linear)
    T = eltype(x)
    ω[j] = T(1)
    for k=1:N+1
        if (k ≠ j)
            ω[j] = ω[j]*(x[j] - x[k])
        end
    end
    ω[j] = T(1)/ω[j]

end


function LagrangeInterpolation(x,x_j,f,w)
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
    num::Float64=0.0
    den::Float64=0.0
    N=size(x_j,1)-1

    for j=1:N+1
        if (AlmostEqual(x,x_j[j]))
            return f[j]
        end
        t=w[j]/(x-x_j[j])+eps(Float64)
        num = num + t*f[j]
        den = den + t
    end
    return num/den
end


function PolynomialInterpolationMatrix(x,w,ξ)
    """
       PolynomialInterpolationMatrix(x,w,ξ)
       x::set of interpolation points
       w::barcyentric weights associated with x
       ξ::set of points to interpolate to
       computes the interpolation matrix needed to interpolate between the points x and ξ
       using algorithm 32 of Kopriva's book
    """
    M=size(ξ,1)
    N=size(x,1)
    mini = minimum(x)
    maxi = maximum(x)
    scale = (maxi - mini)/2
    T = zeros(Float64,M,N)
    for k=1:M
        rowHasMatch = false
        for j=1:N
            T[k,j] = 0.0
            if (AlmostEqual(ξ[k],x[j]))
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
            for j =1:N
                T[k,j] = T[k,j]/s
            end
        end
    end
    return T
end

function PolynomialInterpolationMatrix!(x,w,ξ,T)
    """
       PolynomialInterpolationMatrix(x,w,ξ)
       x::set of interpolation points
       w::barcyentric weights associated with x
       ξ::set of points to interpolate to
       computes the interpolation matrix needed to interpolate between the points x and ξ
       using algorithm 32 of Kopriva's book
    """
    M=size(ξ,1)
    N=size(x,1)
    mini = minimum(x)
    maxi = maximum(x)
    scale = (maxi - mini)/2
    for k=1:M
        rowHasMatch = false
        for j=1:N
            T[k,j] = 0.0
            if (AlmostEqual(ξ[k],x[j]))
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
            for j =1:N
                T[k,j] = T[k,j]/s
            end
        end
    end
    return T
end


function  InterpolateToNewPoints(T,f)
    """
       IterpolateToNewPoints(T,f)
       T::Interpolation matrix
       f::values to use for interpolation
       Interpolates f to new points using the interpolation points T
       uses Algorithm 33 of Kopriva's book
       NOTE use BLAS for this if possible
    """
    M=size(T,1)
    N=size(T,2)
    finterp=zeros(Float64,M)
    for i=1:M
        t=0.0
        for j=1:N
            t=t+T[i,j]*f[j]
        end
        finterp[i]=t
    end
    return finterp
end



function LagrangeInterpolatingPolynomials(x,x_j,w)
    """
       LagrangeInterpolatingPolynomials(x,x_j,w)
       x::point to interpolate to
       x_j::set of N interpolation points
       w::barycentric weights associated with x_j
       compute the value of N lagrangian interpolating polynomials associated 
       with x_j and w at x
       using algorithm 34 of Kopriva's book
    """
    N=size(x_j,1)
    l=zeros(Float64,N)
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
        t=w[j]/(x-x_j[j])
        l[j] = t
        s=s+t
    end
    for j=1:N
        l[j]=l[j]/s
    end
    return l
end

function LagrangeInterpolatingPolynomials!(x,x_j,w,l)
    """
       LagrangeInterpolatingPolynomials(x,x_j,w)
       x::point to interpolate to
       x_j::set of N interpolation points
       w::barycentric weights associated with x_j
       compute the value of N lagrangian interpolating polynomials associated 
       with x_j and w at x
       using algorithm 34 of Kopriva's book
    """
    N=size(x_j,1)
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
        t=w[j]/(x-x_j[j])
        l[j] = t
        s=s+t
    end
    for j=1:N
        l[j]=l[j]/s
    end
    return l
end


function D2CoarseToFineInterpolation(x,y,f,ξ,η)
    """
       2DCoarseToFineInterpolation(x,y,f,ξ,η)
       x::x coordinates on the coarse grid
       y::y coordinates on the coarse grid
       f:: values of f(x,y)
       ξ::x coordinates on the fine grid
       η::y coordinates on the fine grid
       interpolates from coarse to a fine grid in 2D using algorithm 35 of Kopriva's book
    """
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

function D3CoarseToFineInterpolation(x,y,z,f,ξ,η,ζ)
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
            F̅[:,i,j] = InterpolateToNewPoints(Tx,f[:,i,j])
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


function LagrangeInterpolantDerivative(x,x_j,f,w)
    """
       LagrangeInterpolantDerivative(x,x_j,f,w)
       x::evaluation point
       x_j::interpolation points
       f::functions values at x_j
       w::barycentric weights associated with x_j
       evaluates the derivative of the lagrange interpolant of f at x using the point x_j of weight w
       using algorithm 36 of Kopriva's book
    """
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


function PolynomialDerivativeMatrix(x)
"""
       PolynomialDerivativeMatrix(x)
       x::Interpolation points
       computes the Polynomial derivative matrix for the set of point x
       using algorithm 37 of Kopriva's book
       !!!!TODO implement sorting to reduce round-off errors
    """
    N=size(x,1)
    w=BarycentricWeights(x)
    D=zeros(Float64,N,N)
    for i=1:N
        D[i,i] = 0.0
        for j=1:N
            if (j!=i)
                D[i,j] = (w[j]/w[i])/(x[i]-x[j])
            end
        end
        Ds = sort(D[i,:],by=abs)
        for j=1:N
            D[i,i] = D[i,i] - Ds[j]
        end
    end
    return D
end


function CGLDerivativeMatrix(x,N)
"""
       CGLDerivativeMatrix(x,N)
       x::CGL points x_j = cos(jπ/N), j=0…N  (length N+1)
       returns the (N+1)×(N+1) Chebyshev–Gauss–Lobatto first-derivative matrix.

       Corrected from Kopriva's listing: allocate D; the off-diagonals use the
       numerically-stable sin product form of 1/(x_i-x_j) (verified — Kopriva's
       formula); the diagonal is taken from the NEGATIVE-SUM TRICK (rows of an
       exact derivative matrix sum to zero, since d/dx of a constant is 0), which
       reproduces the endpoint values ±(2N²+1)/6 with the correct signs and is
       more robust than the explicit closed forms (the original set both endpoint
       diagonals to +(2N²+1)/6 and used the wrong index in the interior term).
    """
    D = zeros(Float64, N+1, N+1)
    for i = 1:N+1
        c̃i = (i == 1 || i == N+1) ? 2.0 : 1.0
        for j = 1:N+1
            if j != i
                c̃j = (j == 1 || j == N+1) ? 2.0 : 1.0
                D[i,j] = -0.5 * (c̃i/c̃j) * (-1)^(i-1+j-1) /
                         (sinpi((i-1+j-1)/(2*N)) * sinpi((i-j)/(2*N)))
            end
        end
    end
    for i = 1:N+1                       # diagonal: negative-sum trick
        s = 0.0
        for j = 1:N+1
            j != i && (s += D[i,j])
        end
        D[i,i] = -s
    end
    return D
end


function mthOrderPolynomialDerivativeMatrix(m,x)
"""
       mthOrderPolynomialDerivativeMatrix(m,x)
       m::order of the derivative
       x::interpolation points
       computes the mth Polynomial Derivative Matrix using the set of interpolation points x
       using algorithm 38 of Kopriva's book
    """
    w = BarycentricWeights(x)
    Dtemp = PolynomialDerivativeMatrix(x)
    N=size(x,1)
    if (m == 1)
        return Dtemp
    end
    Dm = zeros(Float64,N,N)
    for k=2:m 
        for i=1:N
            Dm[i,i] = 0
            for j=1:N
                if (j!=i)
                    Dm[i,j] = (k/(x[i]-x[j]))*((w[j]/w[i])*Dtemp[i,i] - Dtemp[i,j])
                    Dm[i,i] = Dm[i,i] - Dm[i,j]
                end
            end
        end
        Dtemp.=Dm
    end
    return Dtemp
end


function EOMatrixDerivative(D,f)
"""
       EOMatrixDerivative(D,f)
       D::Derivative matrix
       f:: function array
       computes the first derivative of f using Even Odd decomposition for speedup
       using algorithm 39 of Kopriva's book
    """
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


function FastChebyshevDerivative(f)
"""
       FastChebyshevDerivative(f)
       f:: function array
       computes the Fast Chebysheve Derivative of a function f
       Note that this is not efficient until N>60
       at least not without building a tuned FFT
    """
    N=size(f)-1
    w=InitialzeFFT(N,1)
    C,S = InitializeCosT(N)
    f̃ = FastChebyshevTransform(f,w,C,S,1)
    f̃1 = ChebyshevDerivativeCoefficients(f̃)
    Df = FastChebyshevTransform(f̃1,w,C,S,-1)
    return Df
end

function  FourierCollocationTimeDerivative(Φ,D)
"""
       FourierCollocationTimeDerivative(Φ,D)
       Φ::Approximated function
       D::Derivative matrix
       ν::Diffusion coefficient
       Computes the time derivative for advection diffusion equation
       using the Fourrier Collocation method
       Algorithm 41 of Koopriva's book
    """
    F=MxVDerivative(D,Φ)
    N=size(Φ,1)
    for j=1:N
        F[j] = ν*F[j]-Φ[j]
    end
    Φt = MxVDerivative(D,F)
    return Φt
end

function CollocationStepByRK3(tn,Δt,Φ,D,ν,a,b,g)
"""
       CollocationStepByRK3(tn,Δt,Φ,D,ν,a,b,g)
       tn::current time
       Δt:: time step
       Φ:: Approximated function
       D:: Derivative matrix
       ν:: Diffusion coefficient
       a,b,g:: RK3 tableau
       Advance the advection diffusion equation in time using RK3
       Algorithm 42 of Kopriva's book
    """
    N=size(Φ,1)
    G=zeros(Float64,N)
    for m=1:3
        t=tn + b[m] * Δt
        Φt = FourierCollocationTimeDerivative(Φ,D)
        for j=1:N
            G[j] = a[m]*G[j] + Φt[j]
            Φ[j] = Φ[j] + g[m]*Δt*G[j]
        end
    end
    return Φ
end

function  FourierCollocationDriver(N,NT,T,Φ,ν,a,b,g)
"""
       FourierCollocationDriver(N,NT,T,Φ,ν,a,b,g)
       N::Number of collocation points
       NT::number of time steps
       T:: end time
       Φ::Approximated function
       ν::diffusion coefficients
       a,b,g::RK3 tableau
       A driver for the Fourrier Collocation method using RK3
       Algorithm 43 of Kopriva's book
    """

    N1=size(Φ,1)
    D=FourierDerivativeMatrix(N1)
    Δt = T/NT
    tn =0.0
    for n=0:NT-1
        Φ = CollocationStepByRK3(tn,Δt,Φ,D,ν,a,b,g)
        tn = (n+1)*Δt
    end
    return Φ
end



function ADTimeDerivative(Φ̂,ν)
"""
       ADTimeDerivative(Φ̂)
       Φ̂::set of N+1 fourrier coefficients of Φ
       ν::diffusion coefficient
       Computes the time derivative of the fourrier coefficients for the solution Φ of an advection diffusion equation
       uses Algorithm 44 of Kopriva's book
    """

    N=size(Φ̂,1)-1
    N2=floor(Int64,N/2)
    Φ̂t = zeros(Complex,N+1)
    for k=1:N+1
        iter = k-1-N2
        Φ̂t[k] = -(im*iter + ν*iter^2)*Φ̂[k]
    end
    return Φ̂t
end

function FourrierGalerkinStep(tn,Δt,Φ̂,ν,a,b,g)
"""
       FourrierGalerkinStep(tn,Δt,Φ̂,ν,a,b,g)
       Φ̂::set of N+1 Fourrier coefficients of Φ
       tn::current time
       ν::diffusion coefficient
       Δt::time step
       a,b,g::RK3 tableau
       Advance the advection diffusion equation in time using the RK3 time integrator
       Algorithm 45 of Kopriva's book
    """

    N=size(Φ̂,1)
    G=zeros(Complex,N)
    for m=1:3
        t=tn+b[m]*Δt
        Φ̂t=ADTimeDerivative(Φ̂,ν)
        for k=1:N
            G[k] = a[m]*G[k] + Φ̂t[k]
            Φ̂[k] = Φ̂[k] + g[m]*Δt*G[k]
        end
    end
    return Φ̂
end

function EvaluateFourierGalerkinSolution(x,Φ̂)
"""
       EvaluateFourierGalerkinSolution(x,Φ̂)
       x::evaluation point
       Φ̂::Set of N+1 Fourrier coefficients of Φ
       evaluates the solution Φ at the point x using the Fourrier Coefficients Φ̂
       Algorithm 46 of Kopriva's book
    """

    N=size(Φ̂,1)-1
    N2=floor(Int64,N/2)
    Φ::Complex=0.0
    for k=1:N+1
        iter=k-1-N2
        Φ = Φ + Φ̂[k]*exp(im*iter*x)
    end
    return Φ
end

function FourierGalerkinDriver(Φ̂,N,NT,T,Nout,ν,a,b,g)
"""
       FourierGalerkinDriver(Φ̂,N,NT,T,Nout,ν,a,b,g)
       Φ̂:: initial Fourrier coefficients of Φ
       N:: N number of fourrier coefficients of truncation
       NT:: number of time steps
       T:: end time
       ν:: diffusion coefficient
       a,b,g:: RK3 tableau
       computes the solution of the advection diffusion equation at time T using the RK3
       time integrator following Fourrier Galerkin method
       Algorithm 47 of Kopriva's book
    """

    Δt = T/NT
    tn = 0.0
    for n=0:NT-1
        Φ̂=FourrierGalerkinStep(tn,Δt,Φ̂,ν,a,b,g)
        tn=(n+1)*Δt
    end
    Δx = 2*π/Nout
    x=zeros(Float64,Nout+1)
    Φ=zeros(Complex,Nout+1)
    for j=0:Nout
        x[j+1] = j*Δx
        Φ[j+1] = EvaluateFourierGalerkinSolution(x[j+1],Φ̂)
    end
      @info Φ
    
    return Φ
end
#TODO build driver for scalar advection

function  DirectConvolutionSum(V̂,Ŵ)
"""
       DirectConvolutionSum(V̂,Ŵ)
       V̂::Fourrier Coefficients of V
       Ŵ::Fourrier Coeffiicients of W
       computes the Direct convolutions sum for V and W
       Algorithm48 of Kopriva's book
    """

    N=size(V̂,1)-1
    N2=floor(Int64,N/2)
    VW=zeros(Complex,N+1)
    for k=-N2:N2
        iter=k+1+N2
        VW[iter]=0.0
        for p=max(-N2,k-N2):min(N2,N2+k)
            iter1=p+1+N2
            VW[iter] = VW[iter] + V̂[iter-iter1]*Ŵ[iter1]
        end
    end
    return VW
end

function FastConvolutionSum(V̂,Ŵ)
"""
       FastConvolutionSum(V̂,Ŵ)
       V::Fourrier Coefficients of V
       W::Fourrier Coefficients of W
       computes the Fast Convolution sum for V and W
       Algorithm 49 of Kopriva's book
    """

    N=size(V̂,1)-1
    N2=floor(Int64,N/2)
    M=2*N
    M2 = floor(Int64,M/2)
    Ṽ=zeros(Complex,M+1)
    W̃=zeros(Complex,M+1)
    for k=0:N2
        Ṽ[k+1] = V̂[k+1]
        W̃[k+1] = Ŵ[k+1]
    end
    for k=-1:-N2/2:-1
        Ṽ[M+1+k] = V̂[k+1+N]
        W̃[M+1+k] = Ŵ[k+1+N]
    end
    w=InitializeFFT(M,-1)
    V = Radix2FFT(Ṽ,w)
    W = Radix2FFT(W̃,w)
    Q=zeros(M)
    for j=0:M-1
        Q[j+1] = V[j+1]*W[j+1]
    end
    w=InitializeFFT(M,1)
    Q̃=Radix2FFT(Q,w)
    for k=0:N2
        VW[k+1] = Q̃[k+1]
    end
    for k=-1:-N2:-1
        VW[k+1+N] = Q̃[M+1+k]
    end
    return VW
end

function  CollocationStepByRK3(tn,Δt,Φ,D,D2,gL,gR,a,b,g,o)
"""
       CollocationStepByRK3(tn,Δt,Φ,D,gL,gR)
       tn::current time
       Δt::time step
       Φ::Approximated function
       D::Derivative matrix
       gL::Left side boundary condition
       gR::right side boundary condition
       a,b,g::RK3 tableau
       o::2 for diffusion, 1 for advection, 3 for both
       Advance the diffusion equation in time using RK3 and the collocation method
       Algorithm 50 of Kopriva's book
    """

    N=size(Φ,1)
    G=zeros(Float64,N)
    Φt=zeros(Float64,N)
    for m=1:3
        t=tn+b[m]*Δt
        Φt = MxVDerivative(D2,Φ) - MxVDerivative(D,Φ)
        for j=1:N
            G[j] = a[m]*G[j] + Φt[j]
            Φ[j] = Φ[j] + g[m]*Δt*G[j]
        end
    end
    if (o==1)
        Φ[1] = sinpi(-tn-Δt)
        Φ[N] = sinpi(2-tn-Δt)
    elseif (o==3)
        Φ[1] = sinpi(-tn-Δt)*exp(-π^2*(tn+Δt))
        Φ[N] = sinpi(2-tn-Δt)*exp(-π^2*(tn+Δt))
    else
        Φ[1] = gL
        Φ[N] = gR
    end
    return Φ
end


function  LegendreCollocationIntegrator(N,NT,Nout,T,Φ,a,b,g,gL,gR,ad)
"""
       LegendreCollocationIntegrator(N,NT,Nout,T,Φ,a,b,g)
       N::Number of points
       NT:number of time steps
       Nout::number of output points
       T::end time
       Φ:initial conditions
       a,b,g::RK3 tableau
       gL,gR:: Left and right boundary conditions
       c::Advection coefficient
       m::2 for diffusion, 1 for advection, 3 for both
       ν::diffusion coefficient
       solves the diffusion equation using a Legendre Collocation integrator operating on RK3
       Algorithm 51 of Kopriva's book
    """
    Legendre = St_Legendre{TFloat}(0.0, 0.0, 0.0, 0.0)
    lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    LegendreGaussLobattoNodesAndWeights!(Legendre,lgl,N)
    x=lgl.ξ
    w=lgl.ω
    D=zeros(Float64,N+1,N+1)
    D2=zeros(Float64,N+1,N+1)
    if (ad==3)
        D=PolynomialDerivativeMatrix(x)
        D2 = mthOrderPolynomialDerivativeMatrix(2,x)
    elseif(ad==2)
        D2 = mthOrderPolynomialDerivativeMatrix(2,x)
    else
        D = PolynomialDerivativeMatrix(x)
    end
    Δt = T/NT
    tn=0.0
    for n=0:NT-1
        Φ=CollocationStepByRK3(tn,Δt,Φ,D,D2,gL,gR,a,b,g,ad)
        tn =(n+1)*Δt
    end
    X=zeros(Float64,Nout+1)
    for j=0:Nout
        X[j+1]=-1+2*j/Nout
    end
    w=BarycentricWeights(x)
    T=zeros(Float64,Nout+1,N)
    T=PolynomialInterpolationMatrix(x,w,X)
    ΦI = zeros(Float64,Nout+1)
    ΦI = InterpolateToNewPoints(T,Φ)
    return ΦI
end


function ModifiedLegendreBasis(nop,x)
"""
       ModifiedLegendreBasis(nop,x)
       nop::PolynomialOrder
       x::EvaluationPoint
       Evaluate the LegendreModified basis of order nop at point x
       Algorithm 52 of Kopriva's book
    """

    Legendre = St_Legendre{TFloat}(0.0, 0.0, 0.0, 0.0)
    LegendreAndDerivativeAndQ!(Legendre,nop,x)
    L1 = Legendre
    LegendreAndDerivativeAndQ!(Legendre,nop+2,x)
    L2= Legendre
    Φ=L1.legendre - L2.legendre
    Φ=Φ/sqrt(4*nop+6)
end

function EvaluateLegendreGalerkinSolution(N,x,Φ̂)
"""
       EvaluateLegendreGalerkinSolution(N,x,Φ̂)
       N::polynomial order
       x::Evaluation point
       Φ̂::N-2 Coefficients
       Evaluate the Legendre Galerkin solution at the point x
       Algorithm 53 of Kopriva's book
    """

    Φ=0.0
    for k=0:N-2
        Φ = Φ +Φ̂[k+1]*ModifiedLegendreBasis(N,x)
    end
    return Φ
end

function InitMatrix(N,p)
"""
       InitTMatrix(N,p)
       N::polynomial order
       p::Even or odd index 0 or 1
       Computes even or odd index coefficients of the tridiagonal Legendre Galerkin matrix
       Algorithm 54 of Kopriva's book
    """

    d=zeros(Float64,N+1)
    l=zeros(Float64,N)
    u=zeros(Float64,N)
    for j=1:N+1
        d[j] =(1/sqrt(4*(2*(j-1)+p)+6))^2*(-1)*((-2/(2*(2*(j-1)+p)+1))+(-2/(2*(2*(j-1)+p)+5)))
    end
    for i=2:N+1
        l[i-1] = (-2/(2*(2*(i-1)+p)+1))*(1/sqrt(4*(2*(j-1)+p)+6))*(1/sqrt(4*(2*(j-2)+p)+6))
        u[i-1] = l[i-1]
    end
    return d,l,u
end

function  TriDiagonalSolve(l,d,u,y)
"""
       TriDiagonalSolve(l,d,u,y)
       l::lower diagonal
       d::diagonal
       u::upper diagonal
       y:: rhs vector
       Invert the Tridiagonal matrix of diagonal d, lower diagonal l and u
       Algorithm 141 of Kopriva's book
    """

    N=size(d,1)
    d̂=zeros(Float64,N)
    for j=1:N
        d̂[j] = d[j]
    end
    for j=2:N
        d̂[j] = d̂[j] - l[j-1]/d̂[j-1]*u[j-1]
        y[j] = y[j] - l[j-1]/d̂[j-1]*y[j-1]
    end
    x=zeros(Float64,N)
    x[N]=y[N]/d̂[N]
    for j=N-1:1:-1
        x[j] = (y[j]-u[j]*x[j+1])/d̂[j]
    end
    return x
end

function  ModifiedCoefsFromLegendreCoefs(ϕ̂)
"""
       ModifiedCoefsFromLegendreCoefs(Φ̂)
       ϕ̂::Legendre expansion coefficients
       Computes the modified Legendre coefficients to verify boundary conditions
       Algorithm 55 of Kopriva's book
    """

    N=size(ϕ̂,1)-1
    M=floor(Int64,N/2)
    d,l,u = InitMatrix(M+1,0)
    rhs= zeros(Float64,M+1)
    b=zeros(Float64,M+1)
    for j =1:M+1
        iter = j-1
        rhs[j] = (-2/(4*iter+5))*(1/sqrt(8*iter+6))*ϕ̂[2*iter+3] - (1/sqrt(8*iter+6)) * (-1)*((-2/(4*iter+5))+(-2/(4*iter+1))) * ϕ̂[2*iter+1]
    end
    b=TriDiagonalSolve(l,d,u,rhs)
    for j =0:M
        Φ̂[2*j+1] = b[j+1]
    end
    M=floor(Int64,N/2)-1
    d,l,u = InitMatrix(M,1)
    for j =1:M+1
        iter = j-1
        rhs[j] = (-2/(2*(2*iter+1)+5))*(1/sqrt(4*(2*iter+1)+6))*ϕ̂[2*iter+4] - (1/sqrt(4*(2*iter+1)+6)) * (-1)*((-2/(2*(2*iter+1)+5))+(-2/(2*(2*iter+1)+1))) * ϕ̂[2*iter+1]
    end
    b = TriDiagonalSolve(l,d,u,rhs)
    for j =0:M
        Φ̂[2*j+2] = b[j+1]
    end
    return Φ̂
end

function LegendreGalerkinStep(Δt,Φ̂)
"""
       LegendreGalerkinStep(Δt,Φ̂n)
       Δt::time step
       Φn::Legendre Coefficients at time n
       Advance one time step using the Trapezoidal rule using the Legendre Galerkin method
       Algorithm 56 of Kopriva's book
    """

    N=size(Φ̂,1)
    #Even indexed coefficients
    M=floor(Int64,(N-1)/2)
    d,l,u = InitMatrix(M+1,0)
    rhs= zeros(Float64,M+1)
    rhs[1] = (d[1]-Δt/2)*Φ̂n[1] + u[1]*Φ̂n[3]
    for j=2:M
        iter = j-1
        rhs[j] = l[iter]*Φ̂n[2*(iter-1)+1]+(d[j]-Δt/2)*Φ̂n[2*iter+1] + u[j] * Φ̂n[2*j]
    end
    rhs[M+1] = (d[M+1]-Δt/2)*Φ̂n[2*M+1]+l[M]*Φn[2*(M-1)+1]
    for j=1:M+1
        d[j] = d[j] + Δt/2
    end
    Φ̂ = TriDiagonalSolve(l,d,u,rhs)
    Φ̂n1 = zeros(Float64,N)
    for j=0:M
        Φ̂n1[2*j+1] = Φ̂[j+1]
    end
    #Odd indexed coefficients
    M=floor(Int64,(N-1)/2)-1
    d,l,u = InitMatrix(M+1,1)
    rhs[1] = (d[1]-Δt/2)*Φ̂n[2] + u[1]*Φ̂n[4]
    for j=2:M
        iter = j-1
        rhs[j] = l[iter]*Φ̂n[2*(iter-1)+2]+(d[j]-Δt/2)*Φ̂n[2*iter+2] + u[j] * Φ̂n[2*j+1]
    end
    rhs[M+1] = (d[M+1]-Δt/2)*Φ̂n[2*M+2]+l[M]*Φn[2*(M-1)+2]
    for j=1:M+1
        d[j] = d[j] + Δt/2
    end
    Φ̂ = TriDiagonalSolve(l,d,u,rhs)
    for j=0:M
        Φ̂n1[2*j+2] = Φ̂[j+1]
    end
    return Φ̂n1
end

function LegendreGalerkinDriver(N,NT,T,Nout,Φ̂,a,b,g)
"""
       LegendreGalerkinDriver((N,NT,T,Nout,Φ,a,b,g)
       N::Order of the approximation
       NT::number of time steps
       Nout::number of output points
       T::final time
       Φ̂:: approximated solution coefficients
       a,b,g:: RK3 tableau
       Driver for the Legendre Galerkin method
    """

    Δt = T/Nout
    tn = 0.0
    for n=0:NT-1
        Φ̂=LegendreGalerkinStep(Δt,Φ̂)
        tn = (n+1)*Δt
    end
    Δx=2*π/Nout
    x=zeros(Float64,Nout+1)
    Φ=zeros(Float64,Nout+1)
    for j =0:Nout
        x[j+1] = j*Δx
        Φ[j+1] = EvaluateLegendreGalerkinSolution(x[j+1],Φ̂)
    end
    return Φ
end


function  CGDerivativeMatrix(N)
"""
       CGDerivativeMatrix(N)
       N::Polynomial order
       Generates the Derivative matrix for a Continuous Galerkin method
       Algorithm 57 of Kopriva's book
    """

    Legendre = St_Legendre{TFloat}(0.0, 0.0, 0.0, 0.0)
    lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    LegendreGaussLobattoNodesAndWeights!(Legendre, lgl,N)
    x=lgl.ξ
    w=lgl.ω
    D = PolynomialDerivativeMatrix(x)
    G=zeros(Float64,N+1,N+1)
    for j=1:N+1
        for n=1:N+1
            s=0.0
            for k=1:N+1
                s=s+D[k,n]*D[k,j]*w[k]
            end
            G[j,n]=s#/w[j]
        end
    end
    return G
end

function CGDerivativeMatrixĜ(N)
    Legendre = St_Legendre{TFloat}(0.0, 0.0, 0.0, 0.0)
    lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    LegendreGaussLobattoNodesAndWeights!(Legendre, lgl,N)
    x=lgl.ξ
    w=lgl.ω
    D = PolynomialDerivativeMatrix(x)
    G=zeros(Float64,N+1,N+1)
    for j=1:N+1
        for n=1:N+1
            s=0.0
            for k=1:N+1
                s=s+D[k,n]*D[k,j]*w[k]
            end
            G[j,n]=s/w[j]
        end
    end
    return G
end

function  CGDriver(N,NT,Nout,T,Φ,a,b,g,gL,gR)
"""
       CGDriver(N,NT,Nout,T,Φ,a,b,g,gL,gR)
       N::Order of approximation
       NT::number of time steps
       Nout::Number of outputs
       Φ::approximated solution
       a,b,g::RK3 tableau
       gL,gR::boundary conditions
       CG Driver
    """
    D = zeros(Float64,N+1,N+1)
    D2 = -CGDerivativeMatrixĜ(N)
    Δt = T/NT
    tn=0.0
    for n=0:NT-1
        Φ=CollocationStepByRK3(tn,Δt,Φ,D,D2,gL,gR,a,b,g,2)
        tn =(n+1)*Δt
    end
    X=zeros(Float64,Nout+1)
    for j=0:Nout
        X[j+1]=-1+2*j/Nout
    end
    w=BarycentricWeights(x)
    T=zeros(Float64,Nout+1,N)
    T=PolynomialInterpolationMatrix(x,w,X)
    ΦI = zeros(Float64,Nout+1)
    ΦI = InterpolateToNewPoints(T,Φ)
    return ΦI
end


mutable struct NodalDiscontinuousGalerkin
"""
       NodalDiscontinuousGalerkin
       N::Polynomial order
       D̂::Derivative matrix
       lL:: left boundary lagrange interprolant
       lR:: right boundary lagrange interpolant
       Φ:: Solution Array
       A class for the Nodal Discontinuous Galerkin method
       Algorithm 58 of Kopriva's book
    """
    N::TInt
    D̂::Array{TFloat}
    lL::Array{TFloat}
    lR::Array{TFloat}
    w::Array{TFloat}
    Φ::Array{TFloat}
end

function  buildNodalDiscontinuousGalerkin!(N,DG::NodalDiscontinuousGalerkin)
"""
       BuildNodalDiscontinuousGalerkin(N)
       N::Polynomial order
       Construct the NodalDiscontinuousGalerkin struct
       Algorithm 59 of Kopriva's book
    """

    DG.N = N
    Legendre = St_Legendre{TFloat}(0.0, 0.0, 0.0, 0.0)
    lgl      = St_lgl{TFloat}(zeros(TFloat, N+1),
                              zeros(TFloat, N+1))
    LegendreGaussLobattoNodesAndWeights!(Legendre,lgl,N)
    x=lgl.ξ
    w=lgl.ω
    DG.w = w
    wB = BarycentricWeights(x)
    DG.lL = LagrangeInterpolatingPolynomials(-1.0,x,wB)
    DG.lR = LagrangeInterpolatingPolynomials(1.0,x,wB)
    D= PolynomialDerivativeMatrix(x)
    for j=1:N+1
        for i=1:N+1
            DG.D̂[i,j] = -D[j,i]*w[j]/w[i]
        end
    end
end

function  DGDerivative(DG::NodalDiscontinuousGalerkin,ΦL,ΦR,Φ)
"""
       DGDerivative(DG,ΦL,ΦR,Φ)
       DG::Discontinuous Galerkin class
       ΦL::boundary values left
       ΦR::boundary values right
       Φ::Approximated function
       Computes the first spatial derivative via the DG approximation
       Algorithm 60 of Kopriva's book
    """

    Φ1 = MxVDerivative(DG.D̂,Φ)
    N=DG.N
    for j=1:N+1
        Φ1[j] = Φ1[j] + (ΦR*DG.lR[j] - ΦL*DG.lL[j])/DG.w[j]
    end
    return Φ1
end

function  InterpolateToBoundary(Φ,l)
"""
       InterpolateToBoundary(Φ,l)
       Φ::Nodal values
       l::interpolating polynomial
       Interpolates the nodal values to the boundary
       Part of Algorithm 61 of Kopriva's book
    """
    interpolatedValue = 0.0
    N=size(Φ,1)
    for j=1:N
        interpolatedValue = interpolatedValue + l[j] * Φ[j]
    end
    return interpolatedValue
end

function  DGTimeDerivative(DG::NodalDiscontinuousGalerkin,t,c,g,Δt,tn)
"""
       DGTimeDerivative(DG,t,c)
       DG::Discontinous Galerkin class
       t::time
       c::wave speed
       g::boundary value
       Determines the time derivative for the DG approximation
       Algorithm 61 of Kopriva's book
    """

    if (c>0)
        ΦL= sinpi(-t)
        ΦR = InterpolateToBoundary(DG.Φ,DG.lR)
    else
        ΦR=sinpi(2-t)
        ΦL = InterpolateToBoundary(DG.Φ,DG.lL)
    end
    Φt = -c * DGDerivative(DG,ΦL,ΦR,DG.Φ)
    return Φt
end

#TODO Replace interpolate to boundary with BLAS x dot for efficiency in the future
function  DGStepByRK3!(tn,Δt,DG::NodalDiscontinuousGalerkin,a,b,g,c,gb)
"""
       DGStepByRK3(tn,Δt,DG,a,b,g,c,gb)
       tn::current time
       Δt:: time step
       DG ::Discontinous Galerkin class
       a,b,g::RK3 tableau
       c::wave speed
       gb:: boundary values
       Advances the DG method in time using RK3
       Algorithm 62 of Kopriva's book
    """
    N=size(DG.Φ,1)
    G=zeros(Float64,N)
    for m=1:3
        t = tn + b[m] * Δt
        Φt = DGTimeDerivative(DG,t,c,gb,Δt,tn)
        for j=1:N
            G[j] = a[m] * G[j] + Φt[j]
            DG.Φ[j] = DG.Φ[j] + g[m] * Δt*G[j]
        end
    end
end

function  DGDriver!(N,NT,Nout,T,DG::NodalDiscontinuousGalerkin,a,b,g,c,gb)
"""
       DGDriver(N,NT,Nout,T,DG,a,b,g,c,gb)
       N::Order of approximation
       NT::Number of time steps
       Nout::Number of output points
       DG::DG class
       a,b,g::RK3 tableau
       c::scalar transport coefficient
       gb::boundary values
    """

    Δt = T/NT
    tn=0.0
    for n=0:NT-1
        DGStepByRK3!(tn,Δt,DG,a,b,g,c,gb)
        tn =(n+1)*Δt
    end
    X=zeros(Float64,Nout+1)
    for j=0:Nout
        X[j+1]=-1+2*j/Nout
    end
    w=BarycentricWeights(x)
    T=zeros(Float64,Nout+1,N)
    T=PolynomialInterpolationMatrix(x,w,X)
    ΦI = zeros(Float64,Nout+1)
    ΦI = InterpolateToNewPoints(T,DG.Φ)
    return ΦI
end


"""
    Mass(N)
    N::Order of approximation
    Q::Number of quadrature points

    Algorithm 5.2 Giraldo's book
"""

abstract type Abstract_Mass_Matrix end
struct MassMatrix1D <: Abstract_Mass_Matrix end
struct MassMatrix2D <: Abstract_Mass_Matrix end
struct MassMatrix3D <: Abstract_Mass_Matrix end

function ElementMassMatrix_1D(ψ, ω, N, Q, TFloat)
    
    M = zeros(Float64,N+1,N+1)   
    for k=1:Q+1
        for i=1:N+1
            for j=1:N+1
                M[i,j] = M[i,j] + ω[k]*ψ[i,k]*ψ[j,k]
            end
        end
    end
    
    return M

end

function ElementMassMatrix(N, Q, PT::MassMatrix1D, TFloat)
    
    ψ = zeros(TFloat, N+1, Q+1)
    lgl = St_lgl{TFloat}(zeros(TFloat, N+1),
                         zeros(TFloat, N+1))
    
    build_Integration_points!(lgl, N)
    ξ = lgl.ξ
    ω = BarycentricWeights(ξ)

    for k = 1:Q+1
        L = LagrangeInterpolatingPolynomials(ξ[k],ξ,ω)
        for i = 1:N+1
            ψ[i,k] = L[i]
        end
    end
    
    M = ElementMassMatrix_1D(ψ, ω, N, Q, TFloat)
    return M
end


