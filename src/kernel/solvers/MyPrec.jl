module MyPrecClass

export MyPrec, updateA!, sol!

using LinearSolve: solve

mutable struct MyPrec
    A::AbstractMatrix
    prec::Any
    kwargs::Dict
end

function updateA!(P::MyPrec, A)
    P.A = A
end

function sol!(P::MyPrec, b::AbstractVector, y::AbstractVector)
    maxiter   = P.kwargs[:maxiter]
    abstol    = P.kwargs[:abstol]
    precision = P.kwargs[:precision]

    b_1 = precision.(b)

    sol = solve(P.A, b_1, P.prec, maxiter=maxiter, abstol=abstol)

    y .= sol
end

end
