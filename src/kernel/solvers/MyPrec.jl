module MyPrecClass

export MyPrec, updateA!, sol!

using AlgebraicMultigrid: solve
using IncompleteLU

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

function Jacobisol!(P::MyPrec, b::AbstractVector, y::AbstractVector)
    precision = P.kwargs[:precision]

    b_1 = precision.(b)

    y .= P.prec \ b_1
end

function ilusol!(P::MyPrec, b::AbstractVector, y::AbstractVector)
    precision = P.kwargs[:precision]

    b_1 = precision.(b)

    IncompleteLU.forward_substitution!(P.prec, b_1)
    IncompleteLU.backward_substitution!(P.prec, b_1)
    y .= b_1
end

end
