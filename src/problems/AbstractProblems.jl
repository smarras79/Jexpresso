abstract type AbstractProblem end

#1D
struct Wave1D <: AbstractProblem end
struct AD1D <: AbstractProblem end
struct Burgers1D <: AbstractProblem end

#2D
struct Adv2D <: AbstractProblem end
struct SW2D <: AbstractProblem end
