"""
    ERA5 LES - Source Terms

    Includes gravity and optional ERA5 nudging/forcing.

    Author: Jexpresso Development Team
    Date: 2025-11-30
"""

function user_source!(S::SubArray,
                     q::AbstractArray,
                     qe::AbstractArray,
                     npoin::Int64,
                     ::CL, ::PERT, ::CompEuler,
                     PhysConst;
                     neqs=7)

    g = PhysConst.g

    # Gravity source term
    ρ = q[1] + qe[1]
    S[4] = -ρ * g  # Vertical momentum

    # Additional source terms can be added here
    # For example, ERA5 nudging if enabled
end

function user_source!(S::SubArray,
                     q::AbstractArray,
                     qe::AbstractArray,
                     npoin::Int64,
                     ::CL, ::TOTAL, ::CompEuler,
                     PhysConst;
                     neqs=7)

    g = PhysConst.g

    # Gravity source term
    ρ = q[1]
    S[4] = -ρ * g  # Vertical momentum
end
