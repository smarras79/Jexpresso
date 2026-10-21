"""
Standard Henyey-Greenstein phase function, normalized so ∫Φ dΩ = 4π
(matching Φ_isotropic = 1/(4π) as the g=0 limit).
"""
function phase_function_hg(θ::Float64, θ1::Float64, ϕ::Float64, ϕ1::Float64, g::Float64)
    cosΘ = cos(θ)*cos(θ1) + sin(θ)*sin(θ1)*cos(ϕ-ϕ1)
    return (1 - g^2) / (4π * (1 + g^2 - 2*g*cosΘ)^1.5)
end

function phase_function_mie_like(θ::Float64, θ1::Float64, ϕ::Float64, ϕ1::Float64, ripple_N::Float64;
                                  g1=0.85, g2=-0.2, f=0.9, ripple_ε=0.25)
    cosΘ = cos(θ)*cos(θ1) + sin(θ)*sin(θ1)*cos(ϕ-ϕ1)
    HG(c,g) = (1-g^2) / (4π*(1+g^2-2*g*c)^1.5)
    envelope = f*HG(cosΘ,g1) + (1-f)*HG(cosΘ,g2)
    Θ = acos(clamp(cosΘ, -1.0, 1.0))
    ripple = 1 + ripple_ε*cos(ripple_N*Θ)
    return envelope * ripple
end
