function perfectGasLaw_θPtoρ(PhysConst::PhysicalConst; θ=300.0, Press=100000.0)
    T = typeof(θ)
    return (T(1.0)/θ)*(Press/PhysConst.C0)^(T(1.0)/PhysConst.γ) #ρ
end


function perfectGasLaw_ρTtoP(PhysConst::PhysicalConst; ρ=1.25, Temp=300.0)
    
    return ρ*PhysConst.Rair*Temp #Press
end

function perfectGasLaw_ρPtoT(PhysConst::PhysicalConst; ρ=1.25, Press=100000.0)
    return  Press/(ρ*PhysConst.Rair) #Temp
end


function perfectGasLaw_TPtoρ(PhysConst::PhysicalConst; Temp=300.0, Press=100000.0)
    return  Press/(Temp*PhysConst.Rair) #ρ
end

function perfectGasLaw_ρθtoP(PhysConst::PhysicalConst; ρ=1.25, θ=300.0)
    return PhysConst.C0*(ρ*θ)^PhysConst.γ #Press
end

function perfectGasLaw_ρθtoP!(Press::Float64, PhysConst::PhysicalConst; ρ=1.25, θ=300.0)
    Press = PhysConst.C0*(ρ*θ)^PhysConst.γ #Press
end

function perfectGasLaw_ρθtoP(PhysConst::PhysicalConst, ρ::AbstractArray, θ::AbstractArray)
    return PhysConst.C0 .* (ρ .* θ) .^ PhysConst.γ
end

function perfectGasLaw_ρθtoP!(Press::Array{Float64}, PhysConst::PhysicalConst; ρ=1.25, θ=300.0)
    Press[1] = PhysConst.C0*(ρ*θ)^PhysConst.γ #Press
end

function perfectGasLaw_θPtoT!(T::Float64, PhysConst::PhysicalConst; θ=300.0, Press=100000.0)
    T = θ * (Press / PhysConst.pref)^(PhysConst.Rair / PhysConst.cp)
end

function perfectGasLaw_θPtoT(PhysConst::PhysicalConst; θ=300.0, Press=100000.0)
    return θ * (Press / PhysConst.pref)^(PhysConst.Rair / PhysConst.cp)
end

# Vectorized version
function perfectGasLaw_θPtoT(PhysConst::PhysicalConst, θ::Vector{Float64}, Press::Vector{Float64})
    return θ .* (Press ./ PhysConst.pref).^(PhysConst.Rair / PhysConst.cp)
end

function perfectGasLaw_ρθqtoP(ρ, θ, qt, ql)
    L_v = 2.5008e6
    kappad = 0.28571428571
    cpv = 1859.0
    cpl = 4181.0
    p0 = 101325.0
    R_d = 287.0024
    molmass_ratio = 1.608079364

    cpd = R_d / kappad
    qd = 1.0 - qt
    qv = qt - ql
    Rm = R_d * (1 + (molmass_ratio - 1) * qt - molmass_ratio * ql)
    cpm = cpd * qd + cpv * qv + cpl * ql
    press = (ρ * R_d * θ)^(cpm / (cpm - Rm)) * (p0)^(-Rm / (cpm - Rm))

    return press
end


function perfectGasLaw_θqPtoT(θ, qt, ql, p)
    L_v = 2.5008e6
    kappad = 0.28571428571
    cpv = 1859.0
    cpl = 4181.0
    p0 = 101325.0
    R_d = 287.0024
    molmass_ratio = 1.608079364

    cpd = R_d / kappad
    qd = 1.0 - qt
    qv = qt - ql
    Rm = R_d * (1 + (molmass_ratio - 1) * qt - molmass_ratio * (qt - qv))
    cpm = cpd * qd + cpv * qv + cpl * ql
    theta = R_d / Rm * θ
    temperature = theta * (p / p0) ^ (Rm / cpm)

    return temperature
end

function perfectGasLaw_ρPtoθ(PhysConst::PhysicalConst; ρ=1.25, Press=100000.0)
    T = typeof(ρ)    
    return (T(1.0)/ρ)*(Press/PhysConst.C0)^(T(1.0)/PhysConst.γ) #θ
    
end

function moistPressure(PhysConst::PhysicalConst; ρ=1.25, Tv=300.0, qv = 0.0)
    T = typeof(Tv)
    #@info ρ*Temp*PhysConst.Rair + ρ*Temp*qv*PhysConst.Rvap, Temp, qv, ρ, PhysConst.Rvap
    return (T(ρ*Tv*PhysConst.Rair))
end


"""
Calculate pressure using hydrostatic equation with iterative approach.
dp/dz = -g * rho = -g * p / (R_d * T_v)
"""
function perfectGasLaw_θqvtoP(PhysConst::PhysicalConst,
                              height::Vector{Float64}, 
                              theta::Vector{Float64}, 
                              qv::Vector{Float64}; p0=101325.0)
    
    n = length(height)
    pressure = zeros(Float64, n)
    
    # Surface pressure assumption (can be adjusted)
    pressure[1] = p0  # Standard sea level pressure [Pa]
    
    # Integrate hydrostatic equation upward
    for i in 2:n
        dz = height[i] - height[i-1]
        
        # Use average values for integration
        p_avg = pressure[i-1]  # Initial guess
        
        # Iterative solution for pressure
        for iteration in 1:5
            T_avg   = perfectGasLaw_θPtoT(PhysConst, (theta[i] + theta[i-1]) / 2, p_avg)
            T_v_avg = calculate_virtual_temperature(PhysConst, T_avg, (qv[i] + qv[i-1]) / 2)
            
            # Hydrostatic equation: dp = -g * p / (R_d * T_v) * dz
            dp          = -PhysConst.g * p_avg / (PhysConst.Rair * T_v_avg) * dz
            pressure[i] = pressure[i-1] + dp
            p_avg       = (pressure[i] + pressure[i-1]) / 2
        end
    end
    
    return pressure
end

# Function to update p_ref_theta
function create_updated_TD_Parameters(new_p_ref_theta::FT) where {FT}
    inner_dict = Dict{String, Any}(
                                    "value" => new_p_ref_theta,
                                    "type" => "float", 
                                    "description" => "Reference pressure used in potential temperature definition"
                                  )
    override_file = Dict("potential_temperature_reference_pressure" => inner_dict)
    toml_dict = CP.create_toml_dict(FT;override_file = override_file)
    return TP.ThermodynamicsParameters(toml_dict)
end

"""
Calculate virtual temperature accounting for water vapor.
T_v = T * (1 + qv/epsilon) / (1 + qv/1000)
epsilon = 0.6217504
"""
function calculate_virtual_temperature(PhysConst::PhysicalConst, T::Float64, qv::Float64)
    qv_kg = qv / 1000.0  # Convert g/kg to kg/kg
    return T * (1 + qv_kg / 0.6217504) / (1 + qv_kg)
end

# Vectorized version
function calculate_virtual_temperature(PhysConst::PhysicalConst, T::Vector{Float64}, qv::Vector{Float64})
    qv_kg = qv ./ 1000.0  # Convert g/kg to kg/kg
    return T .* (1 .+ qv_kg ./ 0.6217504) ./ (1 .+ qv_kg)
end
