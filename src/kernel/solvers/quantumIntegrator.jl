include("quantumIntegrator/BldTPoly.jl")
include("quantumIntegrator/Calc_FlowVarResults.jl")
include("quantumIntegrator/Calc_Noz_Area.jl")
include("quantumIntegrator/CalcBCmSW.jl")
include("quantumIntegrator/CalcBCpSW.jl")
include("quantumIntegrator/InitParms.jl")
include("quantumIntegrator/IntegrateGij.jl")
include("quantumIntegrator/IntegrateODE.jl")
include("quantumIntegrator/IPrtn.jl")
include("quantumIntegrator/Set_XGrid.jl")
include("quantumIntegrator/Calc_ExactResultsmSW.jl")
include("quantumIntegrator/SetInCond.jl")

function quantumIntegrator(u, params, inputs)
    N, delta1, n, k = InitParms(16, 0.005, 0.005, 1400) #all from qns_inputdata in jqc (user inputs)
    a = 0 
    d = 3#params.mesh.neqs
    r = 2
    n = 8#Int(n)
    N = 64#Int(N)
    @info n, N
    k = Int(k)
    Gamma = 1.4
    Tot_Int_Pts = params.mesh.npoin - 2
    Tot_X_Pts = params.mesh.npoin
    Shock_Flag = 0
    Exit_Pressure = 1
    ithroat = (Tot_X_Pts + 1)/2
    rho = 1
    x, Del_x = Set_XGrid(params.xmin, params.xmax, Tot_X_Pts)
    A = Calc_Noz_Area(x)
    In_Mass_Flow = 0.579

    Delta_t = inputs[:Δt]
    # Mrho_E = zeros(Float64, 1, Tot_X_Pts) #set to 0 for now?
    # Temp_E =  zeros(Float64, 1, Tot_X_Pts)#?
    # Mach_E =  zeros(Float64, 1, Tot_X_Pts)#?
    # Press_E =  zeros(Float64, 1, Tot_X_Pts)#?
    # Vel_E =  zeros(Float64, 1, Tot_X_Pts)#?
    Mach_E, Mrho_E, Press_E, Temp_E, Vel_E = Calc_ExactResultsmSW(Gamma, Tot_X_Pts, A) #tentative
    InitVal = zeros(Float64, d, Tot_X_Pts) # rhs[ip, ieq] corresponds to InitVal[ieq, ip]...
    for i=1:d
        for j=1:Tot_X_Pts
            InitVal[i, j] = u[params.mesh.npoin*(i-1)+j]
        end
    end
    U2_in = zeros(Float64, Tot_X_Pts, n+1) #initial mass flow rate in col. 1, all other cols are 0
    U2_in[1, :] .= In_Mass_Flow
    #InitVal, Delta_t, In_Mass_Flow_Noisy = SetInCond(shock_flag, In_Mass_Flow, gamma, x, Del_x, A, d, Mrho_E, Temp_E, ICMFlowErrScale, ICrhoErrScale, ICtempErrScale)
    b =  0.0001#1400*Delta_t#inputs[:Δt]
    t, hbar = IPrtn(a, b, n, N)
    ff0_throat_in = zeros(Float64, d, n)
    ff1_throat_in = zeros(Float64, d, n)
    ff2_throat_in = zeros(Float64, d, n)
    #mach_E = mach number at all pts
    #mrho_E = mass density at all grid pts
    # press_E = pressure at all grid-pts
    # temp_E = temperature at all grid-pts
    # vel_E = velocity at all grid-pts
    U2, Mach_D, Mrho_D, Press_D, Temp_D, Vel_D, Rel_MachErr, Rel_MrhoErr, 
    Rel_PressErr, Rel_TempErr, Rel_VelErr, #=AvRelTempErr, AvPlusSDevRelTempErr, 
    AvMinusSDevRelTempErr, AvRelMachErr, AvPlusSDevRelMachErr, AvMinusSDevRelMachErr, 
    AvRelMrhoErr, AvPlusSDevRelMrhoErr, AvMinusSDevRelMrhoErr, AvRelPressErr, 
    AvPlusSDevRelPressErr, AvMinusSDevRelPressErr, AvU2,=# ff0_throat, ff1_throat, 
    ff2_throat, FinalVal = IntegrateODE(d, n, N, hbar, r, Del_x, Gamma, Tot_Int_Pts, k, 
    Tot_X_Pts, Shock_Flag, Exit_Pressure, ithroat, a, delta1, rho, InitVal#=?=#, A, 
    t, U2_in, ff0_throat_in, ff1_throat_in, ff2_throat_in, Mach_E, Mrho_E, Press_E, 
    Temp_E, Vel_E, In_Mass_Flow, params)
    #finalU = zeros(Float64, d*Tot_X_Pts)
    for ieq=1:d
        for ip=1:Tot_X_Pts
            u[(ieq-1)*Tot_X_Pts + ip] = FinalVal[ieq, ip]
        end
    end
    return u
end