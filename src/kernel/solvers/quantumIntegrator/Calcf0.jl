include("CalcFlux.jl"); 
include("CalcSource.jl"); 
include("CalcFunc.jl");

function Calcf0( d, Tot_X_Pts, Tot_Int_Pts, Gamma, Del_x, U, A )
    #CALCF0 evaluates ODE driver function f(U) at flow variables U 
    #   Calcf0 evaluates ODE driver function f(U) at the primary flow variable
    #           U at start of subsubinterval j at each interior grid-point.
    #           Flow is 1D inviscid; compressible flow through a nozzle.
    #   
    #   INPUTS:
    #           d = number of components of ODE driver function f(U)
    #           Tot_X_Pts = number of grid-points
    #           Tot_Int_Pts = number of interior grid-points
    #           Gamma = ratio of specific heats
    #           Del_x = distance between grid-points
    #           U = d x Tot_X_Pts array containing primary flow variables at
    #                   start of subsubinterval j at ALL grid-points
    #           A = 1 x Tot_X_Pts array storing nozzle area at all grid-points
    #
    #   OUTPUT:
    #           f0 = d x Tot_Int_Pts array storing ODE driver function f(U)
    #                   evaluated at start of subsubinterval j at each interior
    #                   grid-point
    #
    #   Support functions: CalcFlux; CalcSource; CalcFunc
    
    #  evaluate flow flux F = d x Tot_X_Pts array
    
    ######evaluate the flux term from eqn 2?######
      F = CalcFlux(U, Gamma, d, Tot_X_Pts)
       
    #  evaluate source current J for momentum equation of motion which is()
    #   1 x Tot_Int_Pts array
    
      ###### evaluatue source matrix from eqn 2 ######
       J = CalcSource(U, A, Gamma, Tot_X_Pts)
       
    #  evaluate f0 which stores ODE driver function values at all interior
    #   grid-points at start of subsubinterval j. f0 is d x Tot_Int_Pts array
    
    ###### evaluate driver function values at grid pts ######
       f0 = CalcFunc(F, J, Del_x, d, Tot_Int_Pts)
   
   return f0;
    end
