export Tables

"""
    Tables()

    This document compares the performance of Jexpresso against a legacy Fortran code for numerical weather prediction, focusing on wall clock times for a rising-thermal-bubble test.

        ---
        

        ### Wall Clock Time Comparison for Simulated 100 Seconds
        
        | Time integrator               | max Δt (s) | Effective resolution (m) | Order | Jexpresso (s) | F90 (s) |
        |-------------------------------|------------|--------------------------|-------|---------------|---------|
        | SSPRK33/RK33                  | 0.2        | $125\times 125$          | 4     | 9.75          | 9.2028  |
        | SSPRK53/RK35                  | 0.3        | $125\times 125$          | 4     | 9.00          | 10.53   |
        | SSPRK54                       | 0.4        | $125\times 125$          | 4     | 10.47         | -       |
        | DP5 (Dormand-Prince RK54)     | 0.6        | $125\times 125$          | 4     | 19.80         | -       |
        | SSPRK73                       | 0.4        | $125\times 125$          | 4     | 12.95         | -       |
        | SSPRK104                      | 0.6        | $125\times 125$          | 4     | 12.50         | -       |
        | CarpenterKennedy2N54          | 0.4        | $125\times 125$          | 4     | 10.57         | -       |
        | Tsit5                         | 2.0 (adaptive) | $125\times 125$      | 4     | 19.08         | -       |
        
        *Note: The wall clock times are to be taken with a ±0.2 due to a small variability from one simulation to the next one.*
        
        ### Comparison for t=1000 Seconds Without Diagnostics or VTK Writing
        
        | Time integrator               | max Δt (s) | Effective resolution (m) | Order | Jexpresso (s) | F90 (s) |
        |-------------------------------|------------|--------------------------|-------|---------------|---------|
        | SSPRK33/RK33                  | 0.2        | $125\times 125$          | 4     | 57.36         | 57.00   |
        
        ### Mass Conservation for Advective vs Flux Forms
        
        | Time integrator | Advection form                | Flux form                    |
        |-----------------|-------------------------------|------------------------------|
        | SSPRK33         | $7.622610626689869 \times 10^{-15}$ | $1.9545155453050947 \times 10^{-16}$ |
        | SSPRK33         | $5.081740417793246 \times 10^{-15}$ | $1.1727093271830568 \times 10^{-15}$ |
        | MSRK5           | $7.818062181220379 \times 10^{-16}$ | $3.9090310906101895 \times 10^{-16}$ |
        
        ### Exact vs Inexact Integration for Advection and Flux Forms
        
        | Time integrator | Integration Type | Advection form                | Flux form                      |
        |-----------------|------------------|-------------------------------|--------------------------------|
        | SSPRK33         | Exact            | $4.104482645140768 \times 10^{-15}$ | $4.495385754201793 \times 10^{-15}$ |
        | SSPRK33         | Inexact          | $5.081740417793246 \times 10^{-15}$ | $1.1727093271830568 \times 10^{-15}$ |
        
        
        ---

"""

function Tables()
end
