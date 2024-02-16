function Calc_dFdt(U, ff_vals, ff_Bvals, Gamma, d, Tot_X_Pts)
  #CALC_DFDT evaluates first time derivative of flow fluxes at grid-points
  #   Calc_dFdt evaluates first time derivative of flow fluxes at
  #   all grid-points for 1D compressible; inviscid flow through a nozzle.
  #
  #   INPUTS:
  #           U = d x Tot_X_Pts array with primary flow variables at all()
  #                   grid-points
  #           ff_vals = d x Tot_Int_Pts array with value of ODE driver
  #                       function at all interior grid-points
  #           ff_Bvals = d x 2 array with values of dU/dt at boundary points
  #           Gamma = ratio of specific heats
  #           d = number of primary flow variables at a grid-point
  #           Tot_X_Pts = number of grid-points
  #
  #   OUTPUT:
  #           dFdt = d x Tot_X_Pts array storing first time derivative of
  #                       flow fluxes at all grid-points

  dFdt = zeros(d, Tot_X_Pts)

  TotXPtm1 = Tot_X_Pts - 1

  fac4 = Gamma - 1
  fac5 = Gamma / 2
  fac6 = fac4 * fac5
  fac7 = (3 - Gamma) / 2
  fac8 = fac4 / Gamma

  # evaluate dFdt at boundaries:

  for ll = 1:2
    if (ll == 1)# at nozzle entrance

      fac1 = U[2, 1] / U[1, 1]
      fac2 = U[3, 1] / U[1, 1]
      fac3 = fac1 * fac2
      fac9 = (fac1)^(2)

      dFdt[1, ll] = ff_Bvals[2, ll]

      dFdt[2, ll] = fac7 * fac1 * (2 * ff_Bvals[2, ll]
                                   -
                                   fac1 * ff_Bvals[1, ll]) + fac8 * ff_Bvals[3, ll]

      dFdt[3, ll] = Gamma * (fac1 * ff_Bvals[3, ll]
                             +
                             fac2 * ff_Bvals[2, ll]
                             -
                             fac3 * ff_Bvals[1, ll])
      -fac6 * fac9 * (3 * ff_Bvals[2, ll]
                      -
                      2 * fac1 * ff_Bvals[1, ll])
    elseif (ll == 2)        # at nozzle exit()
      fac1 = U[2, Tot_X_Pts] / U[1, Tot_X_Pts]
      fac2 = U[3, Tot_X_Pts] / U[1, Tot_X_Pts]
      fac3 = fac1 * fac2
      fac9 = (fac1)^(2)

      dFdt[1, Tot_X_Pts] = ff_Bvals[2, ll]

      dFdt[2, Tot_X_Pts] = fac7 * fac1 * (2 * ff_Bvals[2, ll]
                                          -
                                          fac1 * ff_Bvals[1, ll]) + fac8 * ff_Bvals[3, ll]

      dFdt[3, Tot_X_Pts] = Gamma * (fac1 * ff_Bvals[3, ll]
                                    +
                                    fac2 * ff_Bvals[2, ll]
                                    -
                                    fac3 * ff_Bvals[1, ll])
      -fac6 * fac9 * (3 * ff_Bvals[2, ll]
                      -
                      2 * fac1 * ff_Bvals[1, ll])
    else
      disp(["Unknown switch case label: " ll])
    end
  end

  # evaluate dFdt at interior grid-points

  for ll = 2:TotXPtm1
    fac1 = U[2, ll] / U[1, ll]
    fac2 = U[3, ll] / U[1, ll]
    fac3 = fac1 * fac2
    fac9 = (fac1)^(2)

    IPLabel = ll - 1

    dFdt[1, ll] = ff_vals[2, IPLabel]

    dFdt[2, ll] = fac7 * fac1 * (2 * ff_vals[2, IPLabel]
                                 -
                                 fac1 * ff_vals[1, IPLabel])
    +fac8 * ff_vals[3, IPLabel]

    dFdt[3, ll] = Gamma * (fac1 * ff_vals[3, IPLabel]
                           +
                           fac2 * ff_vals[2, IPLabel]
                           -
                           fac3 * ff_vals[1, IPLabel])
    -fac6 * fac9 * (3 * ff_vals[2, IPLabel]
                    -
                    2 * fac1 * ff_vals[1, IPLabel])

  end

  return dFdt


end
