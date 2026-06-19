function user_inputs()

    # ENV-driven knobs (defaults = DHIT DNS reference: Re_λ=162, constant μ).
    _mu    = parse(Float64, get(ENV, "JEXP_MU",    "3.64e-3"))   # Re_λ=162 molecular viscosity (DNS); SGS intensity (LES)
    _npv   = get(ENV, "JEXP_NPV",  "true")  == "true"
    _ent   = get(ENV, "JEXP_ENT",  "false") == "true"
    _lkep  = get(ENV, "JEXP_LKEP", "true")  == "true"
    _visc  = get(ENV, "JEXP_VISC", "AV")                          # "AV" const μ (DNS) / "SMAG" (LES)
    _flux  = get(ENV, "JEXP_FLUX", "ranocha")
    _mumol = parse(Float64, get(ENV, "JEXP_MUMOL", "3.64e-3"))    # molecular μ floor for LES (Re_λ=162)
    _cs    = parse(Float64, get(ENV, "JEXP_CS", "0.17"))          # Smagorinsky constant (Gassner-style)
    _tend  = parse(Float64, get(ENV, "JEXP_TEND", "5.0"))
    _mesh  = get(ENV, "JEXP_MESH", "./meshes/gmsh_grids/cube_periodic_6.msh")  # Gassner coarse LES = 6 cells, N=7
    _outdir = get(ENV, "JEXP_OUTDIR", "./output-dhit/")

    _visc_model  = _visc == "SMAG" ? SMAG() : AV()
    _volume_flux = _flux == "kennedy_gruber" ? kennedy_gruber() :
                   _flux == "shima"          ? shima()          :
                   _flux == "chandrashekar"  ? chandrashekar()  : ranocha()

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 2.0e-3,
        :tinit                => 0.0,
        :tend                 => _tend,
        :diagnostics_at_times => (0.0:0.25:_tend),
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 7,                    # Gassner polynomial degree N=7
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [0.0, _mu, _mu, _mu, _mu],
        :visc_model           => _visc_model,
        :mu_molecular         => _mumol,
        :C_smag               => _cs,
        :Pr                   => 0.71,
        :energy_equation      => "energy",
        #---------------------------------------------------------------------------
        :lkep                    => _lkep,
        :new_primitive_variables => _npv,
        :entropy_variables       => _ent,
        :volume_flux             => _volume_flux,
        #---------------------------------------------------------------------------
        # Rogallo IC parameters
        #---------------------------------------------------------------------------
        :rogallo_N            => 64,
        :rogallo_kp           => 2.0,
        :rogallo_seed         => 1,
        #---------------------------------------------------------------------------
        # DHIT kinetic-energy spectra (Flad & Gassner 2017): the DGSEM field is
        # super-sampled onto an equidistant N_uni³ grid with N_uni = nel·(N+1)
        # (= "DOF per direction"; 6 cells × N=7 → 48³, 18 cells × N=7 → 144³),
        # FFT'd, and shell-binned. Writes dhit_spectrum_{1D,3D}_<iout>.dat plus
        # the dhit_integrals.dat time series at every :diagnostics_at_times.
        # The 3D file also carries the dissipation-rate spectrum D(k)=2νk²E(k),
        # and dhit_integrals.dat the Fig. 9 energy/dissipation integrals up to
        # :dhit_kcut, using the fixed molecular ν = :dhit_nu (default mu_molecular).
        #---------------------------------------------------------------------------
        :dhit_spectra         => true,
        :dhit_kcut            => 16.0,        # Fourier-space integration cutoff (Fig. 9)
        :dhit_nu              => _mumol,      # fixed kinematic ν for ε=2ν∫k²E dk (ρ_ref=1)
        #---------------------------------------------------------------------------
        :lread_gmsh          => true,
        :gmsh_filename       => _mesh,
        #---------------------------------------------------------------------------
        # Light erf filter (KEP needs little; keep small to avoid contaminating E(k))
        #---------------------------------------------------------------------------
        # Split-form / KEP is the de-aliasing → dissipation-free baseline (no filter), as in Gassner.
        :lfilter             => false,
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => false,
        :lwrite_initial      => true,
        :output_dir          => _outdir,
        :loutput_pert        => false,
        :linitial_refine     => false,
        :ladapt              => false,
    ) #Dict

    return inputs
end
