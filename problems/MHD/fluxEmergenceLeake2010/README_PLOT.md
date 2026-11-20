To create the plots, execute the following command from the Jexpresso home directory,
where NEW=output/MHD/fluxEmergenceLeake2010/merged

python3 tools/plot_leake2010.py dsgs $NEW --which energy --linear

Parameters:

- remove "--linear" if you want log scale (for log(rho) for example)
- "--which energy"    --> plots mu + paper's explicit viscosity (0.013 above the photosphere) for the mass and energy
- "--which momentum"  --> plots mu for rhou, rhov, rhow
- "--which induction" --> plots eta for the induction residuals only, plus the anomalous reistivity where active