################################
#  Constants and IAXO params   #
################################

s2cm = 2.99792458e10 # cm per s
keV2m = 197.327053e-12 # m per keV^-1
T2keVsqr = 1.95353e-4 # keV^2 per T
c_p_conv = (T2keVsqr/keV2m)**2
radius_sol = 6.957e8 # m
distance_sol = 1.49597870700e11 # m
hbar = 6.582119514e-25 # GeV s
distance_factor = (radius_sol/keV2m)**3 / ( (1e2*distance_sol)**2 * (1e6*hbar) )

iaxo = { 'A': 2.26e4, # Area covered by IAXO in cm^2
         't': 3.0 * 365.0*24.0*60.0*60.0,  # Observation time in s
         'L': 20.0, # IAXO length in m
         'eff': 0.7*0.8, # Efficiency of detector and optics
         'B': 2.5 # Magnetic field in T
        }

# xi = kappa_s/T, mean and std
xi_m = 6.1
xi_std = 1.0