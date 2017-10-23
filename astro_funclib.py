# Astronomical Toolkit Function Library
# Created 23.10.17

# Import modules
import numpy as np
import cosmolopy.distance as cd

# Cosmology
cosmo = { 'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72 }
cosmo = cd.set_omega_k_0( cosmo )
