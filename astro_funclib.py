# Astronomical Toolkit Function Library
# Created 23.10.17

# Import modules
import numpy as np
import cosmolopy.distance as cd
from numba import njit

# Cosmology
cosmo = { 'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72 }
cosmo = cd.set_omega_k_0( cosmo )


# Degrees to Hours
@njit
def deg2hr( deg ):
	hr = deg * ( 24. / 360. )
	return hr

# Accretion Disk Color
@njit
def accretion_color( wav1, wav2 ):
	return - 2.5 * ( - 1. / 3. ) * np.log10( wav1 / wav2 )