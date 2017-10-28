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
	color = - 2.5 * ( - 1. / 3. ) * np.log10( wav1 / wav2 )
	return color

# Magnitude (AB) to Flux (mJy)
# Cannot Njit
def ab2mJy( m, msig ):
	f = 10. ** ( - 0.4 * ( m - 16.4 ) )
	fsig = msig * f * 0.4 * np.log( 10. )

	return f, fsig

# Flux (mJy) to Magnitude (AB)
# Cannot Njit
def mJy2ab( f, fsig ):
	m = 16.4 - 2.5 * np.log10( abs( f ) )
	msig = fsig / ( f * 0.4 * np.log( 10. ) )

	return m, msig

# Observed to rest wavelength
@njit
def obs2restwav( wav, z ):
	return wav / ( 1. + z )

# Rest to observed wavelength
@njit
def rest2obswav( wav, z ):
	return wav * ( 1. + z )

# Extinction correction
# Cannot njit
def extcorr( intable, extmag_ref, extmag_coeff, verbose = False ):
	if verbose:
		print 'Extcorr: Corrected extinction'

	table = intable.copy()
	for colname in extmag_coeff.keys():
		col = table[colname]
		coeff = extmag_coeff[colname]
		extmag = coeff * extmag_ref
		table[colname] = col + extmag

		if verbose:
		    print '* %s ... %2.2f mag' %( i, extmag )

	return table

