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

# Accretion Disk Spectrum in Magnitudes
# Cannot Njit
def accretion_magspec( wav, fidwav, fidmag = 0. ):
	mag = -2.5 * np.log10( ( wav / fidwav )**(-1./3.) ) + fidmag
	return mag

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
	wavz = np.array( [ w / ( 1. + z ) for w in wav ] )
	return wavz

# Rest to observed wavelength
@njit
def rest2obswav( wavz, z ):
	wav = np.array( [ w * ( 1. + z) for w in wavz ] )
	return wav

# Apparent to absolute magnitude
def mag2amag( mag, magsig, z ):
	CoMDist = cd.comoving_distance( z, **cosmo )
	cterm = - 5 * np.log10( CoMDist ) - 25.
	return mag + cterm, magsig

# Absolute to apparent magnitude
def amag2mag( amag, amagsig, z ):
	CoMDist = cd.comoving_distance( z, **cosmo )
	cterm = - 5 * np.log10( CoMDist ) - 25.
	return amag - cterm, amag_sig

# Extinction correction
# Cannot njit
def extcorr( intable, extmag_ref, extmag_coeff, verbose = False, filters = None ):
	if verbose:
		print 'Extcorr: Corrected extinction'

	table = intable.copy()
	for i, colname in enumerate( extmag_coeff.keys() ):
		col = table[colname]
		coeff = extmag_coeff[colname]
		extmag = coeff * extmag_ref
		table[colname] = col + extmag

		if filters is None:
			pass
		elif verbose:
		    print '* %s ... %2.2f mag' %( filters[i], extmag )

	return table

# Chemical line plotter
def chemplot( ax, wavdict, xlo, xhi, ylo, yhi, color = 'k', tcolor = 'k' ):
	name = np.array( wavdict.keys() )
	wav = np.array( wavdict.values() )
	print wav
	mask = ( ( wav > 0.98 * xlo ) & ( wav < 0.98 * xhi) )
	print mask
	ax.vlines( wav[mask], ylo, yhi, colors = color )
	yrange = yhi - ylo
	for i in xrange( len( wav[mask] ) ):
		pwav = wav[mask][i]
		pname = name[mask][i]
		ax.text( pwav, yhi +  0.1 * yrange, pname, color = tcolor )


	
