# Advanced Function Library
# Created 23.10.17

# Imports
from __future__ import division
import numpy as np
from numba import njit

# Optimal Average
def optavg( x, sig ):
	wt = 1. / sig**2
	top = np.nansum( x * wt )
	bot = np.nansum( wt )
	
	avg = top / bot
	var = 1. / bot
	std = np.sqrt( var )

	return avg, std

# Optimal Scaling
def optscl( x, sig, pfl ):
	wt = pfl / sig**2
	top = np.nansum( x * wt )
	bot = np.nansum( pfl * wt )

	scl = top / bot
	var = 1. / bot
	std = np.sqrt( var )

	return scl, std

# Fit Line y = A + B * ( x - x0 )
def fitline( x, y, sig ):

	# centroids and scales
	x0 = optavg( x, sig )[0]
	y0 = optavg( y, sig )[0]

	xmin, xmax = np.nanmin( x ), np.nanmax( x )
	ymin, ymax = np.nanmin( y ), np.nanmax( y )
	xscale = abs( xmax - xmin )
	yscale = abs( ymax - ymin )

	# Hessian Matrix
	wt = ( yscale / sig )**2
	xh = ( x - x0 ) / xscale
	yh = ( y - y0 ) / yscale 

	H11 = np.nansum( wt )
	H12 = - np.nansum( xh * wt )
	H22 = np.nansum( xh**2 * wt )

	C1 = np.nansum( xh * yh * wt )
	C2 = np.nansum( yh * wt )

	det = H11 * H22 - H12**2
	A = ( C2 * H22 + H12 * C1 ) / det
	B = ( H11 * C1 + H12 * C2 ) / det

	'''
	H = [ [ H11, H12 ], [ H12, H22 ] ]
	C = [ C1, C2 ]
	Hinv = np.linalg.inv( H )

	A, B = Hinv.dot( C )
	'''

	Asig = 1. / np.sqrt( H11 )
	Bsig = 1. / np.sqrt( H22 )
	
	# de-scale
	A = A * yscale + y0
	Asig = Asig * yscale
	Bscale = yscale / xscale
	B = B * Bscale
	Bsig = Bsig * Bscale

	return A, B, Asig, Bsig, x0


