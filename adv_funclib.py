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
	
	if bot == 0:
		avg = np.nan
		std = np.nan

	else:
		avg = top / bot
		var = 1. / bot
		std = np.sqrt( var )

	return avg, std

# Optimal Scaling
def optscl( x, sig, pfl ):
	wt = pfl / sig**2
	top = np.nansum( x * wt )
	bot = np.nansum( pfl * wt )

	if bot == 0:
		scl = np.nan
		std = np.nan

	else:
		scl = top / bot
		var = 1. / bot
		std = np.sqrt( var )

	return scl, std

# Standard Deviation from model
def std( x, modx, sig ):
	dlen = float( len( x ) )
	var = np.nansum( ( x - modx )**2 / dlen )
	return np.sqrt( var )

# Chisq
def chisq( x, modx, sig ):
	return np.nansum( ( ( x - modx ) / sig )**2 )

# Linear Interpolation
# Does not provide extrapolated sigma
def linint( x, y, ysig, x0, extrapolation = False ):
	yvar = ysig**2
	if x0 < min( x ):
		if extrapolation:
			scl  = ( y[1] - y[0] ) / ( x[1] - x[0] )
			b = y[0] - x[0] * scl
			val =  scl * x0 + b

			valvar = 0.

		else:
			val = np.nan
			valvar = np.nan


	elif x0 > max( x ):
		if extrapolation:
			scl = ( y[-1] - y[-2] ) / ( x[-1] - x[-2] )
			b = y[-1] - x[-1] * scl
			val = scl * x0 + b

			valvar = 0.

		else:
			val = np.nan
			valvar = np.nan


	elif x0 in x:
		val = y[ x==x0 ][0]
		valvar = yvar[ x==x0 ][0]


	else:
		xlo = max( x[ x < x0 ] )
		xhi = min( x[ x > x0 ] )

		ylo = y[ x == xlo ][0]
		yhi = y[ x == xhi ][0]
		
		yvarlo = yvar[ x == xlo ][0]
		yvarhi = yvar[ x == xhi ][0]

		scl = ( x0 - xlo ) / ( xhi - xlo )

		val = ylo + scl * ( yhi - ylo )

		valvar = scl * yvarhi + ( 1 - scl ) * yvarlo

	return val, np.sqrt(valvar)


# Fit Line y = A + B * ( x - x0 )
def fitline( x, y, sig, orthogonalize = True ):

	# centroids and scales
	x0 = optavg( x, sig )[0]
	y0 = optavg( y, sig )[0]

	xmin, xmax = np.nanmin( x ), np.nanmax( x )
	ymin, ymax = np.nanmin( y ), np.nanmax( y )
	xscale = abs( xmax - xmin )
	yscale = abs( ymax - ymin )

	# Hessian Matrix
	if orthogonalize:
		wt = ( yscale / sig )**2
		xh = ( x - x0 ) / xscale
		yh = ( y - y0 ) / yscale 
	else:
		wt = 1. / sig**2
		xh = x
		yh = y

	H11 = np.nansum( wt )
	H12 = - np.nansum( xh * wt )
	H22 = np.nansum( xh**2 * wt )

	C1 = np.nansum( xh * yh * wt )
	C2 = np.nansum( yh * wt )

	det = H11 * H22 - H12**2
	A = ( C2 * H22 + H12 * C1 ) / det
	B = ( H11 * C1 + H12 * C2 ) / det

	Asig = 1. / np.sqrt( H11 )
	Bsig = 1. / np.sqrt( H22 )
	
	if H12 == 0.:
		ABcov = 0.
	else:
		ABcov = 1 / H12

	p1 = np.sqrt( H11 / det )
	p2 = np.sqrt( H22 / det )
	corr =  ( H12 / det ) / ( p1 * p2 )

	if orthogonalize:
		# de-scale
		A = A * yscale + y0
		Asig = Asig * yscale
		Bscale = yscale / xscale
		B = B * Bscale
		Bsig = Bsig * Bscale

	return A, B, Asig, Bsig, ABcov, corr, x0


