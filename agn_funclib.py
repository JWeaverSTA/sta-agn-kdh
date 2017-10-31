# AGN Project Function Library
# Mainly subroutine like functions
# Created 27.10.17

# Imports
from __future__ import division
import numpy as np
from astropy.table import Table


# Mask values given logical condition
def maskval( intable, val, test = 'eq', str_pattern = '_mag', verbose = False ):
	if test not in ( 'gt', 'gte', 'lt', 'lte', 'eq' ):
		print '"test" arguement "%s" not valid' %( test )
		return None

	if intable.masked:
		table = intable
	else:
		table = Table( intable, masked = True )
	tablelen = len( table )
	count = 0.
	mcount = np.zeros( tablelen )
	size = tablelen * len( table.colnames )
	for col in table.colnames:
		if str_pattern[1:] in col.split( str_pattern[0] ):
			if test == 'gt':
				mask = ( intable[col] > val )
			if test == 'gte':
				mask = ( intable[col] >= val )
			if test == 'lt':
				mask = ( intable[col] < val )
			if test == 'lte':
				mask = ( intable[col] <= val )
			if test == 'eq':
				mask = ( intable[col] == val )

			table[col].mask = np.logical_or( table[col].mask, mask )
			count += np.sum( mask )
			mcount[mask] += 1

	rcount = np.sum( mcount > 0 )

	if verbose:
		pc = count / size * 100.
		pce = rcount / tablelen * 100
		print 'Maskval: Masked values %s %3.2f with pattern "%s"' %(test, val, str_pattern)
		print '* %i of %i elements affected (%2.1f %%)' %( count, size, pc )
		print '* %i of %i epochs affected (%2.1f %%)' %( rcount, tablelen, pce )
		
		masklist = table.mask.as_array().tolist()
		tcount = np.sum( masklist ) 
		trcount = np.sum( [ np.sum( masklist[i] ) > 0  for i in xrange( tablelen ) ] )

		tpc = tcount / size * 100.
		tpce = trcount / tablelen * 100.
		print '* Total of %i of %i elements masked (%2.1f %%)' %( tcount, size, tpc )
		print '* Total of %i of %i epochs affected (%2.1f %%)' %( trcount, tablelen, tpce)

	return table


# Estimate x when lower sig passes y=0 for a set of model
def estzero( A, Asig, B, Bsig, x0, filters = None ):
	Avar = Asig**2
	Bvar = Bsig**2
	
	p2 = B**2 - Bvar
	p1 = A * B
	p0 = A**2 - Avar

	x = - p1 / p2
	q = x * x - p0 / p2

	xp = x + np.sqrt( q )
	xm = x - np.sqrt( q )

	fz = np.nanmax( [ xp, xm ], 0 )
	mask = np.greater( x, x0 )
	if mask.sum() > 0:
		print xm, xp, mask
		fz[mask] = np.nanmin(  xm, xp, 0 )[mask]
	sel = fz == np.nanmax( fz )

	fzero = ( fz[ sel ] + x0[ sel ] )[0]

	if filters is None:
		return fzero

	else:
		filterzero = filters[ sel ][0]
		return fzero, filterzero


