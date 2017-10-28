# AGN Project Function Library
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
