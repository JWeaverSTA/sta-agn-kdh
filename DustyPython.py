# DustyPyton
# For decomposing AGN lightcurves
# Created 23.10.17

# Imports
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
from matplotlib.ticker import ScalarFormatter, NullFormatter, MaxNLocator
import astro_funclib as astrotools
import adv_funclib as advtools
import agn_funclib as agntools
import dust_funclib as dusttools
import time

# Input Parameters
lightcurvedir = 'QSO_S82/'
masterdir = 'QSO_Master/'
masterfile = 'DB_QSO_S82.dat'
outputdir = 'Output/'
outputfile = 'DustyOutput.trial.cat'
outcolnamesfile = 'colnames.txt'

# SHOULD INCLUDE WARNING CHECK THAT OUTPUT FILE DOES NOT ALREADY EXIST!


fidwav = 2400. # AA
conlimit = 4
concrit = 1E-10
conthresh = 0.9

markersize = 10
emarkersize = 10

# Debugging tools
# Extra verbosity in lc fitting
verbose1 = False
# Convergence plotter
convplotter = False


# Data Associated Arrays
# Dustlaws
dustlaws = ( 'smc', 'lmc', 'mw', 'agn' )
dustlaws = np.array( dustlaws )
# Filters (in order)
filters = ( 'u', 'g', 'r', 'i', 'z' )
filters = np.array( filters )
flen = len( filters )
# Wavelengths
wav = ( 3543., 4770., 6231., 7625., 9134. )
wav = np.array( wav )
# Wavlength colors
wavcolor = { 'u' : 'b',
			 'g' : 'g',
			 'r' : 'orange',
			 'i' : 'r',
			 'z' : 'k'}
wcolor = ( 'b', 'g', 'orange', 'r', 'k' )
# Extinction Coefficients (Mag)
extmag_coeff = { 'u_mag' : 1.0,
				 'g_mag' : 0.736,
				 'r_mag' : 0.534,
				 'i_mag' : 0.405,
				 'z_mag' : 0.287 }
# Emission lines
chem = { 'Ly$_{\infty}$'  : 912.00,
		 'Ly$_{\\alpha}$' : 1215.67,
		 'CIV'            : 1549.05,
		 '[CIII]'         : 1908.73,
		 '[MgII]'         : 2797.92,
		 'H$_{\infty}$'   : 3560.00,
		 'H$_{\\beta}$'   : 4861.32,
		 'H$_{\\alpha}$'  : 6562.80}

# Load dustlaw interpolations
tk_smc = dusttools.smccubic()
tk_lmc = dusttools.lmccubic()

# Load Masterfile
masterdata = ascii.read( table = masterdir + masterfile )
masterdata['ra'] = astrotools.deg2hr( masterdata['ra'] )
nagn = np.arange( 1, len( masterdata ) + 1 )
nagn_init = len( nagn ) + 1

# Main Menu Navigation
# Inital Sample
iiagn1 = 1
iiagn2 = len( nagn ) + 1
izlo = 0
izhi = 5
igbrt = 16
igfnt = 23
icblu = -1.
icred = 2.
iralo = 0
irahi = 24.
ideclo = -1.27
idechi = 1.27
inumlook = 0
ioutput = 'OFF'
iverbose = 'ON'

# Colors
gi_color = masterdata['g'] - masterdata['i']

# Transfer inital values for use in loop
iagn1 = iiagn1
iagn2 = iiagn2
zlo = izlo
zhi = izhi
gfnt = igfnt
gbrt = igbrt
cblu = icblu
cred = icred
ralo = iralo
rahi = irahi
declo = ideclo
dechi = idechi
numlook = inumlook
output = ioutput
verbose = iverbose

# Interface
truth_menu = 0
truth_nav = 0
while truth_nav not in ('end','exit','q'):

	while truth_menu not in ('P','L'):
	
		# Define Indices
		ind_agn = ( nagn >= iagn1 )             & ( nagn < iagn2 )
		ind_z   = ( masterdata['reds'] > zlo )  & ( masterdata['reds'] < zhi )
		ind_g   = ( masterdata['g'] > gbrt )    & ( masterdata['g'] < gfnt )
		ind_c   = ( gi_color > cblu )           & ( gi_color < cred )
		ind_ra  = ( masterdata['ra'] > ralo )   & ( masterdata['ra'] < rahi )
		ind_dec = ( masterdata['dec'] > declo ) & ( masterdata['dec'] < dechi )
		ind_nl  = ( masterdata['dbID'] == numlook )
	
		index = np.logical_and.reduce( [ ind_agn, ind_z, ind_g,
									ind_c, ind_ra, ind_dec ] )
							 
		# Override index
		if numlook != 0:
			index = ind_nl
		if iagn1 == iagn2:
			index = iagn1 - 1
		
		# Redefine table    
		selectdata = masterdata[index]
		
		if type( index ) == int:
			nselect = 1
		else:
			nselect = len( selectdata )

		# Report to user
		print ' '
		print 'SAMPLE SIZE: ', nselect
		print ' '
		print 'SELECTION CRITERIA:'
		print '  AGN:', iagn1, iagn2, 'OF', nagn_init
		print '  Z  :', zlo, zhi, ' REDSHIFT'
		print '  MAG:', gbrt, gfnt, ' SLOAN g'
		print '  COL:', cblu, cred, ' SLOAN g-i'
		print '  RA :', ralo, rahi, ' HR'
		print '  DEC:', declo, dechi, ' DEG'
		print 'SDSSN:', numlook
		print ' '
		print 'VERBOSE...'+verbose
		print 'OUTPUT....'+output
		if output == 'ON':
			print 'FILE NAME: ', outputfile
		print ' '
		if len( selectdata )>0:
			print 'P.....Print'
			print 'L.....List'
		if len( selectdata )==0:
			print 'ZERO SIZE SAMPLE - TRY AGAIN.'
		print ' '
	  
		# Ask for changes
		truth_menu = raw_input('CHANGE: ')
		if truth_menu in ('end','exit','q'):
			sys.exit(0)

		if truth_menu == 'AGN':
			iagn1 = int(raw_input('  lower: '))
			iagn2 = int(raw_input('  upper: '))
		if truth_menu == 'Z':
			zlo = float(raw_input('  lower: '))
			zhi = float(raw_input('  upper: '))
		if truth_menu == 'MAG':
			gbrt = float(raw_input('  lower: '))
			gfnt = float(raw_input('  upper: '))
		if truth_menu == 'COL':
			cblu = float(raw_input('  lower: '))
			cred = float(raw_input('  upper: '))
		if truth_menu == 'RA':
			ralo = float(raw_input('  lower: '))
			rahi = float(raw_input('  upper: '))
		if truth_menu == 'Dec':
			declo = float(raw_input('  lower: '))
			dechi = float(raw_input('  upper: '))
		if truth_menu == 'SDSSN':
			numlook = int(raw_input('  number: '))
		if truth_menu == 'VERBOSE':
			verbose = str(raw_input('  ON/OFF: '))      
			verbose == verbose.upper() 
		if truth_menu == 'OUTPUT':
			output = str(raw_input('  ON/OFF: '))
			output = output.upper()
		if truth_menu == 'reset':
			iagn1 = iiagn1
			iagn2 = iiagn2
			zlo = izlo
			zhi = izhi
			gfnt = igfnt
			gbrt = igbrt
			cblu = icblu
			cred = icred
			ralo = iralo
			rahi = irahi
			declo = ideclo
			dechi = idechi
			numlook = inumlook
			output = ioutput
			verbose = iverbose

		print ' '


		# Display list of AGN
		if truth_menu == 'L':
			print selectdata.more()
			truth_nav = raw_input('BACK: ')
			truth_menu = 0

		# Further Plotting
		if truth_menu == 'P':
			
			# Plotting Options
			truth_plot = 0
			while truth_plot not in ('S','M','end','exit','q'):
			
				print 'PRINTING OPTIONS:'
				print '   [S] Sample Properties'
				print '   [M] Lightcurve Model'
				print '   [A] ... Automatic'
				print ' '

				truth_plot = raw_input( 'CHOICE: ' )
				print ' ' 
				
				plt.ion()
				
				# --- SAMPLE PLOT ---
				if truth_plot == 'S':

					# Create Figure 1
					fig1 = plt.figure( figsize = ( 11, 6 ) )
					fig1.suptitle( '%i | SDSS Quasars' % len( selectdata ) )
					nbins = [50, 50]

					gi_model = astrotools.accretion_color( wav[1], wav[3] )

					# ax1: Redshift-Colour
					ax11 = fig1.add_subplot( 221 )
					ax11.hist2d( selectdata['reds'], selectdata['g'],
								 range = [ [ izlo, izhi ], [ igbrt, igfnt ] ],
								 bins = nbins )
					ax11.set_ylim( igfnt, igbrt )
					ax11.set_xlim( izlo, izhi )
					ax11.set_xlabel( 'redshift z' )
					ax11.set_ylabel( 'AB(g) (mag)' )

					# ax2 : Colour-Magnitude diagram
					ax12 = fig1.add_subplot( 222 )
					ax12.hist2d( selectdata['g'] - selectdata['i'], selectdata['g'],
								 range = [ [ icblu, icred ], [ igbrt, igfnt ] ],
								 bins = nbins )
					ax12.vlines( gi_model, igfnt, igbrt, colors='r', linestyles='--')
					ax12.set_ylim( igfnt, igbrt )
					ax12.set_xlim( icblu, icred )
					ax12.set_xlabel( 'AB(g-i) (mag)' )
					ax12.set_ylabel( 'AB(g) (mag)' )

					# ax3 : Redshift-Colour diagram
					ax13 = fig1.add_subplot( 223 )
					ax13.hist2d( selectdata['reds'], selectdata['g'] - selectdata['i'],
								 range = [ [ izlo, izhi ], [ icblu, icred ] ],
								 bins = nbins )
					ax13.hlines( gi_model, izlo, izhi, colors='r', linestyles='--' )
					ax13.set_ylim( icred, icblu )
					ax13.set_xlim( izlo, izhi ) 
					ax13.set_ylabel( 'AB(g-i) (mag)' )
					ax13.set_xlabel( 'redshift z' )

					# ax4 : RA-Dec diagram
					ax14 = fig1.add_subplot( 224 )
					ax14.hist2d( selectdata['ra'], selectdata['dec'],
								 range = [ [ iralo, irahi ], [ ideclo, idechi ] ],
								 bins = nbins )
					ax14.set_ylim( ideclo, idechi )
					ax14.set_xlim( iralo, irahi )
					ax14.set_ylabel( 'Dec (deg)' )
					ax14.set_xlabel( 'RA (h)' )

					# Must specifty layout POST plot
					fig1.tight_layout( rect=[0,0,1,0.9] )

					# Draw figure
					fig1.canvas.draw_idle()
					truth_plot = ''

					# Close figure
					truth_plot = raw_input( 'CLOSE: ' )
					if truth_plot == '':
					  plt.close('all')
					  print ' '


				# Lightcurves
				if truth_plot in ( 'M', 'A' ):

					# Remind user that output is off in automatic mode
					if ( truth_plot in 'A' ) & ( output == 'OFF' ):
						truth_output = raw_input( 'WARNING: OUTPUT IS OFF - CHANGE? [Y/n] ')
						if truth_output in ( 'Y', 'y' ):
							output = 'ON'
							print '* OUTPUT SET TO "ON"'

					# Remind user that verbose is on in automatic mode
					if ( truth_plot in 'A' ) & ( verbose == 'ON' ):
						truth_output = raw_input( 'WARNING: VERBOSE IS ON - CHANGE? [Y/n] ')
						if truth_output in ( 'Y', 'y' ):
							verbose = 'OFF'
							print '* VERBOSE SET TO "OFF"'

					# Translate verbosity
					if verbose == 'ON':
						verbose = True
					elif verbose == 'OFF':
						verbose = False     

					# init convergence plotter
					if convplotter:
						figconv = plt.figure()
						axconv = figconv.add_subplot(111)


					# init figures
					if truth_plot in 'M':
						fig2, ( ( ax21, ax22 ),
								( ax23, ax24 ),
								( ax25, ax26 ) ) = plt.subplots(nrows=3 ,ncols=2, figsize = (10,10), dpi = 100 )
						ax23t = ax23.twiny()
						ax24t = ax24.twiny()

						fig3, ( ( ax31, ax32 ),
								( ax33, ax34 ) ) = plt.subplots( nrows=2, ncols=2, figsize = (9.25, 5) )
						ax31t = ax31.twiny()
						ax32t = ax32.twiny()
						ax33t = ax33.twiny()
						ax34t = ax34.twiny()
					
					# Convert output
					if output == 'ON':
						output = True
					if output == 'OFF':
						output = False

					# Output file
					if output:
						outcolnames = ascii.read( outcolnamesfile, data_start = 0  )
						outputtable = Table( names = [ i[0] for i in outcolnames ],
											 dtype = ['i8', 'f8', 'f8', 'i8'] + ['f8' for i in range(109)] )

					# Loop over LC
					numtotal = len( selectdata )
					for lc_index, num in enumerate( xrange( numtotal ) ):
						
						# Get master data info
						lc_info = selectdata[lc_index]

						# Open raw lightcurve
						lc_id = lc_info['dbID']
						colnames = [ [ i + '_mjd', i + '_mag', i + '_merr' ] for i in filters ]
						colnames = sum( colnames, [] ) + [ 'RA', 'Dec' ]
						rawmag = ascii.read( table = lightcurvedir + str( lc_id ),
											 names = colnames )
						rlen = len( rawmag )

						# Update user
						num = num + 1
						pc = num / numtotal * 100.
						print ''
						print '%3.2f%% (%i/%i) - No. %i' %( pc, num, numtotal, lc_id )
						if verbose:
							print ''
							print 'Read: Opened file %i containing %i epochs' %( lc_id, rlen )
							print '-------------------------------------------------------------------'

						# Mask bad values
						cleanmag = agntools.maskval( intable = rawmag,
													  val = -99.99,
													  test = 'eq',
													  str_pattern = '_mag',
													  verbose = verbose )

						cleanmag = agntools.maskval( intable = cleanmag,
													  val = 30.0,
													  test = 'gt',
													  str_pattern = '_mag',
													  verbose = verbose )

						cleanmag = agntools.maskval( intable = cleanmag,
													  val = 10.00,
													  test = 'lt',
													  str_pattern = '_mag',
													  verbose = verbose )

						cleanmag = agntools.maskval( intable = cleanmag,
													  val = 1.0,
													  test = 'gt',
													  str_pattern = '_merr',
													  verbose = verbose )

						# Correct for dust extinction
						extmag_ref = lc_info['Au']
						corrmag = astrotools.extcorr( intable = cleanmag,
													  extmag_ref = extmag_ref,
													  extmag_coeff = extmag_coeff,
													  verbose = verbose,
													  filters = filters )

						# Create table of fluxes
						fluxtable = Table()
						if verbose:
							print 'Fluxtable: Created table of fluxes'

						avgepoch = np.mean([ corrmag[ i + '_mjd'] for i in filters ], 0)
						fluxtable.add_column( MaskedColumn( data = avgepoch,
															name = 'mjd' ) )

						if verbose:
							print '* Added MJD column'
						for i in filters:
							f, fsig = astrotools.ab2mJy( corrmag[ i + '_mag' ], corrmag[ i + '_merr'] )
							fcol = MaskedColumn( data = f,
												 name = i + '_f' )
							fsigcol = MaskedColumn( data = fsig,
													name = i + '_fsig' )
							fluxtable.add_columns( [ fcol, fsigcol ] ) 

							if verbose:
								avgf, avgfsig = advtools.optavg( fcol, fsigcol )
								print '* Average %s flux ... %1.4f +/- %1.4f mJy' %( i, avgf, avgfsig )

						
						# Fill with nans and remove all bad rows
						fluxtable = fluxtable.filled( np.nan )
						corrmag = corrmag.filled( np.nan )
						remrows = agntools.removenan( fluxtable )
						if remrows is not None:
							fluxtable.remove_rows( remrows )
							corrmag.remove_rows( remrows )
						elen = len( fluxtable )
						nrem = rlen - elen
						rpc = nrem / elen * 100

						# Update user
						if verbose:
							print 'Fluxtable: Removed %i all-NaN epochs ( %3.1f%% ) - %i remaining' %( nrem, rpc, elen ) 

						# Fit lightcurves
						# Update user
						if verbose:
							print 'Lightcurve: Fitting for model -> A(w) + B(w) * L(t)'

						# inital values
						dlen = len( filters )
						xlen = len( fluxtable )

						A = np.zeros( dlen )
						B = np.zeros( dlen )
						Asig = np.zeros( dlen )
						Bsig = np.zeros( dlen )

						L = np.ones( xlen )
						Lsig = np.ones( xlen )
						L0 = np.ones( dlen )

						loopnum = 0
						conval = 999.
						contrack = 0
						conabandon = False

						# Setup convergence plotter
						if convplotter:
							axconv.clear()
							axconv.set_xlabel( 'Trial' )
							axconv.set_ylabel( 'Criterion' )
							axconv.set_yscale( 'log' )
							axconv.set_xlim( 0, 20 )
							axconv.hlines( concrit, 0, 20,
										   linestyles = 'dotted', colors = 'k' )

						# Loop over filters for inital guesses for A(w), B(w)
						if verbose1:
							print 'Optimal Average: Initalising A(w), B(w) ...'
						for i, j in enumerate( filters ):
							x = fluxtable[ j +'_f' ]
							xsig = fluxtable[ j + '_fsig' ]
							A[i], Asig[i] = advtools.optavg( x, xsig )
							B[i], Bsig[i] = np.sqrt( advtools.optavg( ( x - A[i] )**2, xsig ) )

							if verbose1:
								print '* A(%s) = %2.5f +/- %2.5f | B(%s) = %2.5f +/- %2.5f' %( j, A[i], Asig[i], j, B[i], Bsig[i] )

						# Loop until convergence
						while ( conval > concrit ) & ( not conabandon ):

							# Add to loop:
							loopnum += 1

							# Update user                            
							if verbose1:
								print 'Loop %i | Criterion ... %2.1e > %2.1e' %( loopnum, conval, concrit )

							# Store old L
							oldL = np.mean( L )

							# Scale to find L(t)
							if verbose1:
								print 'Optimal Scaling: Calculating L(t) ...'
							for i, fluxrow in enumerate( fluxtable ):
								t = np.array( [ fluxrow[ j +'_f' ] for j in filters  ] )
								tsig = np.array( [ fluxrow[ j + '_fsig' ] for j in filters ] )
								L[i], Lsig[i] = advtools.optscl( t - A, tsig, B )

								if verbose1:
									mjd = fluxrow['mjd']
									print '* L(%.1f) = %2.1f +/- %2.1f' %( mjd, L[i], Lsig[i] )

							# Normalise L(t)
							if verbose1:
								print 'Normalisation: Fixing L(t)'
							avgL = advtools.optavg( L, Lsig )[0]
							L = L - avgL
							avgL_squ = advtools.optavg( L**2, Lsig )[0]
							L = L / np.sqrt( avgL_squ )
							if verbose1:
								newavgL = advtools.optavg( L, Lsig )[0]
								print '* <L> = %2.1f | <L^2> = %2.1f ... <L> = %2.1f' %( avgL, avgL_squ, newavgL )

							# Solve for A(w), B(w)
							if verbose1:
								print ' Fitline: Calculating A(w), B(w) ...'
							for i, j in enumerate( filters ):
								x = fluxtable[ j +'_f' ]
								xsig = fluxtable[ j + '_fsig' ]
								A[i], B[i], Asig[i], Bsig[i], _, _, L0[i] = advtools.fitline( L, x, xsig )

								if verbose1:
									print '* A(%s) = %2.5f +/- %2.5f | B(%s) = %2.5f +/- %2.5f' %( j, A[i], Asig[i], j, B[i], Bsig[i] )

							# Compute convergence value
							Ldiff = abs( np.mean( L ) - oldL )
							avgLsig = np.mean( Lsig )
							oldconval = conval
							conval = Ldiff / avgLsig

							# Track convergence is OK
							dconval = conval / oldconval
							if dconval > conthresh:
								contrack += 1
								if verbose:
									print oldconval, conval, dconval, conthresh
									print 'Fitline: convergence did not progress ( C(i)/C(i-1) = %2.1e ) in loop %i, %i times' %( dconval, loopnum, contrack )


							if contrack > conlimit:
								if verbose:
									print 'Fitline: convergence limit %i exceeded - converged to %2.1e > %2.1e ' %( conlimit, conval, concrit )
								conabandon = True
								continue

							# Update plot
							if convplotter:
								axconv.scatter( loopnum, conval, c='k' )
								plt.pause(0.0001)

						# Abandon analysis for object
						if conabandon:
							print 'Fitline: analysis abandoned for %i' %lc_id
							continue

						# Update user
						if verbose:
							print 'Lightcurve: Convergence succeeded on loop %i' %( loopnum )
							print 'Lightcurve: Final Values for %i' %( lc_id ) 
							for i, j in enumerate( filters ):
								print '* A(%s) = %2.5f +/- %2.5f | B(%s) = %2.5f +/- %2.5f' %( j, A[i], Asig[i], j, B[i], Bsig[i] )

						# Calculate Fzero
						F0, filterzero = agntools.estzero( A, Asig, B, Bsig, L0, filters )

						# Update user
						if verbose:
							print 'Lightcurve: Flux zeropoint found at L(t) = %3.2f in %s' %( F0, filterzero )

						# Calculate component spectra
						# Update user
						if verbose:
							print 'Lightcurve: Calculating component spectra'

						# Create table
						spectable = Table()

						# Calculate variances
						Avar = Asig**2
						Bvar = Bsig**2

						# Precalculations
						Lmax = np.nanmax( L )
						Lmin = np.nanmin( L )

						# Update user
						if verbose:
							print '* L(fnt) = %5.2f , L(brt) = %5.2f' %(Lmin, Lmax)

						# Mean spectrum
						mean_spec = A
						mean_col = Column( name = 'mean', data = mean_spec )

						meansig_spec = Asig
						meansig_col = Column( name = 'mean_sig', data = meansig_spec )

						# RMS spectrum
						rms_spec = B
						rms_col = Column( name = 'rms', data = rms_spec )

						rmssig_spec = Bsig
						rmssig_col = Column( name = 'rms_sig', data = rmssig_spec )

						# Galaxy spectrum
						gal_spec = A + B * ( F0 - L0 )
						gal_col = Column( name = 'gal', data = gal_spec )

						galsig_spec = np.sqrt( Avar + Bvar * ( F0 - L0 )**2 )
						galsig_col = Column( name = 'gal_sig', data = galsig_spec )

						# Bright spectrum
						brt_spec = A + B * ( Lmax - L0 )
						brt_col = Column( name = 'brt', data = brt_spec )

						brtsig_spec = np.sqrt( Avar + Bvar * ( Lmax - L0 )**2 )
						brtsig_col = Column( name = 'brt_sig', data = galsig_spec )

						# Faint spectrum
						fnt_spec = A + B * ( Lmin - L0 )
						fnt_col = Column( name = 'fnt', data = fnt_spec )

						fntsig_spec = np.sqrt( Avar + Bvar * ( Lmin - L0 )**2 )
						fntsig_col = Column( name = 'fnt_sig', data = fntsig_spec )

						# Difference spectrum
						diff_spec = B * ( Lmax - Lmin )
						diff_col = Column( name = 'diff', data = diff_spec )

						diffsig_spec = np.sqrt( Bvar * ( Lmax - Lmin )**2 )
						diffsig_col = Column( name = 'diff_sig', data = diffsig_spec )

						# Variable spectrum
						var_spec = A + B * ( L0 - F0 )
						var_col = Column( name = 'var', data = var_spec )

						varsig_spec = np.sqrt( Avar + Bvar * ( L0 - F0 ) )
						varsig_col = Column( name = 'var_sig', data = varsig_spec )

						# disc spectrum
						disc_spec = B * ( ( brt_spec - A ) / B - F0 )
						disc_col = Column( name = 'disc', data = disc_spec )

						discsig_spec = np.sqrt( Bvar ) * ( ( brt_spec - A ) / B - F0 )
						discsig_col = Column( name = 'disc_sig', data = discsig_spec )

						# Append to table
						spectable.add_columns( [ mean_col, meansig_col,
												 rms_col, rmssig_col,
												 gal_col, galsig_col,
												 brt_col, brtsig_col,
												 fnt_col, fntsig_col,
												 diff_col, diffsig_col,
												 var_col, varsig_col,
												 disc_col, discsig_col ] )

						# Update user
						if verbose1:
							for i in spectable.colnames[::2]:
								print '* %s spectrum ' %i
								for j, k in enumerate( filters ):
									val = spectable[ i ][ j ]
									valsig = spectable[ i + '_sig'][ j ]
									print '** F(%s) = %5.5f +/- %5.5f mJy' %( k, val, valsig )
						elif verbose:
							for i in ('mean', 'disc', 'gal'):
								print '* %s spectrum ' %i
								for j, k in enumerate( filters ):
									val = spectable[ i ][ j ]
									valsig = spectable[ i + '_sig'][ j ]
									print '** F(%s) = %5.5f +/- %5.5f mJy' %( k, val, valsig )

						# Fit Dust Models
						# Calculate restwavelengths
						z = lc_info['reds']
						wavz = astrotools.obs2restwav( wav, z )

						# Update user
						if verbose:
							print 'Dustfit: Found rest wavelengths'
							print '* Redshift = %3.2f' %z
							for i, j in enumerate( filters ):
								print '* %s observed light at %5.2f A' %(j, wavz[i])


						# Convert to mag
						brt_mag, brtsig_mag = astrotools.mJy2ab( spectable['brt'], spectable['brt_sig'] )
						fnt_mag, fntsig_mag = astrotools.mJy2ab( spectable['fnt'], spectable['fnt_sig'] )
						diff_mag, diffsig_mag = astrotools.mJy2ab( spectable['diff'], spectable['diff_sig'] )
						rms_mag, rmssig_mag = astrotools.mJy2ab( spectable['rms'], spectable['rms_sig'] )

						# Covert to absolute mag
						mean_mag, meansig_mag = astrotools.mJy2ab( spectable['mean'], spectable['mean_sig'])
						mean_amag, meansig_amag = astrotools.mag2amag( mean_mag, meansig_mag, z )

						disc_mag, discsig_mag = astrotools.mJy2ab( spectable['disc'], spectable['disc_sig'])
						disc_amag, discsig_amag = astrotools.mag2amag( disc_mag, discsig_mag, z )

						gal_mag, galsig_mag = astrotools.mJy2ab( spectable['gal'], spectable['gal_sig'])
						gal_amag, galsig_amag = astrotools.mag2amag( gal_mag, galsig_mag, z )

						# Update user
						if verbose:
							print 'Dustfit: Prepared disc spectrum flux -> absolute mag'
							if verbose1:
								for i, j in enumerate( filters ):
									print '* Amag(%s) = %5.5f +/- %5.5f mag' %( j, disc_amag[i], discsig_amag[i] )


						# Calculate accretion spectrum
						acc_amag = astrotools.accretion_magspec( wav = wavz,
																fidwav = fidwav )
						residuals = disc_amag - acc_amag

						# Update user
						if verbose:
							print 'Dustfit: Calculated accretion spectrum'
							if verbose1:
								for i, j in enumerate( filters ):
									print '* Amag(%s) = %5.5f mag' %( j, acc_amag[i] )

							print 'Dustfit: Calculated residual spectrum'
							if verbose1:
								for i, j in enumerate( filters ):
									print '* Amag(%s) = %5.5f +/- %5.5f mag' %( j, residuals[i], discsig_amag[i] )

						# Output
						if output:
							outvals = [ i for i in lc_info ]
							outvals = outvals + mean_amag.tolist() + meansig_amag.tolist()
							outvals = outvals + disc_amag.tolist() + discsig_amag.tolist()
							outvals = outvals + gal_amag.tolist() + galsig_amag.tolist()

						# Plotting
						if truth_plot == 'M':

							# Update user
							if verbose:
								print 'Plot: Beginning plotting...'

							# Preparation

							xlo1 = - abs( 1.1 * F0 )
							xhi1 = abs( 0.5 * F0 ) + max( L )
							ylo1 = - 0.05 * np.max( A + B * xhi1 )
							yhi1 = np.max( A + B * xhi1 )

							xlen = 1000
							model_Lx = np.linspace( xlo1, xhi1, xlen )


							lc_mag = np.array( [ corrmag[ i + '_mag' ] for i in filters ] )
							lcsig_mag = np.array( [ corrmag[ i + '_merr' ] for i in filters ] )
							lc_mag_max = np.nanmax( lc_mag + lcsig_mag )
							lc_mag_min = np.nanmin( lc_mag - lcsig_mag ) 

							sc_mag_max = np.nanmax( [ brt_mag + brtsig_mag,
													  fnt_mag + fntsig_mag,
													  diff_mag + diffsig_mag,
													  rms_mag + rmssig_mag ] )
							sc_mag_min = np.nanmin( [ brt_mag - brtsig_mag,
													  fnt_mag - fntsig_mag,
													  diff_mag - diffsig_mag,
													  rms_mag - rmssig_mag ] )

							gd_mag_max = np.nanmax( [ disc_mag + discsig_mag,
													  gal_mag + galsig_mag ] )
							gd_mag_min = np.nanmin( [ disc_mag - discsig_mag,
													  gal_mag - galsig_mag ] )

							wrange = wav[-1] - wav[0]
							xlo34 = wav[0] - 0.05 * wrange
							xhi34 = wav[-1] + 0.15 * wrange
							xlo34z, xhi34z = astrotools.obs2restwav( [ xlo34, xhi34 ], z )

							sc_range = sc_mag_max - sc_mag_min
							ylo3 = sc_mag_min - 0.15 * sc_range
							yhi3 = sc_mag_max + 0.10 * sc_range

							gd_range = gd_mag_max - gd_mag_min
							ylo4 = gd_mag_min - 0.15 * gd_range
							yhi4 = gd_mag_max + 0.10 * gd_range


							# Ax21 - Model
							ax21.set_title( 'Model')
							ax21.set_xlabel( 'L(t)' )
							ax21.set_ylabel( 'f($\lambda$) (mJy)' )
							ax21.text( 0.05, 1.07, '(a)', transform = ax21.transAxes,
									   horizontalalignment = 'center', verticalalignment = 'center' )

							ax21.set_xlim( xlo1, xhi1 )
							ax21.set_ylim( ylo1, yhi1 )

							ax21.vlines( F0, ylo1, yhi1, linestyles = 'dashed' )
							ax21.vlines( 0.0, ylo1, yhi1, zorder = -1 )
							ax21.hlines( 0.0, xlo1, xhi1 )
							ax21.vlines( max( L ), ylo1, yhi1, linestyle = 'dotted' )
							ax21.vlines( min( L ), ylo1, yhi1, linestyle = 'dotted' )

							ax21.text( 0.05, 1.07, '(a)', transform=ax21.transAxes,
									  horizontalalignment='center', verticalalignment='center' )                        
							ax21.text( 0.08, 0.93, 'ID: %i' %lc_id, transform=ax21.transAxes,
									  horizontalalignment='left', verticalalignment='center' )
							ax21.text( 0.08, 0.85, 'z=%.2f' %z, transform=ax21.transAxes,
									  horizontalalignment='left', verticalalignment='center' )

							# Ax22 - Lightcurve
							ax22.set_title( 'Lightcurve' )
							ax22.set_xlabel( 'MJD' )
							ax22.set_ylabel( 'M($\lambda$) (mag)')
							ax22.text( 0.05, 1.07, '(b)', transform = ax22.transAxes,
									   horizontalalignment = 'center', verticalalignment = 'center' )


							# Ax23 - Spectral Components
							ax23.set_title( 'Spectral Components', y=1.1)
							ax23.set_xlabel( '$\lambda$ = $\lambda_{o}$(1 + z) ($\AA$)' )
							ax23.set_ylabel( 'M($\lambda$) (mag)' )
							ax23.set_xscale( 'log' )
							ax23.xaxis.set_major_formatter( ScalarFormatter() )
							ax23.xaxis.set_minor_formatter( ScalarFormatter() )
							ax23.text( 0.05, 1.07, '(c)', transform = ax23.transAxes,
									   horizontalalignment = 'center', verticalalignment = 'center' )

							ax23t.set_xlim( xlo34z, xhi34z )
							ax23t.set_xscale( 'log' )
							ax23t.xaxis.set_major_formatter( ScalarFormatter() )
							ax23t.xaxis.set_minor_formatter( ScalarFormatter() )

							ax23.set_xlim( xlo34, xhi34 )
							ax23.set_ylim( ylo3, yhi3 )	
							ax23.invert_yaxis()



							# Plot spectra
							color23 = ( 'k', 'b', 'b', 'b', 'k', )
							name23 = ( 'mean', 'brt', 'fnt', 'diff', 'rms' )
							i_mag = ( mean_mag, brt_mag, fnt_mag, diff_mag, rms_mag )
							sig_mag = ( meansig_mag, brtsig_mag, fntsig_mag, diffsig_mag, rmssig_mag )
							for i in xrange( len( i_mag ) ):
								color = color23[i]
								name = name23[i]
								p_mag, psig_mag = i_mag[i], sig_mag[i]
								ax23.fill_between( wav, p_mag - psig_mag, p_mag + psig_mag,
												   alpha = 0.2, color = color )
								ax23.plot( wav, p_mag, marker = 'None', c = color, zorder = 1 )
								ax23.scatter( wav, p_mag, c = wcolor, s = markersize, zorder = 2 )
								ax23.text( 1.01 * xhi34, p_mag[-1], name, color = color )

							# Plot chemical lines
							chemname = np.array( chem.keys() )
							chemwavz = np.array( chem.values() )

							mask = ( ( chemwavz > 1.07 * xlo34z ) & ( chemwavz < 0.93 * xhi34z ) )
							ax23t.vlines( chemwavz[mask], ylo3 + 0.15 * sc_range, yhi3,
										 colors = 'k', linestyles = 'dotted', zorder = 0 )
							chemwavz = chemwavz[mask]
							chemname = chemname[mask]
							for i in xrange( len( chemwavz ) ):
								pwav = chemwavz[i]
								pname = chemname[i]
								ax23t.text( pwav, ylo3 + 0.1 * sc_range, pname, size = 10,
											color = 'k', horizontalalignment = 'center' )


							# Ax24 - Gal-Disc
							ax24.set_title( 'Galaxy-Disc Separation', y=1.1 )
							ax24.set_xlabel( '$\lambda$ = $\lambda_{o}$(1 + z) ($\AA$)' )
							ax24.set_ylabel( 'M($\lambda$) (mag)' )
							ax24.set_xscale( 'log')
							ax24.xaxis.set_major_formatter( ScalarFormatter() )
							ax24.xaxis.set_minor_formatter( ScalarFormatter() )
							ax24.text( 0.05, 1.07, '(d)', transform = ax24.transAxes,
									   horizontalalignment = 'center', verticalalignment = 'center' )

							ax24t.set_xlim( xlo34z, xhi34z )
							ax24t.set_xscale( 'log' )
							ax24t.xaxis.set_major_formatter( ScalarFormatter() )
							ax24t.xaxis.set_minor_formatter( ScalarFormatter() )

							ax24.set_xlim( xlo34, xhi34 )
							ax24.set_ylim( ylo4, yhi4 )	
							ax24.invert_yaxis()

							# Plot spectra
							color24 = ( 'k', 'r' )
							name24 = ( 'disc', 'gal' )
							i_mag = ( disc_mag, gal_mag )
							sig_mag = ( discsig_mag, galsig_mag )
							for i in xrange( len( i_mag ) ):
								color = color24[i]
								name = name24[i]
								p_mag, psig_mag = i_mag[i], sig_mag[i]
								ax24.fill_between( wav, p_mag - psig_mag, p_mag + psig_mag,
												   alpha = 0.2, color = color )
								ax24.plot( wav, p_mag, marker = 'None', c = color, zorder = 1 )
								ax24.scatter( wav, p_mag, c = wcolor, s = markersize, zorder = 2 )
								ax24.text( 1.01 * xhi34, p_mag[-1], name, color = color )

							# Plot chemical lines
							chemname = np.array( chem.keys() )
							chemwavz = np.array( chem.values() )

							mask = ( ( chemwavz > 1.07 * xlo34z ) & ( chemwavz < 0.93 * xhi34z ) )
							ax24t.vlines( chemwavz[mask], ylo4 + 0.15 * gd_range, yhi4,
										 colors = 'k', linestyles = 'dotted', zorder = 0 )
							chemwavz = chemwavz[mask]
							chemname = chemname[mask]
							for i in xrange( len( chemwavz ) ):
								pwav = chemwavz[i]
								pname = chemname[i]
								ax24t.text( pwav, ylo4 + 0.1 * gd_range, pname, size = 10,
											color = 'k', horizontalalignment = 'center' )
						

							# Ax25 - Time Sequence
							ax25.set_title( 'Time Sequence' )
							ax25.set_xlabel( 'Time Sequence Number' )
							ax25.set_ylabel( 'f($\lambda$) (mJy)' )
							ax25.text( 0.05, 1.07, '(e)', transform = ax25.transAxes,
									   horizontalalignment = 'center', verticalalignment = 'center' )




							# Ax26 - Date Sequence
							ax26.set_title( 'Date Sequence' )
							ax26.set_xlabel( 'MJD' )
							ax26.set_ylabel( 'f($\lambda$) (mJy)' )
							ax26.text( 0.05, 1.07, '(f)', transform = ax26.transAxes,
									   horizontalalignment = 'center', verticalalignment = 'center' )


							# Loop over colors
							for i, j in enumerate( filters ):

								color = wavcolor[j]

								epoch = fluxtable[ 'mjd' ]
								lc_flux_col = fluxtable[ j + '_f' ]
								lcsig_flux_col = fluxtable[ j + '_fsig' ]
								model_mask = lc_flux_col == lc_flux_col

								lc_mag_col = lc_mag[i]
								lcsig_mag_col = lcsig_mag[i]

								model_flux = A[i] * np.ones( xlen ) + B[i] * ( model_Lx - L0[i] )
								modelsig_flux = np.sqrt( Avar[i]  * np.ones( xlen ) + Bvar[i] * ( model_Lx - L0[i] )**2 )

								model_mag, modelsig_mag = astrotools.mJy2ab( model_flux, modelsig_flux )

								L_col = L - L0[i]

								modelev_flux = A[i] * np.ones( elen ) + B[i] * L_col
								modelevsig_flux = np.sqrt( Avar[i] * np.ones( elen ) + Bvar[i] * L_col )									

								modelev_flux = modelev_flux[model_mask]
								modelevsig_flux = modelevsig_flux[model_mask]

								number = np.arange( 1, len( epoch ) + 1 )
								model_number = number[model_mask]

								model_epoch = epoch[model_mask]

								# Ax21
								ax21.errorbar( x = L_col, y = lc_flux_col, yerr = lcsig_flux_col, alpha = 0.8,
											   fmt = '.', color = color, ms = emarkersize, zorder = - i )
								ax21.plot( model_Lx, model_flux, c = color )
								ax21.fill_between( x = model_Lx,
												   y1 = model_flux - modelsig_flux,
												   y2 = model_flux + modelsig_flux,
												   alpha = 0.2,
												   color = color, zorder = - 10 )

								# Ax22
								ax22.errorbar( x = epoch, y = lc_mag_col, yerr = lcsig_mag_col, alpha = 0.8,
											   fmt = '.', color = color, ms = emarkersize, zorder = - i )

								# Ax25
								model_mask = ( lc_flux_col != lc_flux_col )
								ax25.errorbar( x = number, y = lc_flux_col, yerr = lcsig_flux_col, alpha = 0.8,
											   fmt = '.', color = color, ms = emarkersize, zorder = - i )
								ax25.plot( model_number, modelev_flux, c = color, zorder = 1 )
								ax25.fill_between( x = model_number,
												   y1 = modelev_flux - modelevsig_flux,
												   y2 = modelev_flux + modelevsig_flux,
												   alpha = 0.2,
												   color = color, zorder = - 10 )

								# Ax26
								ax26.errorbar( x = epoch, y = lc_flux_col, yerr = lcsig_flux_col, alpha = 0.8,
											   fmt = '.', color = color, ms = emarkersize, zorder = - i )
								ax26.plot( model_epoch, modelev_flux, c = color, zorder = 1 )
								ax26.fill_between( x = model_epoch,
												   y1 = modelev_flux - modelevsig_flux,
												   y2 = modelev_flux + modelevsig_flux,
												   alpha = 0.2,
												   color = color, zorder = - 10 )	

							# Update plots
							fig2.subplots_adjust( top = 0.95, bottom = 0.05, left = 0.08, right = 0.95, hspace = 0.5, wspace = 0.3 )
						  
							fig2.canvas.draw_idle() 


						# Prepare for dustlaws

						wrange = wav[-1] - wav[0]
						xlo = wav[0] - 0.05 * wrange
						xhi = wav[-1] + 0.15 * wrange
						xloz, xhiz = astrotools.obs2restwav( [ xlo, xhi ], z )

						model_wavz = np.linspace( xloz, xhiz, 1000 )

						# Loop over dustlaws
						for dust in dustlaws:	
							if dust == 'smc':
								extmag = dusttools.dustlaw( wavz, name = 'smc', tk = tk_smc )
								model_extmag = dusttools.dustlaw( model_wavz, name = 'smc', tk = tk_smc )
								axi = ax31
								axit = ax31t
								axi.tick_params( labelbottom='off' )
								axi.set_ylabel( 'M($\lambda$) (mag)' )
							if dust == 'lmc':
								extmag = dusttools.dustlaw( wavz, name = 'lmc', tk = tk_lmc )
								model_extmag = dusttools.dustlaw( model_wavz, name = 'lmc', tk = tk_lmc )
								axi = ax32
								axit = ax32t
								axi.tick_params( labelbottom='off' )
								axi.set_ylabel( 'M($\lambda$) (mag)' )
							if dust == 'mw':
								extmag = dusttools.dustlaw( wavz, name = 'mw' )
								model_extmag = dusttools.dustlaw( model_wavz, name = 'mw' )
								axi = ax33
								axit = ax33t
								axi.set_ylabel( 'M($\lambda$) (mag)' )
								axi.set_xlabel( '$\lambda$ = $\lambda_{o}$(1 + z) ($\AA$)' )
							if dust == 'agn':
								extmag = dusttools.dustlaw( wavz, name = 'agn' )
								model_extmag = dusttools.dustlaw( model_wavz, name = 'agn' )
								axi = ax34
								axit = ax34t
								axi.set_ylabel( 'M($\lambda$) (mag)' )
								axi.set_xlabel( '$\lambda$ = $\lambda_{o}$(1 + z) ($\AA$)' )

							# Update user
							if verbose:
								print 'Dustfit: Loaded %s extinction spectrum' %(dust)
								if verbose1:
									for i, j in enumerate( filters ):
										print '** Ext(%s) = %5.5f mag' %( j, extmag[i] )

							# Fit line
							fidmag, ebmv, fidmagsig, ebmvsig, covar, corr, w0 = advtools.fitline( x = extmag,
																								  y = residuals,
																								  sig = discsig_mag,
																								  orthogonalize = False )

							# Update user
							if verbose:
								print '* Calculated fit parameters'
								print '** E(B-V) = %5.2f +/- %5.2f | Mfid(%5.1f) = %5.2f +/- %5.2f' %( ebmv, ebmvsig, fidwav, fidmag, fidmagsig )


							# dereddened spectrum
							dered_amag = disc_amag - extmag * ebmv
							deredsig_amag = discsig_amag

							# Update user
							if verbose:
								print '* Calculated dereddened spectrum'
								if verbose1:
									for i, j in enumerate( filters ):
										print '* Amag(%s) = %5.5f +/- %5.5f mag' %( j, dered_amag, deredsig_amag[i] )


							# Calculate std, chisq/N
							mod_dered_amag = astrotools.accretion_magspec( wavz, fidwav ) + fidmag
							std = advtools.std( dered_amag, mod_dered_amag, deredsig_amag )
							chisq = advtools.chisq( dered_amag, mod_dered_amag, deredsig_amag )
							reduced_chisq = chisq / 5.

							# Update user
							if verbose:
								print '* Calculated model agreement parameters'
								print '** std = %5.2f mag | corr = %5.2f | chisq/N = %5.2f' %( std, corr, reduced_chisq )


							# Output
							if output:
								outvals = outvals + dered_amag.tolist() + deredsig_amag.tolist()
								outvals = outvals + [ ebmv, ebmvsig, fidmag, fidmagsig, std, corr, reduced_chisq ]

						
						
							# Plotting
							if truth_plot:

								# Preparation
								ax_mag_max = np.nanmax( [ disc_amag + discsig_amag,
														   dered_amag + deredsig_amag ] )
								ax_mag_min = np.nanmin( [ disc_amag - discsig_amag,
														   dered_amag - deredsig_amag ] )

								ax_range = ax_mag_max - ax_mag_min
								ylo = ax_mag_min - 0.10 * ax_range
								yhi = ax_mag_max + 0.10 * ax_range

								model_acclaw = astrotools.accretion_magspec( model_wavz, fidwav )
								model_disc_amag = model_acclaw + fidmag + model_extmag * ebmv
								model_discsig_amag = np.sqrt( model_extmag**2 * ebmvsig**2 + fidmagsig**2 + 2 * model_extmag * covar )
								model_dered_amag = model_acclaw + fidmag
								model_deredsig_amag = model_discsig_amag

								# Chemical lines
								chemname = np.array( chem.keys() )
								chemwavz = np.array( chem.values() )

								mask = ( ( chemwavz > 1.07 * xloz ) & ( chemwavz < 0.93 * xhiz ) )
								axit.vlines( chemwavz[mask], ylo + 0.15 * ax_range, yhi,
											 colors = 'k', linestyles = 'dotted', zorder = 0 )
								chemwavz = chemwavz[mask]
								chemname = chemname[mask]
								for i in xrange( len( chemwavz ) ):
									pwav = chemwavz[i]
									pname = chemname[i]
									axit.text( pwav, ylo + 0.1 * ax_range, pname, size = 10,
												color = 'k', horizontalalignment = 'center' )


								# Plot
								axi.set_xlim( xlo, xhi )
								axi.set_ylim( ylo, yhi )
								axi.set_xscale( 'log' )
								axi.xaxis.set_major_formatter( ScalarFormatter() )
								axi.xaxis.set_minor_formatter( ScalarFormatter() )
								axi.invert_yaxis()

								axit.set_xlim( xloz, xhiz )
								axit.set_xscale( 'log' )


								axit.plot( model_wavz, model_disc_amag, color = 'b' )
								axit.fill_between( x = model_wavz,
												   y1 = model_disc_amag - model_discsig_amag,
												   y2 = model_disc_amag + model_discsig_amag,
												   alpha = 0.1,
												   color = color, zorder = - 10 ) 

								axit.plot( model_wavz, model_dered_amag, color = 'r' )
								axit.fill_between( x = model_wavz,
												   y1 = model_dered_amag - model_deredsig_amag,
												   y2 = model_dered_amag + model_deredsig_amag,
												   alpha = 0.1,
												   color = color, zorder = - 10 ) 


								for i, j in enumerate( filters ):
									color = wavcolor[j]
									axi.errorbar( x = wav[i], y = disc_amag[i], yerr = discsig_amag[i], alpha = 0.8,
												   fmt = '.', color = color, ms = emarkersize, zorder = - i )
									axi.errorbar( x = wav[i], y = dered_amag[i], yerr = deredsig_amag[i], alpha = 0.8,
												   fmt = '.', color = color, ms = emarkersize, zorder = -i )

								axi.text( 0.75, 0.2, 'z = %.2f' %z,  
										  transform = axi.transAxes, horizontalalignment='left', fontsize = 6 )    
								axi.text( 0.75, 0.13, '%s | $\chi^{2}_{N}$ = %.2f'%(dust.swapcase(),reduced_chisq),
										  transform = axi.transAxes, horizontalalignment='left', fontsize = 6  )
								axi.text( 0.75, 0.06, 'E(B-V) = %.2f$\pm$%.2f'%(ebmv,ebmvsig),
										  transform = axi.transAxes, horizontalalignment='left', fontsize = 6  )
				

								if dust in ( 'mw', 'agn' ):
									axit.xaxis.set_major_formatter( NullFormatter() )
									axit.xaxis.set_minor_formatter( NullFormatter() )
									# axi.yaxis.set_major_locator( MaxNLocator( prune = 'lower' ) )
								else:
									axit.xaxis.set_major_formatter( ScalarFormatter() )
									axit.xaxis.set_minor_formatter( ScalarFormatter() )
									# axi.yaxis.set_major_locator( MaxNLocator( prune = 'upper' ) )

								# Update plots
								fig3.subplots_adjust( left = 0.09, right = 0.99, hspace = 0, wspace = 0.25 )
								fig3.canvas.draw_idle()     


						# Update summary
						if verbose:
							print ''
							print '***** Object %i *****' %lc_id
							print 'z = %2.1f | Mg = %5.2f ' %( z, mean_mag[1] )
							print 'Removed %i epochs (%2.1f%%) - %i remaining' %( nrem, rpc, elen )
							print 'Loop succeed on %i iteration' %( loopnum )


						# Output
						if output:
							outputtable.add_row( outvals )

							# Update user
							if verbose:
								print 'Output: Updated output table for %s' %(outputfile)

						if truth_plot == 'M':
							raw_input( '...' )
							[ i.clear() for i in ( ax21, ax22, ax23, ax24, ax25, ax26 ) ]
							[ i.clear() for i in ( ax23t, ax24t ) ]
							[ i.clear() for i in ( ax31, ax32, ax33, ax34 ) ]
							[ i.clear() for i in ( ax31t, ax32t, ax33t, ax34t ) ]

					# Save output
					if output:
						outputtable.write( outputdir + outputfile, format = 'ascii' )
						
						# Update user
						if verbose:
							print 'Output: Saved output table to %s' %( outputdir + outputfile )