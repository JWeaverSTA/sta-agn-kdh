# Composite spectrum

# Imports
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import adv_funclib as advtools
import astro_funclib as astrotools
from matplotlib.ticker import ScalarFormatter, NullFormatter

showplot = False


# Parameters
dustlaw = 'smc' # make multiple composites?
fidwav = 2400.
wavdict = { 'u' : 3543.,
		'g' : 4770.,
		'r' : 6231.,
		'i' : 7625.,
		'z' : 9134. 
		}
filt = ( 'u', 'g', 'r', 'i', 'z')
wav = np.array( [ 3543., 4770., 6231., 7625., 9134. ] )
emlines = { 'Ly$_{\infty}$'  : 912.00,
		 'Ly$_{\\alpha}$' : 1215.67,
		 'CIV'            : 1549.05,
		 '[CIII]'         : 1908.73,
		 '[MgII]'         : 2797.92,
		 'H$_{\infty}$'   : 3560.00,
		 'H$_{\\beta}$'   : 4861.32}


N = 1

# Read in data
dat = ascii.read( 'Output/DustyOutput.trial.cat', format = 'fixed_width' )

z = dat['reds']
dalen = len( z )

wavz = np.array( [ astrotools.obs2restwav( wav, z[i] ) for i in xrange( len( z ) ) ] )
dered_amag = np.array( [ [ dat[j][i+'_dered_'+dustlaw+'_amag'] for i in filt ] for j in xrange(dalen) ] )
dered_amerr = np.array( [ [ dat[j][i+'_dered_'+dustlaw+'_amerr'] for i in filt ] for j in xrange(dalen) ] )

lum, lumsig = astrotools.amag2lum( dered_amag, dered_amerr )

print 'Data loaded - length: ', dalen
dlen = dalen

# Mask out nan
'''
nanmask = np.sum( lum == lum, 1 )
nlost = np.sum(~nanmask) / 5
nkept = np.sum(nanmask) / 5
lum = lum[nanmask].reshape( nkept, 5 )
lumsig = lumsig[nanmask].reshape(nkept,5)
datmask = np.sum(nanmask, axis = 1 ) / 5 == 1
wavz = wavz[datmask]
dat = dat[datmask]
dlen = nkept
'''
#print 'NaN masked - removed: ', nlost


# FAKE DATA - CANNOT DO THIS! accpowerlaw is in MAGNITUDES!
'''
print 'FAKE DATA'
np.random.seed(1)
fidwav = 2400.
fakelo = np.nanmin(lum)
fakehi = np.nanmax(lum)
randlum = fakelo + np.random.random( dlen ) * (fakehi - fakelo)
lum = np.array( [ [ mag2lum(accpowerlaw( wavz[i,j], fidwav ) ) + randlum[i] for j in range(len(filt)) ] for i in xrange(dlen) ] )
lumsig = np.ones(shape=lum.shape)
'''

# Setup wavelength grid  (should probs add .5 of bin)
Nbins = 1000
wavmin = wavdict['u'] / ( 1. + max( z ) ) 
wavmax = wavdict['z'] / ( 1. + min( z ) ) 
wavgrid = np.logspace( np.log10(wavmin), np.log10(wavmax), Nbins )

# Plotting
if showplot:
	fig1, ( ( ax11), (ax12) ) = plt.subplots( nrows = 2, ncols = 1,
	                                gridspec_kw = {'height_ratios':[4,1] } )



# First Setup
f = np.ones( dlen )
fsig = np.ones( dlen )
Lmod_avg = np.zeros( Nbins )
Lmod_sig = np.zeros( Nbins )
Lmod_psig = np.zeros( Nbins )
Nbin_contrib = np.zeros( Nbins )

# loop over grid
for i in xrange( Nbins ):
	hold_intavg = np.zeros( dlen )
	hold_intsig = np.zeros( dlen )

	# for each grid bin, optimally average interpolated value
	for j in xrange( dlen ):
		hold_intavg[j], hold_intsig[j] = advtools.linint( wavz[j], lum[j], lumsig[j], wavgrid[i], extrapolation = False)

	Nbin_contrib[i] = np.sum( ~ np.isnan( hold_intavg ) )
	Lmod_avg[i], Lmod_sig[i] = advtools.optavg( hold_intavg, hold_intsig )

ylo = 1.3 * np.nanmin(Lmod_avg)
yhi = 1.3 * np.nanmax(Lmod_avg)
xlo = wavgrid[0]
xhi = wavgrid[-1]	

if showplot:
	plum = np.array( [ lum[m]/f[m] for m in range(len(lum)) ] )
	xbins = np.logspace( np.log10(xlo), np.log10(xhi), 100)
	ybins = np.logspace( np.log10(ylo), np.log10(yhi), 50)
	H, xedges, yedges = np.histogram2d( wavz.flatten(), plum.flatten(), 
				 bins = [xbins, ybins],
				 range = ( (xlo, xhi), (ylo, yhi) ) )
	ax11.pcolormesh( xedges, yedges, H.T,
				 cmap = plt.cm.Reds )
	#xbins = np.linspace( wavmin, wavmax, 100 )
	#ybins = np.linspace( np.nanmin(plum), np.nanmax(plum), 100 )
	#counts, _, _ = np.histogram2d( wavz.flatten(), plum.flatten(),
	#                               bins=(xbins, ybins) )
	#ax11.pcolormesh( xbins, ybins, counts.T, cmap = 'Greys' )

	ax11.plot( wavgrid, Lmod_avg, 'navy' , zorder = 10)
	ax11.fill_between( wavgrid, Lmod_avg - Lmod_sig, Lmod_avg + Lmod_sig, alpha = 0.2, color = 'royalblue', zorder = 10)
	ax11.set_yscale('log')
	ax11.set_xscale('log')
	ax11.set_ylim( 1.3*np.nanmin(Lmod_avg), 1.3*np.nanmax(Lmod_avg) )

# loop over objects
count = 0
#figt,axt = plt.subplots()
fdiff = 99.
ln_f_sig = 1.
while fdiff > 1E-3: #np.mean( ln_f_sig ):

	# Save old f
	f_old = f
	ln_f_old = np.log(f_old)


	# Calculate f
	for k in xrange( dlen ):
		x, sig = lum[k], lumsig[k]
		wz = wavz[k]
		xlen = len( x )
		Lmod_int = np.zeros( xlen )
		for l in xrange( xlen ):
			Lmod_int[l], _ = advtools.linint( wavgrid, Lmod_avg, Lmod_sig, wz[l] )
		f[k], fsig[k] = advtools.optscl( x, sig, Lmod_int )
		

	# Impose <log(f)> = 0
	ln_f_avg, ln_f_sig = advtools.optavg( np.log(f), ( fsig / f ) )

	ln_f= np.log(f) - ln_f_avg
	f = np.exp( ln_f )

	ylo = 1.3 * np.nanmin(Lmod_avg)
	yhi = 1.3 * np.nanmax(Lmod_avg)
	xlo = wavgrid[0]
	xhi = wavgrid[-1]
	

	if showplot:
		ax11.clear()
		ax11.set_yscale('log')
		ax11.set_xscale('log')
		ax11.set_xlim( xlo, xhi)
		ax11.set_ylim( ylo, yhi )
		plum = np.array( [ lum[m]/f[m] for m in range(len(lum)) ] )
		plumsig = np.array( [ lumsig[m]/f[m] for m in range( len(lum )) ] )
		xbins = np.logspace( np.log10(xlo), np.log10(xhi), 100)
		ybins = np.logspace( np.log10(ylo), np.log10(yhi), 50)
		H, xedges, yedges = np.histogram2d( wavz.flatten(), plum.flatten(), 
				 bins = [xbins, ybins],
				 range = ( (xlo, xhi), (ylo, yhi) ) )
		ax11.pcolormesh( xedges, yedges, H.T,
				 cmap = plt.cm.Reds )
			

	# recalculate L
	for i in xrange( Nbins ):
		hold_intavg = np.zeros( dlen )
		hold_intsig = np.zeros( dlen )

		# for each grid bin, optimally average interpolated value
		for j in xrange( dlen ):
			hold_intavg[j], hold_intsig[j] = advtools.linint( wavz[j], lum[j] / f[j], lumsig[j] / f[j], wavgrid[i], extrapolation = False)

		Nbin_contrib[i] = np.sum( ~ np.isnan( hold_intavg ) )
		Lmod_avg[i], Lmod_sig[i] = advtools.optavg( hold_intavg, hold_intsig )
		Lmod_psig[i] = np.sqrt( Lmod_sig[i]**2 + advtools.std(hold_intavg, Lmod_avg[i], hold_intsig )**2 )

	if showplot:
		ax11.plot( wavgrid, Lmod_avg, 'navy' , zorder = 10)
		ax11.fill_between( wavgrid, Lmod_avg - Lmod_sig, Lmod_avg + Lmod_sig, alpha = 0.2, color = 'royalblue', zorder = 10)
		plt.pause(0.001)


	fdiff = abs( ln_f.mean() - ln_f_old.mean() )	
	count += 1
	print count, '||', fdiff, '||', 1E-2 *np.mean( ln_f_sig )
	#axt.scatter( count, fdiff, c='k' )
	#axt.scatter( count, 1E-8 * np.mean( ln_f_sig), c='r')
	#axt.set_yscale('log')

plt.close('all')

fig1, ( ( ax11), (ax12) ) = plt.subplots( nrows = 2, ncols = 1, figsize = (14,10),
                                gridspec_kw = {'height_ratios':[4,1] } )	
ax11.set_yscale('log')
ax11.set_xscale('log')
ylo = 1.3 * np.nanmin(Lmod_avg)
yhi = 1.3 * np.nanmax(Lmod_avg)
xlo = wavgrid[0]
xhi = wavgrid[-1]
ax11.set_ylim( ylo, yhi )
ax11.set_xlim( xlo, xhi )
plum = np.array( [ lum[m]/f[m] for m in range(len(lum)) ] )
plumsig = np.array( [ lumsig[m]/f[m] for m in range( len(lum )) ] )
xbins = np.logspace( np.log10(xlo), np.log10(xhi), 100)
ybins = np.logspace( np.log10(ylo), np.log10(yhi), 50)
H, xedges, yedges = np.histogram2d( wavz.flatten(), plum.flatten(), 
			 bins = [xbins, ybins], normed = True, weights = 1./plumsig.flatten(),
			 range = ( (xlo, xhi), (ylo, yhi) ) )
ax11.pcolormesh( xedges, yedges, H.T,
			 cmap = plt.cm.Reds )

fidwav = 2400.
wcalc = abs( wavgrid - fidwav )
Lfid = Lmod_avg[ wcalc == min( wcalc ) ]
fidmag = astrotools.lum2amag( Lfid, 0. )[0]
accmag = astrotools.accretion_magspec( wavgrid, fidwav ) + fidmag
lumt, _ = astrotools.amag2lum( accmag, 0. )
ax11.plot( wavgrid, lumt.flatten(), c='k', ls='dotted' )
ax11.plot( wavgrid, Lmod_avg, 'navy' , zorder = 10)
ax11.fill_between( wavgrid, Lmod_avg - Lmod_psig, Lmod_avg + Lmod_psig,
                   alpha = 0.2, color = 'royalblue', zorder = 10)

ax11.vlines( emlines.values(), ylo, ylo + 0.55 * (yhi-ylo),
		 color = 'k', alpha = 0.6, linestyles = 'dotted' )
for i in range(len(emlines)):
	ax11.text( emlines.values()[i], ylo + 0.65 * (yhi-ylo), emlines.keys()[i],
			   fontsize = 10, alpha = 1, horizontalalignment='center' )

ax12.set_yscale('log')
ax12.set_xscale('log')
ax12.set_xlim( xlo, xhi )
ax12.plot( wavgrid, Nbin_contrib, 'navy')

ax11.xaxis.set_minor_formatter( NullFormatter() )
ax12.xaxis.set_minor_formatter( ScalarFormatter() )

ax12.set_xlabel( '$\lambda$ = $\lambda_{o}$(1 + z) ($\AA$)' )
ax11.set_ylabel( 'Luminosity (erg/s)' )
ax12.set_ylabel( 'N/bin')
fig1.subplots_adjust( hspace = 0. )


output = Table( data = [dat['dbID'], ln_f], names = ['dbID', 'scl'] )
ascii.write( output, 'composite_output.txt', format = 'fixed_width')
