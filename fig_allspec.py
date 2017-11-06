# Mother of all spectra

# NOTES: add colorbar ?
#		 normalisation technique?

# Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import AutoMinorLocator
from astropy.io import ascii
from astropy.table import Column
import adv_funclib as advtools
import astro_funclib as astrotools



# Parameters
dustlaw = 'smc'
wavdict = { 'u' : 3543.,
		'g' : 4770.,
		'r' : 6231.,
		'i' : 7625.,
		'z' : 9134. 
		}
filt = ( 'u', 'g', 'r', 'i', 'z')
wav = np.array( [ 3543., 4770., 6231., 7625., 9134. ] )
emlines = { 'Ly$\\alpha$'  : 1215.24,
			'C IV'  : 1549.48,
			'C III' : 1908.734,
			'Mg II' : 2799.117,
			'H$\\alpha$'   : 6564.61 
			}

# Read in data
dat = ascii.read( 'Output/DustyOutput.trial.cat', format = 'fixed_width' )
#dat = dat_raw[ dat_raw['dustname'] == dustlaw ]

z = dat['reds']
dalen = len( z )

sel = 'lred_'
sigsel = 'lredsig_'

outputfilename = 'Figures/motherplot_smc.png'
dust = 'smc'

pname = 'Dust Corrected Disc Specta | smc'
#---- TO GRAB NON LUM DAT (This is crappy coding...)

#rsel = 'amm_'
#rsel = 'amd_'
#rsel = ''

# RESIDUALS NOT CORRECT and NEED TO USE mult factor!

# Mult factor
datf = ascii.read('composite_output.txt', format = 'fixed_width')
sclf = np.e**(datf['scl'])

for i in filt:
	cdat, cdat_err = astrotools.amag2lum( dat[i + '_dered_'+dust+'_amag'], dat[i+'_dered_'+dust+'_amerr'] )
	cdatacc, _  = astrotools.amag2lum( astrotools.accretion_magspec( wavdict[i] / (1+z), 2400. ), 0. )
	cdatr = cdat - cdatacc
	ccol = Column( cdat, name = i+'_dered_lum' )
	ccolsig = Column( cdat_err, name =i+'_dered_lerr' )
	cres = Column( cdatr, name = i+'_residual_lum' )
	dat.add_column( ccol )
	dat.add_column( ccolsig )
	dat.add_column( cres )


#----
plt.ioff()

lum = np.array( [ [ dat[ i + '_dered_lum'][j] for i in filt ] for j in range( dalen ) ] )
lumsig = np.array( [ [ dat[i + '_dered_lerr'][j] for i in filt ] for j in range( dalen ) ] )

print 'Data loaded - length: ', dalen
dlen = dalen

# Mask out nan
'''
nanmask = ( lum == lum )
nlost = np.sum(~nanmask) / 5
nkept = np.sum(nanmask) / 5
lum = lum[nanmask].reshape( nkept, 5 )
z = z[np.sum(nanmask, axis = 1 ) / 5 == 1]

dlen = nkept
# NOT WORKING - NEED 2D -> 2D

print 'NaN masked - removed: ', nlost
'''
# Setup figure
fig1 = plt.figure( figsize = (15,10), dpi = 300 )
ax11 = fig1.add_subplot(111)
ax11.set_xlabel('Rest Wavelength ($\AA$)', fontsize = 15)
ax11.set_ylabel('Object', fontsize = 15)

xlo = wavdict['u']/( 1 + max( z ) ) 
xhi = wavdict['z']/( 1 + min( z ) )
ax11.set_xlim( xlo, xhi )
ax11.set_ylim( 0, dlen )

# Setup grid
inds = z.argsort()
zsort = z[inds]
lumsort = lum[inds]
lumsigsort = lumsig[inds]
sclf = sclf[inds]
xbins = 2000
ybins = dlen
wavgrid = np.linspace( xlo, xhi, xbins )
imagegrid = np.zeros( shape = ( xbins, ybins ) )

for i in xrange( dlen ):
	wavz = astrotools.obs2restwav( wav, zsort[i] )
	intspec = ( [ advtools.linint( wavz, lumsort[i], lumsigsort[i], j, extrapolation = False )[0] for j in wavgrid ] )
	norm_factor = sclf[i]
	norm_intspec = intspec / norm_factor
	#norm_intspec[norm_intspec!=norm_intspec] = 0.
	imagegrid[:,i] = np.log10(norm_intspec)
	print '%3.1f %%'%(100.*i/float(dlen))

# Plot
cmap = plt.cm.gnuplot2
sig = np.nanstd( imagegrid )
mean = np.nanmean( imagegrid )
vmin = mean - 2*sig
vmax = mean + 2*sig
ax11.imshow( imagegrid.T, cmap = cmap,
			 vmin = vmin, vmax = vmax,
			 extent = [xlo, xhi, 0, dlen],
			 origin = 'lower' )

'''
# Interpolate spectra
Nbins = 30
inds = z.argsort()
zsort = z[inds]
lumsort = lum[inds]
px = np.zeros( Nbins * dlen )
py = np.zeros( Nbins * dlen )
pz = np.zeros( Nbins * dlen ) * np.nan
dw = np.zeros( Nbins * dlen )
dh = np.zeros( Nbins * dlen )
for i in xrange( dlen ):
	# setup grid
	wavmin = wavdict['u'] / (1 + zsort[i])
	wavmax = wavdict['z'] / (1 + zsort[i])
	
	
	wavz = restwav( wav, zsort[i])

	# reinterpolate
	intspec = ( [ linint( wavz, lumsort[i], j, extrapolation = False ) for j in wavgrid ] )
	norm_factor = np.nanmax(intspec)
	norm_intspec = intspec / norm_factor

	px[i*Nbins:(i*Nbins)+Nbins] = wavgrid
	py[i*Nbins:(i*Nbins)+Nbins] = (i+1) * np.ones( Nbins )
	pz[i*Nbins:(i*Nbins)+Nbins] = norm_intspec
	
	rangewav = wavgrid[-1] - wavgrid[0]
	dw[i*Nbins:(i*Nbins)+Nbins] = rangewav / float( Nbins ) * np.ones( Nbins )
	dh[i*Nbins:(i*Nbins)+Nbins] = np.ones( Nbins )

# Plot
cmap = plt.cm.gnuplot2
checker = 0
for x, y, c, w, h in zip( px, py, pz, dw, dh ):
	ax11.add_artist( Rectangle( xy = (x,y), color = cmap(c**2), 
					  width = w, height = 0.9*h, edgecolor = None, rasterized = True ) )
	checker += 1
	print '%3.1f'%(100* checker / float(len(px)))
'''	

# More plotting
ax11.vlines( emlines.values(), 0, dlen,
			 color = 'k', alpha = 0.6, linestyles = 'dotted' )
for i in range(len(emlines)):
	ax11.text( emlines.values()[i], 1.02*dlen, emlines.keys()[i],
			   fontsize = 10, alpha = 1, horizontalalignment='center' )

zloc = np.array( (1.,2.,3.) )
znum = [ advtools.linint( zsort, np.arange(dlen), 0. * zsort, zloc[i] )[0] for i in xrange( len( zloc ) ) ]	
znum = np.round( znum )

for i in range( len( znum ) ):
	ax11.hlines( znum[i], xlo, xlo + 0.91 * (xhi - xlo),
				 color = 'k', linestyles = 'dotted', alpha = 0.6)
	ax11.text( xlo + 0.95 * (xhi - xlo), znum[i], '$z=%i$'%zloc[i], 
				fontsize = 10, verticalalignment = 'center', horizontalalignment = 'center' )

fig1.subplots_adjust( left = 0.05, right = 0.97, bottom = 0.05, top = 0.9 )

ax11.minorticks_on()
ax11.tick_params( axis = 'both',
                  direction = 'in',
                  width = 2,
                  right = 'on',
                  length = 6 )
ax11.tick_params( which = 'minor',
                  axis = 'both',
                  right = 'on',
                  direction = 'in',)
fig1.suptitle( pname + ' | N = ' + str(dlen), y = 0.98, fontsize = 15 )

ax11.set_xscale('log')
ax11.xaxis.set_ticks( [1000,2000,3000,4000,5000,6000,7000] )
ax11.xaxis.set_ticklabels( ['%i'%i for i in ax11.xaxis.get_majorticklocs() ] )
minor_locator = AutoMinorLocator(5)
ax11.xaxis.set_minor_locator(minor_locator)
ax11.tick_params( which = 'minor', labelbottom = 'off')
ax11.set_xlim(xlo,xhi)
fig1.savefig(outputfilename)