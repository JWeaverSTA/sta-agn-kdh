# 2x3 analysis evolution plot

# Imports
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.io import ascii
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import astro_funclib as astrotools


# Read in data
pref = 'dustout/dustout_'
sp = [ ('amg_', 'galaxy'), ('amd_', 'disk'), ('amm', 'mean' ) ]
posf = '.txt'
files = [ ['smc','cyan'], ['lmc', 'r'], ['mw','g'], ['gask', 'purple' ]  ]
wav = np.array([3543., 4770., 6231., 7625., 9134.])
emlines = { 'Ly$\\alpha$'  : 1215.24,
			'C IV'  : 1549.48,
			'C III' : 1908.734,
			'Mg II' : 2799.117,
			'H$\\alpha$'   : 6564.61 
			}



dat = ascii.read( 'Output/DustyOutput.trial.cat', format = 'fixed_width' )

# Mg vs. z, Mg vs. g-i, g-i vs. z
# Mean magnitude, disk mags, dust corr mags

labelsize = 15
textsize = 13

# Setup figure
fig1, ( ( ax12, ax13 ),
        ( ax22, ax23 ),
        ( ax32, ax33 ) ) = plt.subplots( nrows = 3, ncols = 2, figsize=(10,10) )
        

     

#ax11.set_ylabel( 'M$_{g}$', fontsize = labelsize )
ax12.set_ylabel( 'M$_{g}$', fontsize = labelsize )  
ax13.set_ylabel( '$g-i$', fontsize = labelsize )

#ax11.set_xlabel( 'z', fontsize = labelsize )
ax12.set_xlabel( '$g-i$', fontsize = labelsize )
ax13.set_xlabel( 'z', fontsize = labelsize )

#ax21.set_ylabel( 'M$_{g}$', fontsize = labelsize )
ax22.set_ylabel( 'M$_{g}$', fontsize = labelsize )  
ax23.set_ylabel( '$g-i$', fontsize = labelsize )

#ax21.set_xlabel( 'z', fontsize = labelsize )
ax23.set_xlabel( '$g-i$', fontsize = labelsize )
ax23.set_xlabel( 'z', fontsize = labelsize )

#ax31.set_ylabel( 'M$_{g}$', fontsize = labelsize )
ax32.set_ylabel( 'M$_{g}$', fontsize = labelsize )  
ax33.set_ylabel( '$g-i$', fontsize = labelsize )

#ax31.set_xlabel( 'z', fontsize = labelsize )
ax32.set_xlabel( '$g-i$', fontsize = labelsize )
ax33.set_xlabel( 'z', fontsize = labelsize )


for i in ( ax12, ax13, ax22, ax23, ax32, ax33):
    i.minorticks_on()
    i.tick_params( axis = 'both',
                      direction = 'in',
                      width = 2,
                      right = 'on',
                      length = 6,
                      color = 'black',
                      labelcolor = 'black' )
    i.tick_params( which = 'minor',
                      axis = 'both',
                      right = 'on',
                      direction = 'in',
                      color = 'black')
    nbins = len(i.get_yticklabels()) # added
    i.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))

for i in ( ax12, ax13, ax22, ax23 ):
    i.tick_params( labelbottom = 'off' )
    



fig1.subplots_adjust( hspace = 0.0, wspace = 0.25,
                      left = 0.08, right = 0.98,
                      top = 0.99, bottom = 0.05 )

colormap = plt.cm.get_cmap( 'gnuplot2')
colormap.set_under('w')


# Bouding box values
z_lo = 0
z_hi = 4
col_lo = -0.6
col_hi = 1
mag_lo = -19
mag_hi = -27

xbins = 50
ybins = 40

xedges_z = np.linspace( z_lo, z_hi, xbins )
yedges_z = np.linspace( z_lo, z_hi, ybins )

xedges_col = np.linspace( col_lo, col_hi, xbins )
yedges_col = np.linspace( col_lo, col_hi, ybins )

xedges_mag = np.linspace( mag_hi, mag_lo, xbins )
yedges_mag = np.linspace( mag_hi, mag_lo, ybins )


# Mask for chosen law

# Mask for NaN
meanmag_g_mask = dat['g_mean_amag'] == dat['g_mean_amag']
dat = dat[meanmag_g_mask]

diskmag_g_mask = dat['g_disc_amag'] == dat['g_disc_amag']
dat = dat[diskmag_g_mask]

corrmag_g_mask = dat['g_dered_smc_amag'] == dat['g_dered_smc_amag']
dat = dat[corrmag_g_mask]

meanmag_i_mask = dat['i_mean_amag'] == dat['i_mean_amag']
dat = dat[meanmag_i_mask]

diskmag_i_mask = dat['i_disc_amag'] == dat['i_disc_amag']
dat = dat[diskmag_i_mask]

corrmag_i_mask = dat['i_dered_smc_amag'] == dat['i_dered_smc_amag']
dat = dat[corrmag_i_mask]


# Grab relevant columns
redshift = dat['reds']

meanmag_g = dat['g_mean_amag']
meanmag_i = dat['i_mean_amag']
meancol_gi = meanmag_g - meanmag_i

diskmag_g = dat['g_disc_amag']
diskmag_i = dat['i_disc_amag']
diskcol_gi = diskmag_g - diskmag_i

corrmag_g = dat['g_dered_smc_amag']
corrmag_i = dat['i_dered_smc_amag']
corrcol_gi = corrmag_g - corrmag_i


# Calculate nu^{1/3} power law colours
diskmodel_col = astrotools.accretion_color( wav[1] , wav[3] )

# Create plotting function
def plot_smooth( x, y, bins, cmap, ax ):
    h, x, y, p = plt.hist2d( x, y, bins = bins )
    ax.imshow( h.T, aspect = 'auto', interpolation = 'gaussian',
               vmin = 8,
               extent = (x[0], x[-1], y[0], y[-1]), cmap = cmap )
  
# Populate row 1 - mean magnitudes
#plot_smooth( redshift, meanmag_g,
#             bins = (xedges_z, yedges_mag),
#             cmap = colormap, ax = ax11 )        
                             
plot_smooth( meancol_gi, meanmag_g,
             bins = (xedges_col, yedges_mag),
             cmap = colormap, ax = ax12 )  
             
ax12.vlines( diskmodel_col, mag_hi, mag_lo,
             colors='r',linestyles='--' )
plot_smooth( redshift, meancol_gi,
             bins = (xedges_z, yedges_col),
             cmap = colormap, ax = ax13 )
ax13.hlines( diskmodel_col, z_lo, z_hi,
             colors='r', linestyles='--' )                     
             
# Populate row 2 - disk magnitudes
#plot_smooth( redshift, diskmag_g,
#             bins = (xedges_z, yedges_mag),
#             cmap = colormap, ax = ax21 ) 
plot_smooth( diskcol_gi, diskmag_g,
             bins = (xedges_col, yedges_mag),
             cmap = colormap, ax = ax22 )
ax22.vlines( diskmodel_col, mag_hi, mag_lo,
             colors='r',linestyles='--' )
plot_smooth( redshift, diskcol_gi,
             bins = (xedges_z, yedges_col),
             cmap = colormap, ax = ax23 )
ax23.hlines( diskmodel_col, z_lo, z_hi,
             colors='r', linestyles='--' )

# Populate row 3 - corrected magnitudes
#plot_smooth( redshift, corrmag_g,
#             bins = (xedges_z, yedges_mag),
#             cmap = colormap, ax = ax31) 
plot_smooth( corrcol_gi, corrmag_g,
             bins = (xedges_col, yedges_mag),
             cmap = colormap, ax = ax32 )
ax32.vlines( diskmodel_col, mag_hi, mag_lo,
            colors='r',linestyles='--' )
plot_smooth( redshift, corrcol_gi,
             bins = (xedges_z, yedges_col),
             cmap = colormap, ax = ax33 )
legp = ax33.hlines( diskmodel_col, z_lo, z_hi,
             colors='r', linestyles='--' )


# Histrograms middle column

ax12tx = plt.axes([0,0,1,1])
ip1x = InsetPosition( ax12, [0,0,1,0.1] )
ax12tx.set_axes_locator(ip1x)
ax12tx.axis('off')
ax12tx.hist( meancol_gi, color = '#3224d4', alpha = 1, histtype='step')

ax12ty = plt.axes([0,0,.999,1])
ip1y = InsetPosition( ax12, [0.9,0,0.1,1] )
ax12ty.set_axes_locator(ip1y)
ax12ty.axis('off')
ax12ty.hist( meanmag_g, color = '#3224d4', alpha = 1, histtype='step',
             orientation = 'horizontal')
ax12ty.invert_xaxis()
             
ax22tx = plt.axes([0,0,0.998,1])
ip2x = InsetPosition( ax22, [0,0,1,0.1] )
ax22tx.set_axes_locator(ip2x)
ax22tx.axis('off')
ax22tx.hist( diskcol_gi, color = '#3224d4', alpha = 1, histtype='step')

ax22ty = plt.axes([0,0,0.9999,1])
ip2y = InsetPosition( ax22, [0.9,0,0.1,1] )
ax22ty.set_axes_locator(ip2y)
ax22ty.axis('off')

ax22ty.hist( diskmag_g, color = '#3224d4', alpha = 1, histtype='step',
             orientation = 'horizontal')
ax22ty.invert_xaxis()

ax32tx = plt.axes([0,0,0.999998,1])
ip3x = InsetPosition( ax32, [0,0,1,0.1] )
ax32tx.set_axes_locator(ip3x)
ax32tx.axis('off')
ax32tx.hist( corrcol_gi, color = '#3224d4', alpha = 1, histtype='step')

ax32ty = plt.axes([0,0,0.99,1])
ip3y = InsetPosition( ax32, [0.9,0,0.1,1] )
ax32ty.set_axes_locator(ip3y)
ax32ty.axis('off')
ax32ty.hist( corrmag_g, color = '#3224d4', alpha = 1, histtype='step',
             orientation = 'horizontal')
ax32ty.invert_xaxis()

ax33tx = plt.axes([0,0,0.99998,1])
ip3ax = InsetPosition( ax33, [0,0,1,0.1] )
ax33tx.set_axes_locator(ip3ax)
ax33tx.axis('off')
ax33tx.hist( redshift, color = '#3224d4', alpha = 1, histtype='step')
            
            
# Indicates emission line peaks
# g - band
wav_g = wav[1]
for em in emlines.keys():
	emwav = emlines[em]
	em_z = ( wav_g / emwav ) - 1.
	if em_z < 0:
		continue
	legg = ax13.vlines( em_z, col_hi, col_hi - 0.80 * (col_hi-col_lo),
            colors = 'g', alpha = 1, linestyles = 'dotted' )
	ax13.text( em_z, col_hi - 0.83 * ( col_hi - col_lo ),
           	em, fontsize = 10, alpha = 0.9, horizontalalignment='center', color = 'k')

# i - band             
wav_i = wav[3]
for em in emlines.keys():
	if em == 'C IV':
		continue
	emwav = emlines[em]
	em_z = ( wav_i / emwav ) - 1.
	legi = ax13.vlines( em_z, col_hi, col_hi - 0.90 * (col_hi-col_lo),
            colors = 'orange', alpha = 1, linestyles = 'dotted' )
	ax13.text( em_z, col_hi - 0.93 * ( col_hi - col_lo ),
           	em, fontsize = 10, alpha = 0.9, horizontalalignment='center', color = 'k')   

        
fig1.legend( (legp, legg, legi),
			 ('F$_{\\nu}$ ~ $\\nu^{1/3}$', '$g$ centre', '$i$ centre'),
			  bbox_to_anchor = (0.97, 0.06), loc = 'lower right' )

# Labels
ax12.text(0.4, 0.9, 'Mean Spectrum', fontsize = textsize, transform=ax12.transAxes)
ax22.text(0.4, 0.9, 'Disk Spectrum', fontsize = textsize, transform=ax22.transAxes)
ax32.text(0.4, 0.9, 'Dust Corrected (SMC)', fontsize = textsize, transform=ax32.transAxes)

             
# Tweak bounding boxes

#ax11.set_xlim( z_lo, z_hi )
ax12.set_xlim( col_lo, col_hi )
ax13.set_xlim( z_lo, z_hi )

#ax21.set_xlim( z_lo, z_hi )
ax22.set_xlim( col_lo, col_hi )
ax23.set_xlim( z_lo, z_hi )

#ax31.set_xlim( z_lo, z_hi )
ax32.set_xlim( col_lo, col_hi )
ax33.set_xlim( z_lo, z_hi )

#ax11.set_ylim( mag_lo, mag_hi )
ax12.set_ylim( mag_lo, mag_hi )
ax13.set_ylim( col_hi, col_lo )

#ax21.set_ylim( mag_lo, mag_hi )
ax22.set_ylim( mag_lo, mag_hi )
ax23.set_ylim( col_hi, col_lo )

#ax31.set_ylim( mag_lo, mag_hi )
ax32.set_ylim( mag_lo, mag_hi )
ax33.set_ylim( col_hi, col_lo )


ax12tx.set_xlim( col_lo, col_hi )
ax22tx.set_xlim( col_lo, col_hi )
ax32tx.set_xlim( col_lo, col_hi )

ax12ty.set_ylim( mag_lo, mag_hi )
ax22ty.set_ylim( mag_lo, mag_hi )
ax32ty.set_ylim( mag_lo, mag_hi )

