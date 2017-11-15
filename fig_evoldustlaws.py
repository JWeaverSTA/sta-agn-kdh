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
import adv_funclib as advtools


# Read in data
pref = 'dustout/dustout_'
sp = [ ('amg_', 'galaxy'), ('amd_', 'disk'), ('amm', 'mean' ) ]
posf = '.txt'
files = [ ['smc','cyan'], ['lmc', 'r'], ['mw','g'], ['gask', 'purple' ]  ]
wav = np.array([3543., 4770., 6231., 7625., 9134.])
emlines = { 'Ly$_{\infty}$'  : 912.00,
	 'Ly$_{\\alpha}$' : 1215.67,
	 'CIV'            : 1549.05,
	 '[CIII]'         : 1908.73,
	 '[MgII]'         : 2797.92,
	 'H$_{\infty}$'   : 3560.00,
	 'H$_{\\beta}$'   : 4861.32,
	 'H$_{\\alpha}$'  : 6562.80}


pw1 = 2
pw2 = 3
pwr = pw1 / float(pw2)

dat = ascii.read( 'Output/DustyOutput.powerlaw.%i.%i.cat'%(pw1,pw2), format = 'fixed_width' )
dlaw = 'agn'
# Mg vs. z, Mg vs. g-i, g-i vs. z
# Mean magnitude, disk mags, dust corr mags

labelsize = 15
textsize = 13

# Setup figure
fig1, ( ( ax12, ax13 ),
		( ax22, ax23 ),
		( ax32, ax33 ),
		( ax42, ax43) ) = plt.subplots( nrows = 4, ncols = 2, figsize=(10,10) )
		

	 

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

#ax31.set_ylabel( 'M$_{g}$', fontsize = labelsize )
ax42.set_ylabel( 'M$_{g}$', fontsize = labelsize )  
ax43.set_ylabel( '$g-i$', fontsize = labelsize )

#ax31.set_xlabel( 'z', fontsize = labelsize )
ax42.set_xlabel( '$g-i$', fontsize = labelsize )
ax43.set_xlabel( 'z', fontsize = labelsize )


for i in ( ax12, ax13, ax22, ax23, ax32, ax33, ax42, ax43 ):

	nbins = len(i.get_yticklabels()) # added
	i.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))

for i in ( ax12, ax13, ax22, ax23 ):
	i.tick_params( labelbottom = 'off' )
	



fig1.subplots_adjust( hspace = 0.0, wspace = 0.25,
					  left = 0.08, right = 0.98,
					  top = 0.93, bottom = 0.05 )

colormap = plt.cm.magma_r
colormap.set_under('w')


# Bouding box values
z_lo = 0
z_hi = 4.
col_lo = -0.8
col_hi = 0.8
mag_lo = -19
mag_hi = -29

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
corrmag_g_mask = dat['g_dered_smc_amag'] == dat['g_dered_smc_amag']
dat_smc = dat[corrmag_g_mask]

corrmag_i_mask = dat_smc['i_dered_smc_amag'] == dat_smc['i_dered_smc_amag']
dat_smc = dat_smc[corrmag_i_mask]

corrmag_g_mask = dat['g_dered_lmc_amag'] == dat['g_dered_lmc_amag']
dat_lmc = dat[corrmag_g_mask]

corrmag_i_mask = dat_lmc['i_dered_lmc_amag'] == dat_lmc['i_dered_lmc_amag']
dat_lmc = dat_lmc[corrmag_i_mask]

corrmag_g_mask = dat['g_dered_mw_amag'] == dat['g_dered_mw_amag']
dat_mw = dat[corrmag_g_mask]

corrmag_i_mask = dat_mw['i_dered_mw_amag'] == dat_mw['i_dered_mw_amag']
dat_mw = dat_mw[corrmag_i_mask]

corrmag_g_mask = dat['g_dered_agn_amag'] == dat['g_dered_agn_amag']
dat_agn = dat[corrmag_g_mask]

corrmag_i_mask = dat_agn['i_dered_agn_amag'] == dat_agn['i_dered_agn_amag']
dat_agn = dat_agn[corrmag_i_mask]


# Grab relevant columns
redshift_smc = dat_smc['reds']
redshift_lmc = dat_lmc['reds']
redshift_mw = dat_mw['reds']
redshift_agn = dat_agn['reds']

corrmag_g_smc = dat_smc['g_dered_smc_amag']
corrmag_i_smc = dat_smc['i_dered_smc_amag']
corrcol_gi_smc = corrmag_g_smc - corrmag_i_smc

corrmag_g_lmc = dat_lmc['g_dered_lmc_amag']
corrmag_i_lmc = dat_lmc['i_dered_lmc_amag']
corrcol_gi_lmc = corrmag_g_lmc - corrmag_i_lmc

corrmag_g_mw = dat_mw['g_dered_mw_amag']
corrmag_i_mw = dat_mw['i_dered_mw_amag']
corrcol_gi_mw = corrmag_g_mw - corrmag_i_mw


corrmag_g_agn = dat_agn['g_dered_agn_amag']
corrmag_i_agn = dat_agn['i_dered_agn_amag']
corrcol_gi_agn = corrmag_g_agn - corrmag_i_agn


# Calculate nu^{1/3} power law colours
diskmodel_col = astrotools.accretion_color( wav[1] , wav[3], powerlaw = pwr )

# Create plotting function
def plot_smooth( x, y, bins, cmap, ax ):
	h, x, y, p = plt.hist2d( x, y, bins = bins )
	ax.imshow( h.T, aspect = 'auto', interpolation = 'gaussian',
			   vmin = 2,
			   extent = (x[0], x[-1], y[0], y[-1]), cmap = cmap )
  


# Populate row 3 - corrected magnitudes
#plot_smooth( redshift, corrmag_g,
#             bins = (xedges_z, yedges_mag),
#             cmap = colormap, ax = ax31) 
plot_smooth( corrcol_gi_smc, corrmag_g_smc,
			 bins = (xedges_col, yedges_mag),
			 cmap = colormap, ax = ax12 )
ax12.vlines( diskmodel_col, mag_hi, mag_lo,
			colors='r',linestyles='--' )
plot_smooth( redshift_smc, corrcol_gi_smc,
			 bins = (xedges_z, yedges_col),
			 cmap = colormap, ax = ax13 )
legp = ax13.hlines( diskmodel_col, z_lo, z_hi,
			 colors='r', linestyles='--' )

plot_smooth( corrcol_gi_lmc, corrmag_g_lmc,
			 bins = (xedges_col, yedges_mag),
			 cmap = colormap, ax = ax22 )
ax22.vlines( diskmodel_col, mag_hi, mag_lo,
			colors='r',linestyles='--' )
plot_smooth( redshift_lmc, corrcol_gi_lmc,
			 bins = (xedges_z, yedges_col),
			 cmap = colormap, ax = ax23 )
legp = ax23.hlines( diskmodel_col, z_lo, z_hi,
			 colors='r', linestyles='--' )

plot_smooth( corrcol_gi_mw, corrmag_g_mw,
			 bins = (xedges_col, yedges_mag),
			 cmap = colormap, ax = ax32 )
ax32.vlines( diskmodel_col, mag_hi, mag_lo,
			colors='r',linestyles='--' )
plot_smooth( redshift_mw, corrcol_gi_mw,
			 bins = (xedges_z, yedges_col),
			 cmap = colormap, ax = ax33 )
legp = ax33.hlines( diskmodel_col, z_lo, z_hi,
			 colors='r', linestyles='--' )

plot_smooth( corrcol_gi_agn, corrmag_g_agn,
			 bins = (xedges_col, yedges_mag),
			 cmap = colormap, ax = ax42 )
ax42.vlines( diskmodel_col, mag_hi, mag_lo,
			colors='r',linestyles='--' )
plot_smooth( redshift_agn, corrcol_gi_agn,
			 bins = (xedges_z, yedges_col),
			 cmap = colormap, ax = ax43 )
legp = ax43.hlines( diskmodel_col, z_lo, z_hi,
			 colors='r', linestyles='--' )

# Histrograms middle column

ax12tx = plt.axes([0,0,0.99987,1])
ip1x = InsetPosition( ax12, [0,0,1,0.1] )
ax12tx.set_axes_locator(ip1x)
ax12tx.axis('off')
ax12tx.hist( corrcol_gi_smc, color = 'purple', alpha = 0.2, histtype='stepfilled')

ax12ty = plt.axes([0,0,0.99979,1])
ip1y = InsetPosition( ax12, [0.9,0,0.1,1] )
ax12ty.set_axes_locator(ip1y)
ax12ty.axis('off')
ax12ty.hist( corrmag_g_smc, color = 'purple', alpha = 0.2, histtype='stepfilled',
			 orientation = 'horizontal')
ax12ty.invert_xaxis()

ax13tx = plt.axes([0,0,0.99997,1])
ip1ax = InsetPosition( ax13, [0,0,1,0.1] )
ax13tx.set_axes_locator(ip1ax)
ax13tx.axis('off')
ax13tx.hist( redshift_smc, color = 'purple', alpha = 0.2, histtype='stepfilled')


ax22tx = plt.axes([0,0,0.99999,1])
pi2x = InsetPosition( ax22, [0,0,1,0.1] )
ax22tx.set_axes_locator(pi2x)
ax22tx.axis('off')
ax22tx.hist( corrcol_gi_lmc, color = 'purple', alpha = 0.2, histtype='stepfilled')

ax22ty = plt.axes([0,0,0.999999,1])
ip2y = InsetPosition( ax22, [0.9,0,0.1,1] )
ax22ty.set_axes_locator(ip2y)
ax22ty.axis('off')
ax22ty.hist( corrmag_g_lmc, color = 'purple', alpha = 0.2, histtype='stepfilled',
			 orientation = 'horizontal')
ax22ty.invert_xaxis()

ax23tx = plt.axes([0,0,0.99998,1])
ip2ax = InsetPosition( ax23, [0,0,1,0.1] )
ax23tx.set_axes_locator(ip2ax)
ax23tx.axis('off')
ax23tx.hist( redshift_lmc, color = 'purple', alpha = 0.2, histtype='stepfilled')


ax32tx = plt.axes([0,0,0.9999998,1])
ip3x = InsetPosition( ax32, [0,0,1,0.1] )
ax32tx.set_axes_locator(ip3x)
ax32tx.axis('off')
ax32tx.hist( corrcol_gi_mw, color = 'purple', alpha = 0.2, histtype='stepfilled')

ax32ty = plt.axes([0,0,0.99888,1])
ip3y = InsetPosition( ax32, [0.9,0,0.1,1] )
ax32ty.set_axes_locator(ip3y)
ax32ty.axis('off')
ax32ty.hist( corrmag_g_mw, color = 'purple', alpha = 0.2, histtype='stepfilled',
			 orientation = 'horizontal')
ax32ty.invert_xaxis()

ax33tx = plt.axes([0,0,0.998,1])
ip3ax = InsetPosition( ax33, [0,0,1,0.1] )
ax33tx.set_axes_locator(ip3ax)
ax33tx.axis('off')
ax33tx.hist( redshift_mw, color = 'purple', alpha = 0.2, histtype='stepfilled')


ax42tx = plt.axes([0,0,0.9999898,1])
ip4x = InsetPosition( ax42, [0,0,1,0.1] )
ax42tx.set_axes_locator(ip4x)
ax42tx.axis('off')
ax42tx.hist( corrcol_gi_agn, color = 'purple', alpha = 0.2, histtype='stepfilled')

ax42ty = plt.axes([0,0,0.999888,1])
ip4y = InsetPosition( ax42, [0.9,0,0.1,1] )
ax42ty.set_axes_locator(ip4y)
ax42ty.axis('off')
ax42ty.hist( corrmag_g_agn, color = 'purple', alpha = 0.2, histtype='stepfilled',
			 orientation = 'horizontal')
ax42ty.invert_xaxis()

ax43tx = plt.axes([0,0,0.999898,1])
ip4ax = InsetPosition( ax43, [0,0,1,0.1] )
ax43tx.set_axes_locator(ip4ax)
ax43tx.axis('off')
ax43tx.hist( redshift_agn, color = 'purple', alpha = 0.2, histtype='stepfilled')
			
			
# Indicates emission line peaks
# g - band
wav_g = wav[1]
for em in emlines.keys():
	emwav = emlines[em]
	if emwav in ( 0, 912.0 ):
		continue
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
	emwav = emlines[em]
	if emwav in ( 1549.05, 1215.67, 912.0 ):
		continue
	em_z = ( wav_i / emwav ) - 1.
	legi = ax13.vlines( em_z, col_hi - 0.12 * (col_hi-col_lo), col_lo,
			colors = 'k', alpha = 1, linestyles = 'dotted' )
	ax13.text( em_z, col_hi - 0.06 * ( col_hi - col_lo ),
			em, fontsize = 10, alpha = 0.9, horizontalalignment='center', color = 'k')   

		
fig1.legend( (legp, legg, legi),
			 ('F$_{\\nu}$ ~ $\\nu^{1/3}$', '$g$ centre', '$i$ centre'),
			  bbox_to_anchor = (0.97, 0.06), loc = 'lower right' )

fig1.suptitle( 'F$_{\\nu}\sim\\nu^{%i / %i }$ Powerlaw' %(pw1,pw2) )

# Labels
ax12.text(0.4, 0.85, 'Dust Corrected (SMC)', fontsize = textsize, transform=ax12.transAxes)
ax22.text(0.4, 0.85, 'Dust Corrected (LMC)', fontsize = textsize, transform=ax22.transAxes)
ax32.text(0.4, 0.85, 'Dust Corrected (MW)', fontsize = textsize, transform=ax32.transAxes)
ax42.text(0.4, 0.85, 'Dust Corrected (AGN)', fontsize = textsize, transform=ax42.transAxes)


chi_smc = np.nansum( 5 * dat['chisq_smc'] )
chi_lmc = np.nansum( 5 * dat['chisq_lmc'] )
chi_mw = np.nansum( 5 * dat['chisq_mw'] )
chi_agn = np.nansum( 5 * dat_agn['chisq_agn'])

acc_col = astrotools.accretion_color( wav[1], wav[3], powerlaw = pwr )
sig_smc =  advtools.std( corrcol_gi_smc, acc_col, np.ones( len(corrcol_gi_smc) ) )
sig_lmc =  advtools.std( corrcol_gi_lmc, acc_col, np.ones( len(corrcol_gi_lmc) ) )
sig_mw =  advtools.std( corrcol_gi_mw, acc_col, np.ones( len(corrcol_gi_mw) ) )
sig_agn =  advtools.std( corrcol_gi_agn, acc_col, np.ones( len(corrcol_gi_agn) ) )

ax12.text(0.4, 0.75, '$\chi^{2}$ = %2.1e' %chi_smc, fontsize = textsize, transform=ax12.transAxes)
ax22.text(0.4, 0.75, '$\chi^{2}$ = %2.1e' %chi_lmc, fontsize = textsize, transform=ax22.transAxes)
ax32.text(0.4, 0.75, '$\chi^{2}$ = %2.1e' %chi_mw, fontsize = textsize, transform=ax32.transAxes)
ax42.text(0.4, 0.75, '$\chi^{2}$ = %2.1e' %chi_agn, fontsize = textsize, transform=ax42.transAxes)

ax12.text(0.4, 0.65, '$\sigma(g-i)$ = %2.4f' %sig_smc, fontsize = textsize, transform=ax12.transAxes)
ax22.text(0.4, 0.65, '$\sigma(g-i)$ = %2.4f' %sig_lmc, fontsize = textsize, transform=ax22.transAxes)
ax32.text(0.4, 0.65, '$\sigma(g-i)$ = %2.4f' %sig_mw, fontsize = textsize, transform=ax32.transAxes)
ax42.text(0.4, 0.65, '$\sigma(g-i)$ = %2.4f' %sig_agn, fontsize = textsize, transform=ax42.transAxes)

			 
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

#ax31.set_xlim( z_lo, z_hi )
ax42.set_xlim( col_lo, col_hi )
ax43.set_xlim( z_lo, z_hi )

#ax11.set_ylim( mag_lo, mag_hi )
ax12.set_ylim( mag_lo, mag_hi )
ax13.set_ylim( col_hi, col_lo )

#ax21.set_ylim( mag_lo, mag_hi )
ax22.set_ylim( mag_lo, mag_hi )
ax23.set_ylim( col_hi, col_lo )

#ax31.set_ylim( mag_lo, mag_hi )
ax32.set_ylim( mag_lo, mag_hi )
ax33.set_ylim( col_hi, col_lo )

#ax31.set_ylim( mag_lo, mag_hi )
ax42.set_ylim( mag_lo, mag_hi )
ax43.set_ylim( col_hi, col_lo )


ax12tx.set_xlim( col_lo, col_hi )
ax22tx.set_xlim( col_lo, col_hi )
ax32tx.set_xlim( col_lo, col_hi )
ax42tx.set_xlim( col_lo, col_hi )

ax12ty.set_ylim( mag_lo, mag_hi )
ax22ty.set_ylim( mag_lo, mag_hi )
ax32ty.set_ylim( mag_lo, mag_hi )
ax42ty.set_ylim( mag_lo, mag_hi )
