# Chisq plots

# Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.io import ascii

# Read in data
posf = '.txt'
files = [ ['smc','cyan'], ['lmc', 'r'], ['mw','g'], ['agn', 'purple' ]  ]
wav = np.array([3543., 4770., 6231., 7625., 9134.])

dat = ascii.read( 'Output/DustyOutput.trial.cat', format = 'fixed_width' )

# Mg vs. z, Mg vs. g-i, g-i vs. z
# Mean magnitude, disk mags, dust corr mags

labelsize = 15
textsize = 15

# Setup figure
fig1, ( ax11 ) = plt.subplots()

ax11.set_xlabel( '$\chi^{2}$', fontsize = labelsize )
ax11.set_ylabel( 'CDF', fontsize = labelsize )


ax11.minorticks_on()
ax11.tick_params( axis = 'both',
                      direction = 'in',
                      width = 2,
                      length = 6 )
ax11.tick_params( which = 'minor',
                      axis = 'both',
                      direction = 'in',)



fig1.subplots_adjust( hspace = 0, wspace = 0.2,
                      left = 0.12, right = 0.95,
                      top = 0.95, bottom = 0.12 )

colormap = plt.cm.Blues


# Bouding box values
z_lo = 0
z_hi = 4
col_lo = -0.5
col_hi = 1
mag_lo = -19
mag_hi = -27

xbins = 60
ybins = 50

xedges_z = np.linspace( z_lo, z_hi, xbins )
yedges_z = np.linspace( z_lo, z_hi, ybins )

xedges_col = np.linspace( col_lo, col_hi, xbins )
yedges_col = np.linspace( col_lo, col_hi, ybins )

xedges_mag = np.linspace( mag_hi, mag_lo, xbins )
yedges_mag = np.linspace( mag_hi, mag_lo, ybins )


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


# Mask for chosen law
for law, col in files: 
  
    # Grab relevant columns
    redshift = dat['reds']
    chisq = dat['chisq_'+law] * 5.
    chisq_sort = np.sort( chisq )
    chisq_cum = np.array( [ len( chisq_sort[0:ii] ) for ii in range( chisq.size ) ] )
    chisq_cum = chisq_cum / float( chisq.size )

    # Plot
    ax11.hist( chisq[chisq < 200], bins = 100, histtype = 'step', normed = True, color = col )
    ax11.plot( chisq_sort, chisq_cum, color = col )


# Tweak frames
ax11.set_xlim( 0, 200 )
ax11.set_ylim( 0, 1 )