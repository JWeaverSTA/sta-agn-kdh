# Figure - L_fid v. E(B-V)

# Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.io import ascii
import adv_funclib as advtools


# Read in data
posf = '.txt'
files = [ ['agn', 'blueviolet' ],
      ['mw','forestgreen'],
      ['lmc', 'orangered'],
      ['smc','royalblue']
      ]
wav = np.array([3543., 4770., 6231., 7625., 9134.])

dat = ascii.read( 'Output/DustyOutput.trial.cat', format = 'fixed_width' )

# Mg vs. z, Mg vs. g-i, g-i vs. z
# Mean magnitude, disk mags, dust corr mags

labelsize = 15
textsize = 15

# Setup figure
fig1, ( ax11 ) = plt.subplots()

ax11.set_xlabel( 'E(B-V)', fontsize = labelsize )
ax11.set_ylabel( 'CDF', fontsize = labelsize )


ax11.minorticks_on()

fig1.subplots_adjust( hspace = 0, wspace = 0.2,
                      left = 0.08, right = 0.98,
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

cbins = np.linspace( -0.5, 2, 100)
# Mask for chosen law
for law, col in files: 
    # Grab relevant columns
    redshift = dat['reds']
    ebmv = dat['ebmv_'+law]
    ebmv_sort = np.sort( ebmv )
    ebmv_cum = np.array( [ len( ebmv_sort[0:ii] ) for ii in range( ebmv.size ) ] ) / float( ebmv.size )
    #ebmv_cum = ebmv_cum / float( ebmv.size )

    lawname = law.upper()

    if law == files[0][0]:
      old_ebmv_cum = np.zeros( len( ebmv_cum ))
      old_ebmv_sort = ebmv_sort 

    old_ebmv_cum = [advtools.linint( old_ebmv_sort, old_ebmv_cum,
                                    np.ones(len(old_ebmv_sort)), o)[0] for o in ebmv_sort ]

    # Plot
    #ax11.hist( ebmv, bins = bins, histtype = 'step', normed = True, color = col, label = None )
    ax11.fill_between( ebmv_sort, old_ebmv_cum, alpha = 0.2, color = col )
    ax11.plot( ebmv_sort, 3.5 *  ebmv_cum, color = col, label = lawname )

    hist, bins_ = np.histogram( ebmv, bins = cbins )
    ax11.step( bins_[:-1], hist/1300., color = col )

    old_ebmv_sort = ebmv_sort
    old_ebmv_cum = ebmv_cum


# Vline
ax11.vlines( 0, 0, 3.5,
             colors='k',linestyles='--' )
             
ax11.legend()             

# Tweak frames
ax11.set_xlim( -0.5, 1.5 )
ax11.set_ylim( 0, 3.5 )