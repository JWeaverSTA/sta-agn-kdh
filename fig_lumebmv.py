# Figure - L_fid v. E(B-V)

# Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from astropy.io import ascii

# Read in data
posf = '.txt'
files = [ ['smc','cyan'], ['lmc', 'r'], ['mw','g'], ['gask', 'purple' ]  ]
wav = np.array([3543., 4770., 6231., 7625., 9134.])

dat = ascii.read( 'Output/DustyOutput.trial.cat', format = 'fixed_width' )

# Mg vs. z, Mg vs. g-i, g-i vs. z
# Mean magnitude, disk mags, dust corr mags

labelsize = 15
textsize = 15

# Setup figure
fig1, ( ax11 ) = plt.subplots()

ax11.set_xlabel( 'E(B-V)', fontsize = labelsize )
ax11.set_ylabel( 'L$_{2400\AA}$ erg s$^{-1}$', fontsize = labelsize )


ax11.minorticks_on()
ax11.tick_params( axis = 'both',
                      direction = 'in',
                      width = 2,
                      length = 6 )
ax11.tick_params( which = 'minor',
                      axis = 'both',
                      direction = 'in',)



fig1.subplots_adjust( hspace = 0, wspace = 0.2,
                      left = 0.15, right = 0.98,
                      top = 0.95, bottom = 0.1 )

colormap = plt.cm.gnuplot2
colormap.set_under('w')


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



# Grab relevant columns
redshift = dat['reds']

corrmag_g = dat['g_dered_smc_amag']
corrmag_i = dat['i_dered_smc_amag']
corrcol_gi = corrmag_g - corrmag_i

ebmv = dat['ebmv_smc']

mstar = dat['fidmag_smc']

# Convert mstar -> Lfid
Lsun = 3.9E33
Lfid = Lsun * 10**( ( 4.77 - mstar ) / 2.5 )


# Plot
H, xedges, yedges = np.histogram2d(
             ebmv, np.log10(Lfid), bins = (50,50),
             range = ( ( -0.2, 1 ), (43.5, 46.5 ) ) )
ax11.imshow( H.T, extent = (xedges[0], xedges[-1], yedges[0], yedges[-1]),
             interpolation = 'gaussian', vmin = 2, cmap = colormap, aspect = 1/5. )
ax11.vlines( 0, 43.5, 46.5,
             colors='r',linestyles='--' )

fig1.subplots_adjust( left = 0.1 )
# Spline optavg
'''
px = yedges[0:-1:2] + 1.
py = np.log10( np.array(Lfid) ) 
inds = np.argsort( py )
pys = py[inds]
pxs = ebmv[inds]
px_avg = [np.mean(pxs[(pys>yedges[i])&(pys<yedges[i+2]) ]) for i in xrange(0,len(yedges)-1,2) ]

ax11.plot( px_avg , px )
             
'''


