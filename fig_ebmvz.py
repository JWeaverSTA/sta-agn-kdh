from astropy.io import ascii
import numpy as np 
import matplotlib.pyplot as plt 

# NEED TO REMOVE THE NON UNIFORM NUMBER DISTRIBUTION FROM EQUATION!

# Get data
dat = ascii.read('Output/DustyOutput.trial.cat', format = 'fixed_width')

# Get columns
ebmv = dat['ebmv_smc']
z = dat['reds']

# Sort
inds = z.argsort()
ebmv = ebmv[inds]
z = z[inds]

ebin1 = 0.5

nodust = z[ebmv<0]
lodust = z[(ebmv<ebin1) & (ebmv>=0)]
hidust = z[ebmv>=ebin1]

zbins = np.linspace( 0, max(z), 10 )
nbins = np.zeros( len(zbins) - 1 )

for i in ( nodust, lodust, hidust ):

	hist, bin_edges = np.histogram( i, zbins )

	nbins = nbins + hist


fig, ax = plt.subplots( nrows = 2, ncols = 1, sharex =True,
						gridspec_kw = {'height_ratios':[2,1] } )
mycmap = plt.cm.magma_r
mycmap.set_under('w')

ax[0].hist2d( z, ebmv, bins = (np.linspace(0,max(z), 50), np.linspace( -1, 4, 30) ),
              cmin = 1, cmap = mycmap )
ax[0].set_ylabel('E(B-V)')

ax[1].set_xlabel('Redshift z')
ax[1].set_ylabel('Fraction of AGN')

names = ('E(B-V)<0', '0$\leq$E(B-V)<%2.1f'%ebin1, 'E(B-V)$\geq$%2.1f'%ebin1 )
color = ('b', 'orange', 'red')
for j,i in enumerate( ( nodust, lodust, hidust ) ):

	hist, bin_edges = np.histogram( i, zbins )
	print bin_edges


	bin_edges = (bin_edges - (bin_edges[1] - bin_edges[0])/2. )[1:]

	phist = hist / nbins

	ax[1].plot( bin_edges, phist, marker = 's', label = names[j], c = color[j] )

ax[1].plot( bin_edges, nbins/nbins.max(), 'k-', alpha = 0.5, label = 'Fraction in bin', zorder = 0)
ax[1].set_ylim(0,1)
ax[0].set_ylim(-0.5, 3)
ax[1].legend(loc=7)
ax[1].set_xlim(0,1.02*max(z))
fig.subplots_adjust( hspace = 0)

