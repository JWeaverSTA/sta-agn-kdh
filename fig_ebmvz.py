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


fig, ax = plt.subplots()
ax.set_xlabel('Redshift z')
ax.set_ylabel('Fraction of AGN')

names = ('E(B-V)<0', '0$\leq$E(B-V)<%2.1f'%ebin1, 'E(B-V)$\geq$%2.1f'%ebin1 )
color = ('b', 'orange', 'red')
for j,i in enumerate( ( nodust, lodust, hidust ) ):

	hist, bin_edges = np.histogram( i, zbins )
	print bin_edges


	bin_edges = (bin_edges - (bin_edges[1] - bin_edges[0])/2. )[1:]

	phist = hist / nbins

	ax.plot( bin_edges, phist, marker = 's', label = names[j], c = color[j] )

ax.legend(loc=7)
ax.set_xlim(0,1.02*max(z))

