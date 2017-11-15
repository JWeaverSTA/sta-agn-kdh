import numpy as np 
from astropy.io import ascii
import matplotlib.pyplot as plt
import astro_funclib as astrotools
import adv_funclib as advtools

# Data Associated Arrays
# Dustlaws
dustlaws = ( 'smc', 'lmc', 'mw', 'agn' )
dustlaws = np.array( dustlaws )
dcolor = ( 'royalblue', 'orangered', 'forestgreen', 'blueviolet' )
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


powerlaw = ( -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 
	          0., 1/8., 0.2, 0.25, 1/3., 0.4, 0.5,
	           0.6, 2/3., 0.75, 0.8, 1.0, 1.2 )
powerfile = ( 'nn.1.2', 'nn.1.0', 'nn.0.8', 'nn.0.6', 'nn.0.4', 'nn.0.2', 
	          'n.0.0', '0.125', 'n.0.2', '0.25', '0.33', 'n.0.4', '0.5',
	          'n.0.6', '0.66', '0.75', 'n.0.8', 'n.1.0', 'n.1.2')

fidwav = 2400.
pwav = np.linspace( 800, 8000, 1000 )
fig1, ax11 = plt.subplots()
ax11.invert_yaxis()
fig2, (ax21,ax22) = plt.subplots(nrows = 1, ncols = 2)
ax22.set( xlabel = 'Accretion Exponent', ylabel = 'Log $\sigma$(g-i)')
ax21.set( xlabel = 'Accretion Exponent', ylabel = 'Log $\chi^{2}$' )

holdchi = np.zeros( shape = ( len( powerlaw ), 4 ) )
holdstd = np.zeros( shape = holdchi.shape )
holdpwr = np.zeros( shape = len(powerlaw) )

holdpy = np.zeros( shape = ( len(powerlaw), len(pwav) ) )

for i, p in enumerate( powerlaw ):

	pwr = p
 	py = astrotools.accretion_magspec( pwav, fidwav, powerlaw = pwr )
	holdpwr[i] = pwr
	holdpy[i] = py

	acc_col = astrotools.accretion_color( wav[1], wav[3], powerlaw = pwr )

	name = 'Output/DustyOutput.powerlaw.'+powerfile[i]+'.cat'

	dat = ascii.read( name, format = 'fixed_width' )

	for j, k in enumerate( dustlaws ):
		totalchisq = np.nansum( 5 * dat['chisq_'+k] )
		holdchi[i, j] = totalchisq

		gi_col = dat['g_dered_'+k+'_amag'] - dat['i_dered_'+k+'_amag']
		std_gi = advtools.std( gi_col, acc_col, np.ones( len(gi_col) ) )
		holdstd[i,j] = std_gi

for i, j in enumerate( holdchi.T ):
	ax21.plot( holdpwr, holdchi.T[i], c = dcolor[i], marker = 'o', label = dustlaws[i] )
	ax22.plot( holdpwr, holdstd.T[i], c = dcolor[i], marker = 'o')
	if i == 0:
		minchi = np.min( holdchi.T[i] )
		maxchi = np.max( holdchi.T[i] )
		for k in range(len(holdpy)):
			ax11.plot( pwav, holdpy[k], color = dcolor[i] )

		a = np.polyfit( holdpwr, holdchi.T[i], 2 )
		b = np.poly1d(a)
		pyx = np.linspace( holdpwr[0], holdpwr[-1], 10000 )
		pyb = b(pyx)
		#ax21.plot(pyx, pyb, 'r--')

ax21.set_yscale('log')
ax22.set_yscale('log')

ax21.legend()
