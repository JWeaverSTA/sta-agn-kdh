import numpy as np 
from astropy.io import ascii
import matplotlib.pyplot as plt
import astro_funclib as astrotools
import adv_funclib as advtools

# Data Associated Arrays
# Dustlaws
dustlaws = ( 'mw', 'lmc', 'agn', 'smc')# 'lmc', 'mw', 'agn' )
dustlaws = np.array( dustlaws )
dcolor = ( 'forestgreen', 'orangered', 'blueviolet', 'royalblue' )
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


powerlaw = ( -1.2, -1.0, -0.8,
 			-0.6, -0.4, -0.2, 
			 0., 1/8.,
			 0.2, 0.25, 
			 1/3.,
			 0.4, 0.5,
			 0.6,
			 2/3., 0.75, #0.8,
			  1.0, 1.2
			  )
powerfile = ( 'n.1.2', 'n.1.0', 
			  'n.0.8',
			  'n.0.6', 'n.0.4', 'n.0.2', 
			  '0.0', '0.125',
			  '0.2', '0.25', 
			  '0.33',
			  '0.4', '0.5',
			  '0.6', 
			  '0.66', '0.75',# '0.8',
			   '1.0', '1.2'
			  )

fidwav = 2400.
pwav = np.linspace( 800, 8000, 1000 )
fig1, ax11 = plt.subplots()
ax11.set( xlabel = 'Wavelength ($\AA$)', ylabel = 'M($\lambda$)')
ax11.plot( pwav, astrotools.accretion_magspec(pwav,fidwav,powerlaw=1/3.), color = 'k', alpha = 0.9, lw = 2, ls='dashed' )

ax11.invert_yaxis()
fig2, (ax21,ax22) = plt.subplots(nrows = 1, ncols = 2)
ax22.set( xlabel = 'Powerlaw index', ylabel = 'Log $\sigma$(g-i)')
ax21.set( xlabel = 'Powerlaw index', ylabel = 'Log $\chi^{2}$' )

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

	name = 'Output/DustyOutput.powerlaw.fakedither.2sig.'+powerfile[i]+'.cat'

	dat = ascii.read( name, format = 'fixed_width' )

	for j, k in enumerate( dustlaws ):
		totalchisq = np.nansum( 5 * dat['chisq_'+k] )
		holdchi[i, j] = totalchisq

		gi_col = dat['g_dered_'+k+'_amag'] - dat['i_dered_'+k+'_amag']
		std_gi = advtools.std( gi_col, acc_col, np.ones( len(gi_col) ) )
		holdstd[i,j] = std_gi

ylo1, yhi1 = 0.7*np.min(holdchi), 1.5*np.max(holdchi)
ylo2, yhi2 =  0.7*np.min(holdstd), 1.5*np.max(holdstd)
ax21.vlines( 1/3., ylo1, yhi1, 'k', linestyles = 'dotted' )
ax22.vlines( 1/3., ylo2, yhi2, 'k', linestyles = 'dotted' )

for i, j in enumerate( holdchi.T ):
	ax21.scatter( holdpwr, holdchi.T[i], c = dcolor[i], marker = 'o', label = dustlaws[i].upper() )
	ax22.scatter( holdpwr, holdstd.T[i], c = dcolor[i], marker = 'o')

	minchi = np.min( holdchi.T[i] )
	maxchi = np.max( holdchi.T[i] )
	for k in range(len(holdpy)):
		ax11.plot( pwav, holdpy[k], color = 'k', alpha = 0.1, lw = 1 )
		if k in (0, 2, 4, 6, 10, 16, 18 ):
			if k in (0, 2, 4):
				ax11.text( 8400, holdpy[k][-1], '-%2.1f' %powerlaw[k], fontsize = 9 )
			else:
				ax11.text( 8500, holdpy[k][-1], '%2.1f' %powerlaw[k], fontsize = 9 )


	a = np.polyfit( holdpwr, holdchi.T[i], 5 )
	b = np.poly1d(a)
	
	astd = np.polyfit( holdpwr, holdstd.T[i], 5 )
	bstd = np.poly1d(astd)
	
	pyx = np.linspace( holdpwr[0], holdpwr[-1], 10000 )
	
	pyb = b(pyx)
	pybstd = bstd(pyx)
	
	ax21.plot( pyx, pyb, color = dcolor[i] )
	ax22.plot( pyx, pybstd, color = dcolor[i] )
	
	minx = pyx[pyb==min(pyb)][0]

	minxstd = pyx[pybstd == min(pybstd)][0]

	for o in range(np.where(pyb==min(pyb))[0], len(pyb) ):
		x = pyx[o]
		y = pyb[o]
		dof = (5-2) * len(dat)
		if y - min(pyb) > ( min(pyb) / dof ):
			break
	sigx = x - minx

	for o in range(np.where(pybstd==min(pybstd))[0], len(pybstd) ):
		x = pyx[o]
		y = pybstd[o]
		dof = (5-2) * len(dat)
		if y - min(pybstd) > ( min(pybstd) / dof ):
			break
	sigxstd = x - minxstd

	ax21.fill_between( [ minx - sigx, minx + sigx ], [ylo1, ylo1], [yhi1, yhi1],
						color = dcolor[i], alpha = 0.2 )
	ax21.vlines(minx, ylo1, yhi1, linestyles = 'dotted', colors = dcolor[i])

	ax22.fill_between( [ minxstd - sigxstd, minxstd + sigxstd ], [ylo2, ylo2], [yhi2, yhi2],
						color = dcolor[i], alpha = 0.2 )
	ax22.vlines(minxstd, ylo2, yhi2, linestyles = 'dotted', colors = dcolor[i])


	print 'for %s : accexp = %2.5f +/- %2.5f | chisq = %2.3e | dof = %2.5f | chisq/dof = %2.5f '%(dustlaws[i], minx, sigx, min(pyb), dof, min(pyb)/dof )


	pam = astrotools.accretion_magspec( pwav, fidwav, powerlaw = minx )
	pamhi = astrotools.accretion_magspec( pwav, fidwav, powerlaw = minx + sigx )
	pamlo = astrotools.accretion_magspec( pwav, fidwav, powerlaw = minx - sigx )
	ax11.plot(pwav, pam, color=dcolor[i], lw = 3 )
	ax11.fill_between( pwav, pamlo, pamhi, color = dcolor[i], alpha = 0.5)




ax21.set_ylim( ylo1, yhi1 )
ax22.set_ylim( ylo2, yhi2 )
ax21.set_yscale('log')
ax22.set_yscale('log')

ax21.legend()
