# PRINTS DUST LAWS
# USED TO CHECK IF CORRECT

import numpy as np
import matplotlib.pyplot as plt
import dust_funclib as dusttools
from astropy.io import ascii


ls_wav = np.linspace( 2E5, 0.1E4, 1000000 )
ls_x = 1E4 / ( ls_wav ) #inverse microns
vwav = 5550.
bwav = 4361.

vwavz2 = vwav / ( 1 + 2. )

# Alam/Av
ebmv = 1

###
#rv_smc = 5
#rv_lmc = 1
###

smcspl = dusttools.smccubic()
lmcspl = dusttools.lmccubic()

vwav = np.array([vwav,0.1])

smc_vwav = dusttools.dustlaw( vwav, 'smc', smcspl )[0]
lmc_vwav = dusttools.dustlaw( vwav, 'lmc', lmcspl )[0]
mw_vwav = dusttools.dustlaw( vwav, 'mw' )[0]
gask_vwav = dusttools.dustlaw( vwav, 'agn' )[0]

ls_smc  = dusttools.dustlaw( ls_wav, 'smc', smcspl )
ls_lmc  = dusttools.dustlaw( ls_wav, 'lmc', lmcspl )
#ls_lmcss= np.array( [ func_lmcss( i, lmcss_ebmv ) for i in ls_wav ] )
ls_mw   = dusttools.dustlaw( ls_wav, 'mw' )
ls_gask = dusttools.dustlaw( ls_wav, 'agn' )

p_smc = ( ls_smc )  #/ ( smc_vwav + 1 ) # / ( ( smc_vwav / rv_smc ) + 1 )
p_lmc = ( ls_lmc  ) #/ ( lmc_vwav + 1 ) # / ( ( lmc_vwav / rv_lmc ) + 1 )
#p_lmcss = ( ls_lmcss )

p_mw = ls_mw / mw_vwav
p_gask = ls_gask / gask_vwav

#p_smc = ( ls_smc / rv_smc ) #/ ( ( smc_vwav + 1 )  )
#p_lmc = ( ls_lmc + 1 ) #/ ( ( lmc_vwav + 1 )  )


# import Gordon data
gor_smc = ascii.read( 'smcbar_ext.dat' )
gsx = gor_smc['col1']
gsy = gor_smc['col2']
gse = gor_smc['col3']

gor_lmc = ascii.read( 'lmcave_ext.dat' )
glx = gor_lmc['col1']
gly = gor_lmc['col2']
gle = gor_lmc['col3']



# Plotting
fig1, ( ax1 ) = plt.subplots( )

px = ls_x
py = [ p_smc, p_lmc, p_mw, p_gask  ]
plabel = [ 'SMC', 'LMC', 'MW', 'AGN' ]
pcol = [ 'cyan', 'r', 'g', 'purple' ]
#py = [ psmc, plmc, pmw, pgask ]
#py = [ ls_smc, ls_lmc, ls_mw, ls_gask ]
#plabel = [ 'SMC', 'LMC', 'MW', 'GASKELL' ]
#pcol = [ 'cyan', 'r', 'g', 'b' ] 

for i in range( 4 ):
  ax1.plot( px, py[i], c = pcol[i], label = plabel[i] )  
  
plt.legend()
  
plt.errorbar( gsx, gsy, yerr=gse, c='blue', fmt = 'o', ms = 4 )#,# s=2, zorder=5, label = 'LMC-GOR' ) #alpha = 0.5, linestyle = 'dashed' )
plt.errorbar( glx, gly, yerr=gle, c='maroon', fmt = 'o', ms = 4 )#,# s = 2, zorder=5, label = 'LMC-GOR'  ) #, alpha = 0.5, linestyle = 'dashed' )
#plt.plot( gsx, gsy, c='blue', alpha = 0.8, lw=0.5,linestyle = 'solid' )
#plt.plot( glx, gly, c='maroon', alpha = 0.8, lw = 0.5, linestyle = 'solid' )    
  
  
plt.plot( [1/.55,10], [1,10], ls = 'dashed' , c = 'gray' )  
plt.vlines( 1E4 / vwav, 0, 9, colors = 'g', linestyles = 'dashed' )
plt.vlines( 1E4 / vwavz2, 0, 9, colors = 'orange', linestyles = 'dashed' )
plt.vlines( 1E4 / bwav, 0, 9, colors = 'b', linestyles = 'dashed' )
plt.vlines( 3.29, 0, 9, colors = 'k', lw = 0.5, linestyles = 'dotted' )
#plt.vlines( 1E4 / bwav, 0, 9, colors = 'b', linestyles = 'dotted' )
#plt.hlines( 0, px[0], px[-1], colors = 'k', linestyles = 'dotted' )
plt.hlines( 1, px[0], px[-1], colors = 'k', linestyles = 'dotted' )
plt.xlim( px[0], gsx[-1] )

#ls_w = np.array( [3543.0, 4770.0, 6231.0, 7625.0, 9134.0] )
#plt.vlines( 1E4 / ls_w, 0, 9, colors = 'k', linestyles = 'solid' )

plt.xlabel( '$\lambda^{-1}$ ($\mu$m$^{-1}$)', fontsize = 15 )
plt.ylabel( 'A$_{\lambda}$ / A$_{V}$', fontsize = 15 )

#plt.title( 'Spline | SMC + LMC' )


plt.ylim( 0, 9 )

plt.draw()
