# Dust Model Toolkit
# Created 30.10.17

# Imports
import numpy as np
from astropy.io import ascii
from scipy import interpolate
from numba import njit

# --- Gordon Interpolation ---
def smccubic():

	gor_smc = ascii.read( 'Dustdata/smcbar_ext.dat' )
	gsx = gor_smc['col1'][:10]
	gsy = gor_smc['col2'][:10]
	gse = gor_smc['col3'][:10]
	
	gsx[9] = 3.29
	gsy[9] = smc( 1E4/ 3.29, None )
	gse[9] = 99.
	
	k = 3
	

	# for two-point spline...
	gsx = [ 0.43, 1E4/5550., 3.29 ]
	gsy = [ 0., 1., smc( 1E4/ 3.29, None ) ]
	gse = [ 1., 1., 1. ]
	k = 2
	

	tck = interpolate.splrep( gsx, gsy, w = gse, s = 0, k = k )
		
	return tck

def lmccubic():

	gor_lmc = ascii.read( 'Dustdata/lmcave_ext.dat' )
	glx = gor_lmc['col1'][:8]
	gly = gor_lmc['col2'][:8]
	gle = gor_lmc['col3'][:8]
	
	glx[7] = ( 3.29 )
	gly[7] = lmc( 1E4 / 3.29, None )
	gle[7] = 99.
	
	k = 3
	
	
	# for two-point spline...
	glx = [ 0.43, 1E4/5550., 3.29 ]
	gly = [ 0., 1., lmc( 1E4/ 3.29, None ) ]
	gle = [ 1., 1., 1. ]
	k = 2
	

	tck = interpolate.splrep( glx, gly, w = gle, s = 0, k = k )
		
	return tck


# --- D FACTOR ---
@njit
def d( x, x0, gam ):
  a = ( x*x ) - ( x0*x0 )
  b = x*gam
  d = ( x*x ) / ( ( a*a ) + ( b*b ) )

  return d

# --- Fx FACTOR ---
@njit
def fx( x ):
  fa = x - 5.9
  faa = fa * fa
  fx = 0.5392*faa + 0.05644*faa*fa
  if x < 5.9:
    fx = 0.

  return fx

# --- SMC LAW ---
def smc( wav, tck ):

  rv_smc = 2.74

  c1 = -4.959 
  c2 = 2.264  
  c3 = 0.389
  c4 = 0.461  
  x0 = 4.6
  gam = 1.0

  wav_micron = wav/1E4
  x = 1./wav_micron
  
  if x >= 3.29:
  
  	df = d( x, x0, gam )
  	fxf = fx( x )
  
  	eratio = c1 + c2*x + c3*df + c4*fxf

  	return ( eratio / rv_smc + 1 )

  if x < 3.29:
	
  	ynew = interpolate.splev( x, tck, 0 )
		
  	return ynew


  '''
  xvm = 1 / 0.54
  xbm = 1 / 0.442
  xwi = 0.0

  a0 = c1 + c2 * xwi + c3 * d( xwi, x0, gam ) + c4 * fx( xwi )
  avm = c1 + c2 * xvm + c3 * d( xvm, x0, gam ) + c4 * fx( xvm )
  abm = c1 + c2 * xbm + c3 * d( xbm, x0, gam ) + c4 * fx( xbm )
  ebmv_calc = abm - avm

  extmag_smc = ( eratio - a0 ) / ebmv_calc * ebmv  
  '''



# --- LMC LAW ---
def lmc( wav, tck ):

  rv_lmc = 3.41
  
  c1 = -0.890
  c2 = 0.998
  c3 = 2.719
  c4 = 0.400
  x0 = 4.579
  gam = 0.934
  
  wav_micron = wav / 1.E4
  x = 1./ wav_micron
  
  if x >= 3.29:
    
  	df = d( x, x0, gam )
  	fxf = fx( x )

  	eratio = c1 + c2*x + c3*df + c4*fxf
  	
  	return ( eratio /rv_lmc + 1 )
  	
  if x < 3.29:
  
  	ynew = interpolate.splev( x, tck, 0 )
  	
  	return ynew

  '''
  xvm = 1 / 0.54
  xbm = 1 / 0.442
  xwi = 0.0

  a0 = c1 + c2 * xwi + c3 * d( xwi, x0, gam ) + c4 * fx( xwi )
  avm = c1 + c2 * xvm + c3 * d( xvm, x0, gam ) + c4 * fx( xvm )
  abm = c1 + c2 * xbm + c3 * d( xbm, x0, gam ) + c4 * fx( xbm )
  ebmv_calc = abm - avm

  extmag_lmc = ( eratio - a0 ) / ebmv_calc * ebmv
  '''
  

  
# --- LMC Supershell LAW ---
# OBSOLETE! DO NOT USE!!!
def lmcss( wav ):

  rv_lmcss = 2.76
  
  c1 = -1.475
  c2 = 1.132
  c3 = 1.463
  c4 = 0.294
  x0 = 4.558
  gam = 0.945
  
  wav_micron = wav / 1.E4
  x = 1./ wav_micron
    
  d = d( x, x0, gam )
  fx = fx( x )
  
  eratio = c1 + c2*x + c3*d + c4*fx
		
  return ( eratio / rv_lmcss + 1 )

  '''
  xvm = 1 / 0.54
  xbm = 1 / 0.442
  xwi = 0.0

  a0 = c1 + c2 * xwi + c3 * d( xwi, x0, gam ) + c4 * fx( xwi )
  avm = c1 + c2 * xvm + c3 * d( xvm, x0, gam ) + c4 * fx( xvm )
  abm = c1 + c2 * xbm + c3 * d( xbm, x0, gam ) + c4 * fx( xbm )
  ebmv_calc = abm - avm

  extmag_lmc = ( eratio - a0 ) / ebmv_calc * ebmv
  '''




# --- MW LAW ---
@njit
def mw( wav ):

  ntable = 19
  xtable = np.array( [ 0., 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7 ] )
  etable = np.array( [ 0., 1.36, 1.64, 1.84, 2.04, 2.24, 2.44, 2.66, 2.88, 3.14, 3.36, 3.56, 3.77, 3.96, 4.15, 4.26, 4.40, 4.52, 4.64 ] )

 
  x = 1E4 / wav

  # infrared
  if x < 1.0 :
    extmag = etable[1] * x * x

  # optical
  elif x < 2.7 :
    indhi = np.where( xtable > x )[0][0]
    indlo = np.where( xtable < x )[0][-1]
    part = xtable[indhi] - xtable[indlo]
    if part != 0:
      part = ( x - xtable[indlo] ) / part
    extmag = etable[indlo] * ( 1 - part ) + etable[indhi] * part
    

  # ultraviolet
  elif x < 3.65 :
    diff = x - 4.6
    extmag = 1.56 + 1.048 * x + 1.01 / ( diff * diff + 0.280 )

  elif x < 7.14 :
    diff = x - 4.6
    extmag = 2.29 + 0.848 * x + 1.01 / ( diff * diff + 0.280 )

  elif x < 10 :
    extmag = 16.17 + x * ( -3.20 + 0.2975 * x )

  # fuv
  else:
    x = min( x, 50. )
    extmag = 16.17 + x * ( -3.2 + 0.2975 * x )

  return extmag / 3.143963963963964 #normalise to vband


# --- GASKELL LAW ---
@njit
def agn(wav):
  
  x = 1E4 / wav

  if ( x < 1.6 ):
    x0 = 1.6
    y0 = x0 * ( x0 * ( x0 * 0.0296 - 0.377 ) + 1.5848 ) - 0.8175
    dydx = x0 * ( x0 * 3*0.0296 - 2*0.377 ) + 1.5848
    y = y0 + dydx * ( x - x0 )
  elif ( x < 3.69 ):
    y = x * ( x * ( x * 0.0296 - 0.377 ) + 1.5848 ) - 0.8175
  elif ( x < 8. ):
    y = 1.3468 + 0.0087 * x
  else:
    y = 1.3468 + 0.0087 * x

  rv = 5.15
  av = rv * y

  return av / 5.084157174065957 #normalise to vband


# --- DRUDE PROFILE ---
def drude( wav, c1, c2, c3, c4, ebmv ):
  t1 = c1 / ( ( wav / 0.08 )**c2 + ( 0.08 / wav )**c2 + c3 )
  t2t = 233. * ( 1 - c1 / ( 6.88**c2 + 0.145**c2 + c3 ) - c4 / 4.60 )
  t2b = ( wav / 0.046 )**2 + ( 0.046 / wav )**2 + 90.
  t2 = t2t / t2b
  t3 = c4 / ( ( wav / 0.2175 )**2 + ( 0.2175 / wav )**2 - 1.95 )
  
  return t1 + t2 + t3



def dustlaw( wav, name = None, tk = None ):
	if name == 'smc':
		return np.array( [ smc( i, tk ) for i in wav ] )
	elif name == 'lmc':
		return np.array( [ lmc( i, tk ) for i in wav ] )
	elif name == 'mw':
		return np.array( [ mw( i ) for i in wav ] )
	elif name == 'agn':
		return np.array( [ agn( i ) for i in wav ] )
	elif name == 'None':
		return None