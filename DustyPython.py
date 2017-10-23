# DustyPyton
# For decomposing AGN lightcurves
# Created 23.10.17

# Imports
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table, Column

# Input Parameters
lightcurvedir = 'QSO_S82/'
masterdir = 'QSO_Master/'
masterfile = 'DB_QSO_S82.dat'
outputdir = 'DustyPyton/'
outputfile = 'DusyOutput.junk.cat'

# Data Associated Arrays
# Wavelengths
wav = { 'u' : 3543.,
        'g' : 4770.,
        'r' : 6231.,
        'i' : 7625.,
        'z' : 9134.}
# Wavlength colors
wavcolor = { 'u' : 'b',
             'g' : 'g',
             'r' : 'orange',
             'i' : 'r',
             'z' : 'k'}
# Emission lines
chem = { 'Ly$_{\infty}$'  : 912.00,
         'Ly$_{\\alpha}$' : 1215.67,
         'CIV'            : 1549.05,
         '[CIII]'         : 1908.73,
         '[MgII]'         : 2797.92,
         'H$_{\infty}$'   : 3560.00,
         'H$_{\\beta}$'   : 4861.32,
         'H$_{\\alpha}$'  : 6562.80}

# Load Masterfile
masterdata = ascii.read( masterdir + masterfile )
nagn = np.arange( 1, len( masterdata ) + 1 )

# Main Menu Navigation
# Inital Sample
iagn1 = 1
iagn2 = nagn
zlo = 0
zhi = 6
gfnt = 16
gbrt = 23
cblu = -1.
cred = 6.
ralo = 0.
rahi = 24.
declo = -2.
dechi = 2.
numlook = 0
output = 'OFF'

# Colors
gi_color = masterdata['g'] - masterdata['i']

# Inital Values
ind_agn = ( nagn >= iagn1 )             & ( nagn < iagn2 )
ind_z   = ( masterdata['reds'] > zlo )  & ( masterdata['reds'] < zhi )
ind_g   = ( masterdata['g'] > gfnt )    & ( masterdata['g'] < gbrt )
ind_c   = ( gi_color > cblu )           & ( gi_color < cred )
ind_ra  = ( masterdata['ra'] > ralo )   & ( masterdata['ra'] < rahi )
ind_dec = ( masterdata['dec'] > declo ) & ( masterdata['dec'] < dechi )
ind_nl  = ( masterdata['dbID'] == numlook )

ind_init = np.logical_and.reduce( [ ind_agn, ind_z, ind_g,
                                    ind_c, ind_ra, ind_dec ] )

# Interface

