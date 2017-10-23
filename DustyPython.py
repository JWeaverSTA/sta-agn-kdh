# DustyPyton
# For decomposing AGN lightcurves
# Created 23.10.17

# Imports
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table, Column
import astro_funclib as astrotools
import adv_funclib as advtools

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
masterdata['ra'] = astrotools.deg2hr( masterdata['ra'] )
nagn = np.arange( 1, len( masterdata ) + 1 )
nagn_init = len( nagn ) + 1

# Main Menu Navigation
# Inital Sample
iiagn1 = 1
iiagn2 = len( nagn ) + 1
izlo = 0
izhi = 5
igbrt = 16
igfnt = 23
icblu = -1.
icred = 2.
iralo = 0
irahi = 24.
ideclo = -1.27
idechi = 1.27
inumlook = 0
ioutput = 'OFF'
iverbose = 'OFF'

# Colors
gi_color = masterdata['g'] - masterdata['i']

# Transfer inital values for use in loop
iagn1 = iiagn1
iagn2 = iiagn2
zlo = izlo
zhi = izhi
gfnt = igfnt
gbrt = igbrt
cblu = icblu
cred = icred
ralo = iralo
rahi = irahi
declo = ideclo
dechi = idechi
numlook = inumlook
output = ioutput
verbose = iverbose

# Interface
truth_menu = 0
truth_nav = 0
while truth_nav not in ('end','exit','q'):

    while truth_menu not in ('P','L'):
    
        # Define Indices
        ind_agn = ( nagn >= iagn1 )             & ( nagn < iagn2 )
        ind_z   = ( masterdata['reds'] > zlo )  & ( masterdata['reds'] < zhi )
        ind_g   = ( masterdata['g'] > gbrt )    & ( masterdata['g'] < gfnt )
        ind_c   = ( gi_color > cblu )           & ( gi_color < cred )
        ind_ra  = ( masterdata['ra'] > ralo )   & ( masterdata['ra'] < rahi )
        ind_dec = ( masterdata['dec'] > declo ) & ( masterdata['dec'] < dechi )
        ind_nl  = ( masterdata['dbID'] == numlook )
    
        index = np.logical_and.reduce( [ ind_agn, ind_z, ind_g,
                                    ind_c, ind_ra, ind_dec ] )
                             
        # Override index
        if numlook != 0:
            index = ind_nl
        if iagn1 == iagn2:
            index = iagn1 - 1
        
        # Redefine table    
        selectdata = masterdata[index]
        
        if type( index ) == int:
            nselect = 1
        else:
            nselect = len( selectdata )

        # Report to user
        print ' '
        print 'SAMPLE SIZE: ', nselect
        print ' '
        print 'SELECTION CRITERIA:'
        print '  AGN:', iagn1, iagn2, 'OF', nagn_init
        print '  Z  :', zlo, zhi, ' REDSHIFT'
        print '  MAG:', gbrt, gfnt, ' SLOAN g'
        print '  COL:', cblu, cred, ' SLOAN g-i'
        print '  RA :', ralo, rahi, ' HR'
        print '  DEC:', declo, dechi, ' DEG'
        print 'SDSSN:', numlook
        print ' '
        print 'VERBOSE...'+verbose
        print 'OUTPUT....'+output
        if output == 'ON':
            print 'FILE NAME: ', outtab_name
        print ' '
        if len( selectdata )>0:
            print 'P.....Print'
            print 'L.....List'
        if len( selectdata )==0:
            print 'ZERO SIZE SAMPLE - TRY AGAIN.'
        print ' '
      
        # Ask for changes
        truth_menu = raw_input('CHANGE: ')
        if truth_menu in ('end','exit','q'):
            sys.exit(0)

        if truth_menu == 'AGN':
            iagn1 = int(raw_input('  lower: '))
            iagn2 = int(raw_input('  upper: '))
        if truth_menu == 'Z':
            zlo = float(raw_input('  lower: '))
            zhi = float(raw_input('  upper: '))
        if truth_menu == 'MAG':
            gbrt = float(raw_input('  lower: '))
            gfnt = float(raw_input('  upper: '))
        if truth_menu == 'COL':
            cblu = float(raw_input('  lower: '))
            cred = float(raw_input('  upper: '))
        if truth_menu == 'RA':
            ralo = float(raw_input('  lower: '))
            rahi = float(raw_input('  upper: '))
        if truth_menu == 'Dec':
            declo = float(raw_input('  lower: '))
            dechi = float(raw_input('  upper: '))
        if truth_menu == 'SDSSN':
            numlook = int(raw_input('  number: '))
        if truth_menu == 'VERBOSE':
            verbose = str(raw_input('  ON/OFF: '))       
        if truth_menu == 'OUTPUT':
            output = str(raw_input('  ON/OFF: '))
        if truth_menu == 'reset':
            iagn1 = iiagn1
            iagn2 = iiagn2
            zlo = izlo
            zhi = izhi
            gfnt = igfnt
            gbrt = igbrt
            cblu = icblu
            cred = icred
            ralo = iralo
            rahi = irahi
            declo = ideclo
            dechi = idechi
            numlook = inumlook
            output = ioutput
            verbose = iverbose

        print ' '


        # Display list of AGN
        if truth_menu == 'L':
            print selectdata.more()
            truth_nav = raw_input('BACK: ')
            truth_menu = 0

        # Further Plotting
        if truth_menu == 'P':
            
            # Plotting Options
            truth_plot = 0
            while truth_plot not in ('S','M','end','exit','q'):
            
                print 'PRINTING OPTIONS:'
                print '   [S] Sample Properties'
                print '   [M] Lightcurve Model'
                print '   [A] ... Automatic'
                print ' '

                truth_plot = raw_input('CHOICE: ')
                print ' ' 
                
                plt.ion()
                
                # --- SAMPLE PLOT ---
                if truth_plot == 'S':

                    # Create Figure 1
                    fig1 = plt.figure( figsize = ( 11, 6 ) )
                    fig1.suptitle( '%i | SDSS Quasars' % len( selectdata ) )
                    nbins = [50, 50]

                    gi_model = astrotools.accretion_color( wav['g'], wav['i'] )

                    # ax1: Redshift-Colour
                    ax11 = fig1.add_subplot( 221 )
                    ax11.hist2d( selectdata['reds'], selectdata['g'],
                                 range = [ [ izlo, izhi ], [ igbrt, igfnt ] ],
                                 bins = nbins )
                    ax11.set_ylim( igfnt, igbrt )
                    ax11.set_xlim( izlo, izhi )
                    ax11.set_xlabel( 'redshift z' )
                    ax11.set_ylabel( 'AB(g) (mag)' )

                    # ax2 : Colour-Magnitude diagram
                    ax12 = fig1.add_subplot( 222 )
                    ax12.hist2d( selectdata['g'] - selectdata['i'], selectdata['g'],
                                 range = [ [ icblu, icred ], [ igbrt, igfnt ] ],
                                 bins = nbins )
                    ax12.vlines( gi_model, igfnt, igbrt, colors='r', linestyles='--')
                    ax12.set_ylim( igfnt, igbrt )
                    ax12.set_xlim( icblu, icred )
                    ax12.set_xlabel( 'AB(g-i) (mag)' )
                    ax12.set_ylabel( 'AB(g) (mag)' )

                    # ax3 : Redshift-Colour diagram
                    ax13 = fig1.add_subplot( 223 )
                    ax13.hist2d( selectdata['reds'], selectdata['g'] - selectdata['i'],
                                 range = [ [ izlo, izhi ], [ icblu, icred ] ],
                                 bins = nbins )
                    ax13.hlines( gi_model, izlo, izhi, colors='r', linestyles='--' )
                    ax13.set_ylim( icred, icblu )
                    ax13.set_xlim( izlo, izhi ) 
                    ax13.set_ylabel( 'AB(g-i) (mag)' )
                    ax13.set_xlabel( 'redshift z' )

                    # ax4 : RA-Dec diagram
                    ax14 = fig1.add_subplot( 224 )
                    ax14.hist2d( selectdata['ra'], selectdata['dec'],
                                 range = [ [ iralo, irahi ], [ ideclo, idechi ] ],
                                 bins = nbins )
                    ax14.set_ylim( ideclo, idechi )
                    ax14.set_xlim( iralo, irahi )
                    ax14.set_ylabel( 'Dec (deg)' )
                    ax14.set_xlabel( 'RA (h)' )

                    # Must specifty layout POST plot
                    fig1.tight_layout(rect=[0,0,1,0.9])

                    # Draw figure
                    fig1.canvas.draw_idle()
                    truth_plot = ''

                    # Close figure
                    truth_plot = raw_input('CLOSE: ')
                    if truth_plot == '':
                      plt.close('all')
                      print ' '


