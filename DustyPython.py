# DustyPyton
# For decomposing AGN lightcurves
# Created 23.10.17

# Imports
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn
import astro_funclib as astrotools
import adv_funclib as advtools
import agn_funclib as agntools
import dust_funclib as dusttools
import time

# Input Parameters
lightcurvedir = 'QSO_S82/'
masterdir = 'QSO_Master/'
masterfile = 'DB_QSO_S82.dat'
outputdir = 'DustyPyton/'
outputfile = 'DusyOutput.junk.cat'

fidwav = 2400. # AA

# Extra verbosity in lc fitting
verbose1 = False

# Data Associated Arrays
# Filters (in order)
filters = ( 'u', 'g', 'r', 'i', 'z' )
filters = np.array( filters )
# Dustlaws
dustlaws = ( 'smc', 'lmc', 'mw', 'agn' )
dustlaws = np.array( dustlaws )
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

# Load dustlaw interpolations
tk_smc = dusttools.smccubic()
tk_lmc = dusttools.lmccubic()

# Load Masterfile
masterdata = ascii.read( table = masterdir + masterfile )
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
iverbose = 'ON'

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

                truth_plot = raw_input( 'CHOICE: ' )
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
                    fig1.tight_layout( rect=[0,0,1,0.9] )

                    # Draw figure
                    fig1.canvas.draw_idle()
                    truth_plot = ''

                    # Close figure
                    truth_plot = raw_input( 'CLOSE: ' )
                    if truth_plot == '':
                      plt.close('all')
                      print ' '


                # Lightcurves
                if truth_plot in ( 'M', 'A' ):

                    # Remind user that output is off in automatic mode
                    if ( truth_plot in 'A' ) & ( output == 'OFF' ):
                        truth_output = raw_input( 'WARNING: OUTPUT IS OFF - CHANGE? [Y/n] ')
                        if truth_output in ( 'Y', 'y' ):
                            output = 'ON'
                            print '* OUTPUT SET TO "ON"'

                    # Remind user that verbose is on in automatic mode
                    if ( truth_plot in 'A' ) & ( verbose == 'ON' ):
                        truth_output = raw_input( 'WARNING: VERBOSE IS ON - CHANGE? [Y/n] ')
                        if truth_output in ( 'Y', 'y' ):
                            verbose = 'OFF'
                            print '* VERBOSE SET TO "OFF"'

                    # Translate verbosity
                    if verbose == 'ON':
                        verbose = True
                    elif verbose == 'OFF':
                        verbose = False     


                    # init figures
                    if truth_plot in 'M':
                        fig2, ( ( ax21, ax22 ),
                        	    ( ax23, ax24 ),
                        	    ( ax25, ax26 ) ) = plt.subplots(nrows=3 ,ncols=2, figsize = (10,10))
                        ax23t = ax23.twiny()
                        ax24t = ax24.twiny()
                    
                    # output file
                    if output == 'ON':
                        print 'make that table.'
                    	# Setup table now? columns later?


                    # Loop over LC
                    numtotal = len( selectdata )
                    for lc_index, num in enumerate( xrange( numtotal ) ):
                        
                        # Get master data info
                        lc_info = selectdata[lc_index]

                        # Open raw lightcurve
                        lc_id = lc_info['dbID']
                        colnames = [ [ i + '_mjd', i + '_mag', i + '_merr' ] for i in filters ]
                        colnames = sum( colnames, [] ) + [ 'RA', 'Dec' ]
                        rawmag = ascii.read( table = lightcurvedir + str( lc_id ),
                                             names = colnames )

                        # Update user
                        num = num + 1
                        pc = num / numtotal * 100.
                        print '%3.2f%% (%i/%i)' %( pc, num, numtotal )
                        if verbose:
                            print ''
                            print 'Read: Opened file %i containing %i epochs' %( lc_id, len( rawmag ) )
                            print '-------------------------------------------------------------------'

                        # Mask bad values
                        cleanmag = agntools.maskval( intable = rawmag,
                                                      val = -99.99,
                                                      test = 'eq',
                                                      str_pattern = '_mag',
                                                      verbose = verbose )

                        cleanmag = agntools.maskval( intable = cleanmag,
                                                      val = 30.0,
                                                      test = 'gt',
                                                      str_pattern = '_mag',
                                                      verbose = verbose )

                        cleanmag = agntools.maskval( intable = cleanmag,
                                                      val = 10.00,
                                                      test = 'lt',
                                                      str_pattern = '_mag',
                                                      verbose = verbose )

                        cleanmag = agntools.maskval( intable = cleanmag,
                                                      val = 1.0,
                                                      test = 'gt',
                                                      str_pattern = '_merr',
                                                      verbose = verbose )

                        # Correct for dust extinction
                        extmag_ref = lc_info['Au']
                        corrmag = astrotools.extcorr( cleanmag, extmag_ref, extmag_coeff )

                        # Create table of fluxes
                        fluxtable = Table()
                        if verbose:
                            print 'Fluxtable: Created table of fluxes'

                        avgepoch = np.mean([ corrmag[ i + '_mjd'] for i in filters ], 0)
                        fluxtable.add_column( MaskedColumn( data = avgepoch,
                                                            name = 'mjd' ) )

                        if verbose:
                            print '* Added MJD column'
                        for i in filters:
                            f, fsig = astrotools.ab2mJy( corrmag[ i + '_mag' ], corrmag[ i + '_merr'] )
                            fcol = MaskedColumn( data = f,
                                                 name = i + '_f' )
                            fsigcol = MaskedColumn( data = fsig,
                                                    name = i + '_fsig' )
                            fluxtable.add_columns( [ fcol, fsigcol ] ) 

                            if verbose:
                                avgf, avgfsig = advtools.optavg( fcol, fsigcol )
                                print '* Average %s flux ... %1.4f +/- %1.4f mJy' %( i, avgf, avgfsig )

                        fluxtable = fluxtable.filled(np.nan)

                        # Fit lightcurves
                        # Update user
                        if verbose:
                            print 'Lightcurve: Fitting for model -> A(w) + B(w) * L(t)'
                        
                        # inital values
                        dlen = len( filters )
                        xlen = len( fluxtable )

                        A = np.zeros( dlen )
                        B = np.zeros( dlen )
                        Asig = np.zeros( dlen )
                        Bsig = np.zeros( dlen )

                        L = np.ones( xlen )
                        Lsig = np.ones( xlen )
                        L0 = np.ones( dlen )

                        # Loop over filters for inital guesses for A(w), B(w)
                        if verbose1:
                            print 'Optimal Average: Initalising A(w), B(w) ...'
                        for i, j in enumerate( filters ):
                            x = fluxtable[ j +'_f' ]
                            xsig = fluxtable[ j + '_fsig' ]
                            A[i], Asig[i] = advtools.optavg( x, xsig )
                            B[i], Bsig[i] = np.sqrt( advtools.optavg( ( x - A[i] )**2, xsig ) )

                            if verbose1:
                                print '* A(%s) = %2.5f +/- %2.5f | B(%s) = %2.5f +/- %2.5f' %( j, A[i], Asig[i], j, B[i], Bsig[i] )

                        # Loop until convergence
                        loopnum = 0
                        concrit = 1.
                        conval = 1E-10
                        while (concrit > conval ):

                            # Add to loop:
                            loopnum += 1

                            # Update user                            
                            if verbose1:
                                print 'Loop %i | Criterion ... %2.1e > %2.1e' %( loopnum, concrit, conval )

                            # Store old L
                            oldL = np.mean( L )

                            # Scale to find L(t)
                            if verbose1:
                                print 'Optimal Scaling: Calculating L(t) ...'
                            for i, fluxrow in enumerate( fluxtable ):
                                t = np.array( [ fluxrow[ j +'_f' ] for j in filters  ] )
                                tsig = np.array( [ fluxrow[ j + '_fsig' ] for j in filters ] )
                                L[i], Lsig[i] = advtools.optscl( t - A, tsig, B )

                                if verbose1:
                                    mjd = fluxrow['mjd']
                                    print '* L(%.1f) = %2.1f +/- %2.1f' %( mjd, L[i], Lsig[i] )

                            # Normalise L(t)
                            if verbose1:
                                print 'Normalisation: Fixing L(t)'
                            avgL = advtools.optavg( L, Lsig )[0]
                            L = L - avgL
                            avgL_squ = advtools.optavg( L**2, Lsig )[0]
                            L = L / np.sqrt( avgL_squ )
                            if verbose1:
                                newavgL = advtools.optavg( L, Lsig )[0]
                                print '* <L> = %2.1f | <L^2> = %2.1f ... <L> = %2.1f' %( avgL, avgL_squ, newavgL )

                            # Solve for A(w), B(w)
                            if verbose1:
                                print ' Fitline: Calculating A(w), B(w) ...'
                            for i, j in enumerate( filters ):
                                x = fluxtable[ j +'_f' ]
                                xsig = fluxtable[ j + '_fsig' ]
                                A[i], B[i], Asig[i], Bsig[i], L0[i] = advtools.fitline( L, x, xsig )

                                if verbose1:
                                    print '* A(%s) = %2.5f +/- %2.5f | B(%s) = %2.5f +/- %2.5f' %( j, A[i], Asig[i], j, B[i], Bsig[i] )

                            # Compute convergence criterion
                            Ldiff = abs( np.mean( L ) - oldL )
                            avgLsig = np.mean( Lsig )
                            concrit = Ldiff / avgLsig


                        # Update user
                        if verbose:
                            print 'Lightcurve: Convergence succeeded on loop %i' %( loopnum )
                            print 'Lightcurve: Final Values for %i' %( lc_id ) 
                            for i, j in enumerate( filters ):
                                print '* A(%s) = %2.5f +/- %2.5f | B(%s) = %2.5f +/- %2.5f' %( j, A[i], Asig[i], j, B[i], Bsig[i] )

                        # Calculate Fzero
                        F0, filterzero = agntools.estzero( A, Asig, B, Bsig, L0, filters )

                        # Update user
                        if verbose:
                            print 'Lightcurve: Flux zeropoint found at L(t) = %3.2f in %s' %( F0, filterzero )

                        # Calculate component spectra
                        # Update user
                        if verbose:
                            print 'Lightcurve: Calculating component spectra'

                        # Create table
                        spectable = Table()

                        # Calculate variances
                        Avar = Asig**2
                        Bvar = Bsig**2

                        # Precalculations
                        Lmax = np.nanmax( L )
                        Lmin = np.nanmin( L )

                        # Mean spectrum
                        mean_spec = A
                        mean_col = Column( name = 'mean', data = mean_spec )

                        meansig_spec = Asig
                        meansig_col = Column( name = 'mean_sig', data = meansig_spec )

                        # RMS spectrum
                        rms_spec = B
                        rms_col = Column( name = 'rms', data = rms_spec )

                        rmssig_spec = Bsig
                        rmssig_col = Column( name = 'rms_sig', data = rmssig_spec )

                        # Galaxy spectrum
                        gal_spec = A + B * ( F0 - L0 )
                        gal_col = Column( name = 'gal', data = gal_spec )

                        galsig_spec = np.sqrt( Avar + Bvar * ( F0 - L0 )**2 )
                        galsig_col = Column( name = 'gal_sig', data = galsig_spec )

                        # Bright spectrum
                        brt_spec = A + B * ( Lmax - L0 )
                        brt_col = Column( name = 'brt', data = brt_spec )

                        brtsig_spec = np.sqrt( Avar + Bvar * ( Lmax - L0 )**2 )
                        brtsig_col = Column( name = 'brt_sig', data = galsig_spec )

                        # Faint spectrum
                        fnt_spec = A + B * ( Lmin - L0 )
                        fnt_col = Column( name = 'fnt', data = fnt_spec )

                        fntsig_spec = np.sqrt( Avar + Bvar * ( Lmin - L0 )**2 )
                        fntsig_col = Column( name = 'fnt_sig', data = fntsig_spec )

                        # Difference spectrum
                        diff_spec = B * ( Lmax - Lmin )
                        diff_col = Column( name = 'diff', data = diff_spec )

                        diffsig_spec = np.sqrt( Bvar * ( Lmax - Lmin )**2 )
                        diffsig_col = Column( name = 'diff_sig', data = diffsig_spec )

                        # Variable spectrum
                        var_spec = A + B * ( L0 - F0 )
                        var_col = Column( name = 'var', data = var_spec )

                        varsig_spec = np.sqrt( Avar + Bvar * ( L0 - F0 ) )
                        varsig_col = Column( name = 'var_sig', data = varsig_spec )

                        # Disk spectrum
                        disk_spec = B * ( ( brt_spec - A ) / B - F0 )
                        disk_col = Column( name = 'disk', data = disk_spec )

                        disksig_spec = np.sqrt( Bvar ) * ( ( brt_spec - A ) / B - F0 )
                        disksig_col = Column( name = 'disk_sig', data = disksig_spec )

                        # Append to table
                        spectable.add_columns( [ mean_col, meansig_col,
                                                 rms_col, rmssig_col,
                                                 gal_col, galsig_col,
                                                 brt_col, brtsig_col,
                                                 fnt_col, fntsig_col,
                                                 diff_col, diffsig_col,
                                                 var_col, varsig_col,
                                                 disk_col, disksig_col ] )

                        # Update user
                        if verbose:
                            for i in spectable.colnames[::2]:
                                print '* %s spectrum ' %i
                                for j, k in enumerate( filters ):
                                    val = spectable[ i ][ j ]
                                    valsig = spectable[ i + '_sig'][ j ]
                                    print '** F(%s) = %5.5f +/- %5.5f mJy' %( k, val, valsig )


                        # Fit Dust Models
                        # Update user
                        if verbose:
                            print 'Dustfit: Beginning dust fitting'

                        # Calculate restwavelengths
                        z = lc_info['reds']
                        wavz = astrotools.obs2restwav( wav.values(), z )
                        wavz = { k : v for k, v in zip( wav.keys(), wavz ) }

                        # Update user
                        if verbose:
                            print 'Dustfit: Found rest wavelengths'
                            print '* Redshift = %3.2f' %z
                            for k in filters:
                                print '* %s observed light at %5.2f A' %(k, wavz[k])


                        # Loop over dustlaws
                        for dust in dustlaws:
                            if dust == 'smc':
                                extmag = dusttools.dustlaw( wavz.values(), tk_smc, 'smc' )
                            if dust == 'lmc':
                                extmag = dusttools.dustlaw( wavz.values(), tk_lmc, 'lmc' )
                            if dust == 'mw':
                                extmag = dusttools.dustlaw( wavz.values(), 'mw' )
                            if dust == 'agn':
                                extmag = dusttools.dustlaw( wavz.values(), 'agn' )

                            # Calculate residuals
                            #residuals = 




                        if truth_plot == 'M':
                            raw_input( '...' )