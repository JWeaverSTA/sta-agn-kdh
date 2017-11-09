import numpy as np
from astropy.io import ascii

dat = ascii.read('Output/DustyOutput.trial.cat', format = 'fixed_width')

def poppc( stat, law ):
	return len(dat[dat[stat+'_'+law] == np.min( (dat[stat+'_smc'],dat[stat+'_lmc'],dat[stat+'_mw'],dat[stat+'_agn'] ), 0 ) ]) / float(len(dat))

print 'Percentage of objects providing minimal statistic'
for i in ('chisq', 'std', 'ebmv'):
	print '-------'
	for j in ('smc', 'lmc', 'mw', 'agn'):
		print '%2.1f%% ... %s(%s)' %(poppc( i, j) * 100, i, j )

