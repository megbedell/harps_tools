import numpy as np
from PyAstronomy.pyTiming.pyPeriod import *
from scipy.io import readsav
import matplotlib.pyplot as plt
import pdb


def phaseplot(datafile, p=None):

	### set up all possible observing run start/ends:
	brahm1 = [2457637.44, 2457641.95]  # sunset Sept 05, sunrise Sept 10
	brahm2 = [2457706.4, 2457708.9]  # Nov 13 - 16
	#brahm3 = [2457804.4, 2457804.9]  # Feb 19 - 20
    
	melendez1 = [2457587.38, 2457588.94] # sunset July 17, sunrise July 19
	melendez2 = [2457664.4, 2457665.9] # Oct 02 - 04
	melendez3 = [2457763.4, 2457764.9] # Jan 09 - 11
	melendez4 = [2457807.4, 2457809.9] # Feb 22 - 25
    
	kuerster1 = [2457609.39, 2457611.95]  # sunset Aug 08, sunrise Aug 11
    
	gandolfi1 = [2457681.4, 2457687.9] # Oct 20 - 29


	### read in the star and find the period peak:
	data = readsav(datafile)

	if (p == None):
		period_grid = np.logspace(0.3,3,num=10000)  # ~2 days - a couple years
		freq_grid = 1.0/period_grid

		lc = TimeSeries(data.date, data.rv, data.sig)
		lp = Gls(lc, ofac=10, hifac=1, freq=freq_grid)

		p_max = period_grid[np.argsort(lp.power)[-1]]
	else:
		p_max = np.float64(p)

	t_fold = data.date % p_max
	brahm1_fold = brahm1 % p_max
	brahm2_fold = brahm2 % p_max
	melendez1_fold = melendez1 % p_max
	melendez2_fold = melendez2 % p_max
	melendez3_fold = melendez3 % p_max
	melendez4_fold = melendez4 % p_max
	kuerster1_fold = kuerster1 % p_max
	gandolfi1_fold = gandolfi1 % p_max
    
	for x in (brahm1_fold,brahm2_fold,melendez1_fold,melendez2_fold,melendez3_fold,melendez4_fold,kuerster1_fold,gandolfi1_fold):
		if x[0] > x[1]:
		    x[1] = p_max
            

	plt.errorbar(t_fold,data.rv,data.sig,fmt='o')
	plt.axvspan(brahm1_fold[0], brahm1_fold[1], color='#D46A6A', alpha=0.5, lw=0, label='Brahm, Sept') #pink/red
	plt.axvspan(brahm2_fold[0], brahm2_fold[1], color='#b23434', alpha=0.5, lw=0, label='Brahm, Nov')
	plt.axvspan(melendez1_fold[0], melendez1_fold[1], color='#FFCC00', alpha=0.5, lw=0, label='Melendez, July') #yellow/orange
	plt.axvspan(melendez2_fold[0], melendez2_fold[1], color='#FFA200', alpha=0.5, lw=0, label='Melendez, Oct')
	plt.axvspan(melendez3_fold[0], melendez3_fold[1], color='#FF7900', alpha=0.5, lw=0, label='Melendez, Jan')
	plt.axvspan(melendez4_fold[0], melendez4_fold[1], color='#FF4100', alpha=0.5, lw=0, label='Melendez, Mar')
	plt.axvspan(kuerster1_fold[0], kuerster1_fold[1], color='#55AAAA', alpha=0.5, lw=0, label='Kuerster, Aug') #teal
	plt.axvspan(gandolfi1_fold[0], gandolfi1_fold[1], color='#77bb77', alpha=0.5, lw=0, label='Gandolfi, Oct') #green
	plt.legend(loc='best', fancybox=True, framealpha=0.5)
	plt.title('phased to {0:.1f} days'.format(p_max))

	plt.show()
