import numpy as np
from PyAstronomy.pyTiming.pyPeriod import *
from scipy.io import readsav
import matplotlib.pyplot as plt
import pdb


def phaseplot(datafile, p=None):

    ### set up all possible observing run start/ends:

    gandolfi1 = [2458218.5, 2458223.5] # Apr 10 - 15
    gandolfi2 = [2458248.5, 2458252.5] # May 10 - 14 2018
    


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
    gandolfi1_fold = gandolfi1 % p_max
    gandolfi2_fold = gandolfi2 % p_max
    

    for x in [gandolfi1_fold, gandolfi2_fold]:
    	if x[0] > x[1]:
    	    x[1] = p_max
        

    plt.errorbar(t_fold,data.rv,data.sig,fmt='o')
    '''
    plt.axvspan(brahm1_fold[0], brahm1_fold[1], color='#D46A6A', alpha=0.5, lw=0, label='Brahm, Sept') #pink/red
    plt.axvspan(brahm2_fold[0], brahm2_fold[1], color='#b23434', alpha=0.5, lw=0, label='Brahm, Nov')
    plt.axvspan(melendez1_fold[0], melendez1_fold[1], color='#FFCC00', alpha=0.5, lw=0, label='Melendez, July') #yellow/orange
    plt.axvspan(melendez2_fold[0], melendez2_fold[1], color='#FFA200', alpha=0.5, lw=0, label='Melendez, Oct')
    plt.axvspan(melendez3_fold[0], melendez3_fold[1], color='#FF7900', alpha=0.5, lw=0, label='Melendez, Jan')
    plt.axvspan(melendez4_fold[0], melendez4_fold[1], color='#FF4100', alpha=0.5, lw=0, label='Melendez, Mar')
    plt.axvspan(kuerster1_fold[0], kuerster1_fold[1], color='#55AAAA', alpha=0.5, lw=0, label='Kuerster, Aug') #teal
    '''
    plt.axvspan(gandolfi1_fold[0], gandolfi1_fold[1], color='#77bb77', alpha=0.5, lw=0, label='Gandolfi, Apr') #green
    plt.axvspan(gandolfi2_fold[0], gandolfi2_fold[1], color='#FFCC00', alpha=0.5, lw=0, label='Gandolfi, May') #yellow/orange
    
    plt.legend(loc='best', fancybox=True, framealpha=0.5)
    plt.title('phased to {0:.1f} days'.format(p_max))

    plt.show()
