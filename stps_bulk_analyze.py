import numpy as np
from scipy.io.idl import readsav
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming.pyPeriod import TimeSeries, Gls
import analysis
import pdb

starlist = np.genfromtxt('/Users/mbedell/Documents/Career/Thesis/phd-thesis/figures/rv_results/stps_starlist.txt', dtype=None)

linear = ['HIP14501', 'HIP62039', 'HIP87769']

curve = []

sine = ['HIP14614', 'HIP43297', 'HIP54102', 'HIP72043', 'HIP79578', 'HIP81746']

if __name__ == "__main__":
    
    dir = '/Users/mbedell/Documents/Research/HARPSTwins/Results/Bulk/'
    outfile = dir+'summary.csv'
    f = open(outfile, 'w')
    f.write('star, RMS, RMS_corrected, linflag, curveflag, shkflag, fwhmflag, bisflag, \
                rv_peak1, rv_peak2, rv_peak3, rv_peak4, rv_peak5\n')
	
    for starname in starlist:
        print 'beginning star {0}'.format(starname)
        
        linflag, curveflag, sineflag, shkflag, fwhmflag, bisflag = False, False, False, False, False, False
        
        s = analysis.Star(starname)
        s.mask_bad()
        s.bin()
        
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        s.plot_rv(ax=ax1)
        ax1.set_title(starname)

        

        if starname in linear:
            s.subtract_trend(function=analysis.linear)
            linflag = True
        if starname in curve:
            s.subtract_trend(function=analysis.curve)
            curveflag = True
        if starname in sine:
            s.subtract_trend(function=analysis.sine)
            sineflag = True
        rms0 = np.std(s.rv[s.mask])
        
        s.plot_rv(ax=ax1)
        
        for attr in ['shk', 'fwhm', 'bis']:
            if (pearsonr(s.rv, getattr(s,attr))[1] < 0.003): # 3-sigma significant correlation
                s.subtract_activity(attr, errors=False)  # TODO: add error propagation later
                exec(attr+'flag = True')
                fig2 = plt.figure()
                ax2 = fig2.add_subplot(111)
                s.plot_periodogram(y=getattr(s,attr), dy=getattr(s,attr+'_err'), ax=ax2)
                ax2.set_title(starname+' '+attr)
                fig2.savefig(dir+'fig/'+starname+'_'+attr+'.png')
                plt.close(fig2)
        
        rms1 = np.std(s.rv[s.mask])

        s.plot_rv(ax=ax1)
        fig1.savefig(dir+'fig/'+starname+'.png')
        plt.close(fig1)
        
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        rv_peaks = s.plot_periodogram(ax=ax3, return_peaks=5)
        ax3.set_title(starname)
        fig3.savefig(dir+'fig/'+starname+'_periodogram.png')
        plt.close(fig3)
        
        f.write('{name}, {rms0:.3f}, {rms1:.3f}, {lin}, {curve}, {sine}, {shk}, {fwhm}, {bis}, {rv_peaks}\n'.format(name=starname, 
                    rms0=rms0, rms1=rms1, lin=linflag, curve=curveflag, sine=sineflag, shk=shkflag, fwhm=fwhmflag,
                    bis=bisflag, rv_peaks=rv_peaks))
        
    print 'done!'
    f.close()
        