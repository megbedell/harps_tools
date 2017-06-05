import numpy as np
from scipy.io.idl import readsav
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming.pyPeriod import TimeSeries, Gls
import analysis
import pdb

#starlist = np.genfromtxt('/Users/mbedell/Documents/Career/Thesis/phd-thesis/figures/rv_results/stps_starlist.txt', dtype=None)
starlist = np.genfromtxt('/Users/mbedell/Documents/Career/Thesis/phd-thesis/figures/rv_results/harps_starlist.txt', dtype=None)

linear = ['HIP14501', 'HIP62039', 'HIP87769']

curve = []

sine = ['HIP14614', 'HIP43297', 'HIP54102', 'HIP72043', 'HIP79578', 'HIP81746']

if __name__ == "__main__":
    
    dir = '/Users/mbedell/Documents/Research/HARPSTwins/Results/Bulk/'
    outfile = dir+'summary.csv'
    f = open(outfile, 'w')
    f.write('star, RMS, RMS_corrected, n_stps, n_harps, n_bad, baseline (yr), candflag, \
            linflag, curveflag, sineflag, shkflag, fwhmflag, bisflag, rv_peaks, activity_peaks\n')
	
    for starname in starlist:
        print 'beginning star {0}'.format(starname)
        
        candflag, linflag, curveflag, sineflag, shkflag, fwhmflag, bisflag = False, False, False, False, False, False, False
        
        s = analysis.Star(starname)
        s.mask_bad()
        if starname == 'HIP79672':
            s.bin(t=2./24.) # 2-hour bins
        else:
            s.bin() # 20-minute bins
            
        n_harps = len(s.rv)
        n_stps = np.sum(np.char.find(s.files,'archive') == -1) # where files do not have "archive" in the path
        n_bad = np.sum(np.invert(s.mask))
        baseline = (np.max(s.date) - np.min(s.date))/365.
        
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
        
        activity_peaks = []
        for attr in ['shk', 'fwhm', 'bis']:
            if (pearsonr(s.rv, getattr(s,attr))[1] < 0.003): # 3-sigma significant correlation
                s.subtract_activity(attr, errors=False)  # TODO: add error propagation later
                exec(attr+'flag = True')
                fig2 = plt.figure()
                ax2 = fig2.add_subplot(111)
                peaks = s.plot_periodogram(y=getattr(s,attr), dy=getattr(s,attr+'_err'), ax=ax2, return_peaks=2)
                activity_peaks = np.append(activity_peaks, peaks)
                ax2.set_title(starname+' '+attr)
                fig2.savefig(dir+'fig/'+starname+'_'+attr+'.png')
                plt.close(fig2)
        
        rms1 = np.std(s.rv[s.mask])

        s.plot_rv(ax=ax1)
        fig1.savefig(dir+'fig/'+starname+'.png')
        plt.close(fig1)
        
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        rv_peaks, fap = s.plot_periodogram(ax=ax3, return_peaks=5, max_fap=True)
        if fap <= 0.01:
            candflag = True
        ax3.set_title(starname)
        fig3.savefig(dir+'fig/'+starname+'_periodogram.png')
        plt.close(fig3)
        
        f.write('{name}, {rms0:.3f}, {rms1:.3f}, {n_stps}, {n_harps}, {n_bad}, {baseline:.1f}, \
                {cand}, {lin}, {curve}, {sine}, {shk}, {fwhm}, {bis}, {rv_peaks}, {activity_peaks}\n'.format(name=starname, 
                rms0=rms0, rms1=rms1, cand=candflag, lin=linflag, curve=curveflag, sine=sineflag, shk=shkflag, fwhm=fwhmflag,
                bis=bisflag, rv_peaks=np.asarray(rv_peaks, dtype='|S8'), activity_peaks=np.asarray(activity_peaks, dtype='|S8'), 
                n_stps=n_stps, n_harps=n_harps, n_bad=n_bad, baseline=baseline))
        
    print 'done!'
    f.close()
        