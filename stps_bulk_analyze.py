import numpy as np
from scipy.io.idl import readsav
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming.pyPeriod import TimeSeries, Gls
import analysis
import pdb
import read_harps

#starlist = np.genfromtxt('/Users/mbedell/Documents/Career/Thesis/phd-thesis/figures/rv_results/stps_starlist.txt', dtype=None)
starlist = np.genfromtxt('/Users/mbedell/Documents/Career/Thesis/phd-thesis/figures/rv_results/harps_starlist.txt', dtype=None)

linear = ['HIP14501', 'HIP62039', 'HIP87769', 'HIP73241', 'HIP103983', 'HIP108158', 'HIP6407', 'HIP18844', 'HIP64150']

curve = ['HIP79578', 'HIP81746', 'HIP19911', 'HIP83276', 'HIP54582']

sine = ['HIP14614', 'HIP43297', 'HIP54102', 'HIP72043', 'HIP65708', 'HIP67620']

def bootstrap_pearson(x, y, n_trials=10000):
    trials = np.empty(n_trials)
    for i in range(n_trials):
        ind = np.random.choice(len(x), len(x)) # bootstrap indices
        trials[i] = pearsonr(x[ind],y[ind])[0] 
    sigma = np.abs(np.mean(trials))/np.std(trials)  # sigma from zero   
    return np.mean(trials), sigma
    
    

if __name__ == "__main__":
    
    dir = '/Users/mbedell/Documents/Research/HARPSTwins/Results/Bulk/'
    outfile = dir+'summary.csv'
    f = open(outfile, 'w')
    f.write('star, RMS, RMS_corrected, n_stps, n_harps, n_bad, snr, baseline (yr), logrhk, candidate, trend, activity correction(s), rotation period from age (d), rotation period from activity (d)\n')
	
    outfile2 = dir+'rv_bulk.tex'
    tbl = open(outfile2, 'w')
    
    rms_change = []
    rms_percentchange = []
    
    for n,starname in enumerate(starlist):
        print 'beginning star {0}'.format(starname)
        
        candflag, trend, shkflag, fwhmflag, bisflag = False, '', False, False, False
        
        s = analysis.Star(starname)
        s.mask_bad()
        s.calc_rotation() # calculate rotation period from isochrone age
        rotation_age = s.rotation  # save it
        s.calc_age()  # re-estimate age from logrhk
        s.calc_rotation()
        rotation_logrhk = s.rotation
        if starname == 'HIP79672':
            s.bin(t=2./24.) # 2-hour bins
        else:
            s.bin() # 20-minute bins
            
        n_harps = len(s.rv)
        n_stps = np.sum(np.char.find(s.files,'archive') == -1) # where files do not have "archive" in the path
        n_bad = np.sum(np.invert(s.mask))
        baseline = (np.max(s.date) - np.min(s.date))/365.
        
        # estimate the co-added SNR:
        s.snr = np.zeros(len(s.files))
        for i,file in enumerate(s.files):
            snrs = read_harps.read_snr(file)
            s.snr[i] = snrs[59]  # approx 600 nm - see read_harps.read_wavepar(file)
        est_snr = np.sqrt(np.sum(s.snr**2.)) # add in quadrature
        
        
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        s.plot_rv(ax=ax1)
        ax1.set_title(starname)

        

        if starname in linear:
            s.subtract_trend(function=analysis.linear)
            trend = 'linear'
        if starname in curve:
            s.subtract_trend(function=analysis.curve)
            trend = 'curve'
        if starname in sine:
            s.subtract_trend(function=analysis.sine)
            trend = 'sine'
        if starname == 'HIP30037':
            par0 = np.asarray([31.6, 4244., 0.3, 63. * np.pi/180., 2455841.3, -200.]) # [period, K, ecc, omega, tp, offset]
            s.subtract_trend(function=analysis.keplerian, par0=par0)
            trend = 'keplerian'
        rms0 = np.std(s.rv[s.mask])
        
        s.plot_rv(ax=ax1)
        
        activity_peaks = ''
        for attr in ['shk', 'fwhm', 'bis']:
            if (bootstrap_pearson(s.rv, getattr(s,attr))[1] >= 3): # 3-sigma significant correlation
                s.subtract_activity(attr, errors=False)  # TODO: add error propagation later
                exec(attr+'flag = True')
                fig2 = plt.figure()
                ax2 = fig2.add_subplot(111)
                peaks, faps = s.periodogram(y=getattr(s,attr), dy=getattr(s,attr+'_err'), plot=True, ax=ax2, return_peaks=2)
                '''''
                if len(activity_peaks) == 0:
                    activity_peaks += '{0:.2f}, {1:.2f}'.format(peaks[0], peaks[1])
                else:
                    activity_peaks += '; {0:.2f}, {1:.2f}'.format(peaks[0], peaks[1])
                '''
                if len(activity_peaks) == 0:
                    activity_peaks += '{0:.2f}, {1:.2f}'.format(peaks[0], peaks[1])
                else:
                    activity_peaks += ', {0:.1f}'.format(peaks[0])
                ax2.set_title(starname+' '+attr)
                fig2.savefig(dir+'fig/'+starname+'_'+attr+'.png')
                plt.close(fig2)
        
        if shkflag or fwhmflag or bisflag:
            rms1 = '{0:.1f}'.format(np.std(s.rv[s.mask]))
            rms_change = np.append(rms_change, np.std(s.rv[s.mask]) - rms0)
            rms_percentchange = np.append(rms_change, (np.std(s.rv[s.mask]) - rms0)/rms0)            
        else:
            rms1 = ''

        s.plot_rv(ax=ax1)
        fig1.savefig(dir+'fig/'+starname+'.png')
        plt.close(fig1)
        
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        peaks, faps = s.periodogram(ax=ax3, return_peaks=10, plot=True)
        
        p_ind = np.where(faps <= 0.01)[0]
        rv_peaks = ''
        if len(p_ind) > 0:
            rv_peaks += '{0:.1f}'.format(peaks[p_ind[0]])
            for i in range(1, min(len(p_ind),3)):  # max of 3 significant periods
                rv_peaks += ', {0:.1f}'.format(peaks[p_ind[i]])
        if len(p_ind) > 3:
            rv_peaks += ', ...'
                
        ax3.set_title(starname)
        fig3.savefig(dir+'fig/'+starname+'_periodogram.png')
        

        plt.close(fig3)
        
        activity = ''
        if shkflag: activity += 'shk'
        if fwhmflag:
            if len(activity) > 0:
                activity += ', fwhm'
            else:
                activity += 'fwhm'
        if bisflag:
            if len(activity) > 0:
                activity += ', bis'
            else:
                activity += 'bis'
                
        if len(activity) == 0:
            # mark as quiet
            tbl.write('{name} & {n_harps} & {trend} & {rms0:.1f} & {activity} & {rms1} & {activity_peaks} & {rotation_age:.1f} & \
            {rotation_logrhk:.1f} & {rv_peaks}  \\\ \n'.format(
                name=starname[3:], rms0=rms0, rms1=rms1, trend=trend, activity=activity,
                rotation_age=rotation_age, rotation_logrhk=rotation_logrhk, rv_peaks=rv_peaks, 
                activity_peaks=activity_peaks, n_harps=n_harps))
        else:
            tbl.write('{name} & {n_harps} & {trend} & {rms0:.1f} & {activity} & {rms1} & {activity_peaks} & {rotation_age:.1f} & \
            {rotation_logrhk:.1f} & {rv_peaks}  \\\ \n'.format(
                name=starname[3:], rms0=rms0, rms1=rms1, trend=trend, activity=activity,
                rotation_age=rotation_age, rotation_logrhk=rotation_logrhk, rv_peaks=rv_peaks, 
                activity_peaks=activity_peaks, n_harps=n_harps))
        
        f.write('{name}, {rms0:.1f}, {rms1}, {n_stps}, {n_harps}, {n_bad}, {snr:.1f}, {baseline:.1f}, {logrhk}, \
                {cand}, {trend}, {activity}, {rotation_age:.1f}, {rotation_logrhk:.1f}, \n'.format(name=starname, 
                rms0=rms0, rms1=rms1, cand=candflag, trend=trend, activity=activity,
                logrhk=np.mean(s.logrhk), 
                n_stps=n_stps, n_harps=n_harps, n_bad=n_bad, baseline=baseline, rotation_age=rotation_age, 
                rotation_logrhk=rotation_logrhk, snr=est_snr))
        


        '''''
        print('{name} & {n_harps} & {trend} & {rms0:.1f} & {activity} & {rms1} & {activity_peaks} & {rotation_age:.1f} & \
        {rotation_logrhk:.1f} & {rv_peaks}  \\\ \n'.format(
            name=starname[3:], rms0=rms0, rms1=rms1, trend=trend, activity=activity,
            rotation_age=rotation_age, rotation_logrhk=rotation_logrhk, rv_peaks=rv_peaks, 
            activity_peaks=activity_peaks, n_harps=n_harps))
        '''
        
    print 'done!'
    f.close()
    tbl.close()
    
    print "{0} stars have RMS reduced by an average of {1:.1f} m/s or {2:.0f} percent after activity correction".format(np.sum(rms_change < 0.0), \
            - np.mean(rms_change[rms_change < 0.0]), - np.mean(rms_percentchange[rms_percentchange < 0.0] * 100.))
    print "{0} stars have RMS increased by an average of {1:.1f} m/s or {2:.0f} percent after activity correction".format(np.sum(rms_change > 0.0), \
            - np.mean(rms_change[rms_change > 0.0]), np.mean(rms_percentchange[rms_percentchange > 0.0] * 100.))
    

        
        