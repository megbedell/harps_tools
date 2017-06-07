import numpy as np
from scipy.io.idl import readsav
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming.pyPeriod import TimeSeries, Gls
import shk
from utils import linefit_2derror as lf2d
import kepler

def calc_line(x0, x1, y0, y1):
    """" 
    from two data points (x0, y0), (x1, y1), calculate line parameters
    returns: [intercept, slope] as np array
    """
    slope = (y1 - y0)/(x1 - x0)
    intercept = y0 - slope*x0
    return np.asarray([slope, intercept])
   
def linear(par,x):
    return par[0]*x + par[1]

def resid(par,fn,x,y,yerr):
    model = fn(par,x)
    return (y - model)/yerr 
    
def curve(par,x):
    return par[0]*x**2 + par[1]*x + par[2]
     
def sine(par,x):
    return par[0] * np.sin(2.*np.pi*x/par[1]) + par[2]
    
def keplerian(par,x):
    return kepler.calc_rvs(x, par)

class Star:
        def __init__(self, name, data_dir="/Users/mbedell/Documents/Research/HARPSTwins/Results/"):
            self.name = name
            try:
                s = readsav(data_dir+name+'_result.dat')
            except:
                print "Star {0} not found!".format(name)
                return
            self.rv = (np.asarray(s.rv) - np.mean(s.rv)) * 1.e3
            self.sig = np.asarray(s.sig) * 1.e3
            self.mask = np.ones(len(self.rv), dtype=bool)
            self.bis_err = np.asarray(s.sig)  # change these
            self.fwhm_err = np.asarray(s.sig) # change these
            self.files = np.asarray(s.files, dtype=np.string_)
            for a in ['date','shk','bis','fwhm','shk_err','airm','exp','logrhk']:
                setattr(self, a, np.asarray(getattr(s,a)))
                
        def __repr__(self):
            return "Star object with name = {0}, n_data = {1}".format(self.name, len(self.rv))
            
        def plot_rv(self, ax=None):
            # RV time-series plot
            if ax is None:
                fig,ax = plt.subplots(1,1)
            ax.errorbar(self.date[self.mask], self.rv[self.mask], self.sig[self.mask], fmt='o')
            inv = np.invert(self.mask)
            ax.errorbar(self.date[inv], self.rv[inv], self.sig[inv], fmt='o', color='red')
            ax.set_xlabel('BJD')
            ax.set_ylabel(r'RV (m s$^{-1}$)')
            
        def plot_periodogram(self, logpmin=-1, logpmax=4, y=None, dy=None, return_peaks=0, max_fap=False, ax=None):
            # calculate periodogram
            if y is None:
                y = self.rv
            if dy is None:
                dy = self.sig
            period_grid = np.logspace(logpmin, logpmax, num=10000)
            freq_grid = 1.0/period_grid
            lc = TimeSeries(self.date[self.mask], y[self.mask], dy[self.mask])
            lp = Gls(lc, ofac=10, hifac=1, freq=freq_grid)
            if ax is None:
                fig,ax = plt.subplots(1,1)
            ax.plot(period_grid, lp.power)
            ax.axhline(lp.powerLevel(0.1))
            ax.axhline(lp.powerLevel(0.01))
            ax.set_xscale('log')
            ax.set_xlabel('Period (d)')
            ax.set_ylabel('Power') 
            if return_peaks > 0:
                if max_fap:
                    return period_grid[np.argsort(lp.power)[-return_peaks:]], lp.stats(np.max(lp.power))['FAP']
                else:
                    return period_grid[np.argsort(lp.power)[-return_peaks:]]
            
        def bin(self, t=1./3./24.):
            # bin observations by intervals of t days (default 20 minutes)                        
            if any(self.date != np.sort(self.date)):
                print "data not in order, cannot bin"  # TODO: add sorting capability
                return
            delete = []
            i=0
            while (i < len(self.date)):
                j = (self.date[i:] - self.date[i]) <= t
                if np.sum(j) > 1:
                    ind = np.where(j)[0] + i  # indices in this time bin
                    m = self.mask[ind]  # a mini-mask
                    if any(m):  # if there are non-masked elements in the bin...
                        self.date[i] = np.median(self.date[ind[m]])
                        self.exp[i] = np.sqrt(np.sum(self.exp[ind[m]]**2))
                        for (x,dx) in [(self.rv, self.sig), (self.fwhm, self.fwhm_err), 
                                        (self.bis, self.bis_err), (self.shk, self.shk_err)]:
                            if np.sum(m) > 1: # average over good elements
                                x[i] = np.average(x[ind[m]], weights=dx[ind[m]])
                                dx[i] = np.std(x[ind[m]])/np.sqrt(len(ind[m])-1.)  # TODO: is this correct??
                            else: # just keep the good element
                                x[i] = x[ind[m]]
                                dx[i] = dx[ind[m]]
                        self.mask[i] = True
                    delete = np.append(delete, ind[1:])
                    i += len(ind)
                else:
                    i += 1
            delete = np.asarray(delete, dtype=int) # idk why it was float before??
            for attr in ['date', 'rv', 'sig', 'mask', 'bis', 'bis_err', 'fwhm', 'fwhm_err', 'shk', 'shk_err', 'airm', 'exp', 'logrhk', 'files']:
                setattr(self, attr, np.delete(getattr(self,attr), delete))
                    
            
        def mask_bad(self):
            # mask out all epochs where an activity index is 3 stdev away from average
            for marker in [self.shk, self.fwhm, self.bis]:
                resid = marker - np.mean(marker)
                bad = np.where(np.abs(resid) > 3.*np.std(marker))
                for b in bad[0]:
                    if b > -1:
                        self.mask[b] = 0 
                        #print "masked index {0}".format(b)
      
             
        def subtract_trend(self, function=linear, plot=False, par0=np.zeros(1)): 
            if not par0.any():   
                if function == linear:
                    par0 = calc_line(self.date[0],self.date[-1],self.rv[0],self.rv[-1]) # first guess parameters
                elif function == curve:
                    par0 = np.ones(3)
                elif function == sine:
                    # par0 = [K, period, offset]
                    par0 = np.asarray([1.e2, 1.e4, np.mean(self.rv)])
                elif function == keplerian:
                    # par0 = [period, K, ecc, omega, M0, offset]
                    par0 = np.asarray([100., 20., 0.0, 0.0, 0.0, np.mean(self.rv)])
            if function not in [linear, curve, sine, keplerian]:
                print "Function not recognized. Acceptable values: linear, curve, sine, keplerian"
                return
            soln = leastsq(resid, par0, args=(function, self.date[self.mask], self.rv[self.mask], self.sig[self.mask]))
            par = soln[0]
            if plot:
                fig,ax = plt.subplots(1,1)
                self.plot_rv(ax=ax)
                xs = np.arange(min(self.date), max(self.date), 1.)
                ax.plot(xs, function(par, xs))
            self.rv -= function(par, self.date)
            self.trendtype = str(function)
            self.trendpar = par
                             
        def subtract_activity(self, marker_name, errors=True):
            # fit & subtract a linear relation between RV & activity marker marker_name
            valid = ['shk', 'bis', 'fwhm']
            if marker_name not in valid:
                print "Activity index {0} not recognized. Acceptable values: {1}".format(marker_name, valid)
            marker = getattr(self, marker_name)
            marker_err = getattr(self, marker_name+'_err')
            # fit:
            #i0, i1 = np.argmax(marker[self.mask]), np.argmin(marker[self.mask])
            #par0 = calc_line(marker[i0],marker[i1],self.rv[i0],self.rv[i1]) # first guess parameters
            par0 = np.zeros(2)
            par = lf2d.bestfit(self.rv[self.mask], self.sig[self.mask], marker[self.mask], marker_err[self.mask], par0=par0)
            # subtract:
            self.rv -= linear(par, marker)
            setattr(self, marker_name+'par', par)
            if errors:
                par_err = lf2d.emcee_error(self.rv[self.mask], self.sig[self.mask], marker[self.mask], marker_err[self.mask], par0=par)
                # TODO: inflate errors on sig
                setattr(self, marker_name+'par_err', par_err)
                
    
if __name__ == "__main__":
    '''''
    s = Star('HIP14614')
    s.mask_bad()
    s.bin()
    s.subtract_trend(function=sine, plot=True)
    
    s2 = Star('HIP87769')
    s2.mask_bad()
    s2.bin()
    s2.subtract_trend(function=linear, plot=True)

    s3 = Star('HIP101905')
    s3.mask_bad()
    s3.bin()
    print is_correlated(s3.rv, s3.shk)
    plt.errorbar(s3.shk, s3.rv, yerr=s3.sig, xerr=s3.shk_err, fmt='o')
    s3.subtract_activity('shk', errors=False)
    xs = np.arange(min(s3.shk), max(s3.shk), 0.01)
    plt.plot(xs, linear(s3.shkpar, xs))
    plt.errorbar(s3.shk, s3.rv, yerr=s3.sig, xerr=s3.shk_err, fmt='o')
    
    plt.errorbar(s3.bis, s3.rv, yerr=s3.sig, xerr=s3.bis_err, fmt='o')
    print pearsonr(s3.rv, s3.bis)[1]
    s3.subtract_activity('bis', errors=False)
    xs = np.arange(min(s3.bis), max(s3.bis), 0.01)
    plt.plot(xs, linear(s3.bispar, xs))
    plt.errorbar(s3.bis, s3.rv, yerr=s3.sig, xerr=s3.bis_err, fmt='o')
    
    
    s4 = Star('HIP30037')
    s4.mask_bad()
    s4.bin()
    par0 = np.asarray([31.6, 4244., 0.3, 63. * np.pi/180., 2455841.3, -200.]) # [period, K, ecc, omega, tp, offset]
    #s4.plot_rv()
    #xs = np.arange(min(s4.date), max(s4.date), 0.1)
    #plt.plot(xs, kepler.calc_rvs(xs, par0))
    s4.subtract_trend(function=keplerian, plot=True, par0=par0)
    
    '''''
    
    
    s = Star('HIP7585')
    s.mask_bad()
    s.bin()

    
    plt.errorbar(s.fwhm, s.rv, yerr=s.sig, xerr=s.fwhm_err, fmt='o')
    #s.subtract_activity('fwhm', errors=False)
    soln = leastsq(resid, np.ones(2), args=(linear, s.fwhm[s.mask], s.rv[s.mask], s.sig[s.mask]))
    xs = np.arange(min(s.fwhm), max(s.fwhm), 0.01)
    #par = s.fwhmpar
    par = soln[0]
    plt.plot(xs, linear(par, xs))
    #plt.errorbar(s.fwhm, s.rv, yerr=s.sig, xerr=s.fwhm_err, fmt='o')
    
    #soln = leastsq(resid, par0, args=(function, self.date[self.mask], self.rv[self.mask], self.sig[self.mask]))

    
    
    
    #s.plot_rv(ax=plt.gca())    
    