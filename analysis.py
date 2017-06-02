import numpy as np
from scipy.io.idl import readsav
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming.pyPeriod import TimeSeries, Gls
import shk
from utils import linefit_2derror as lf2d

class Star:
        def __init__(self, name, data_dir="/Users/mbedell/Documents/Research/HARPSTwins/Results/"):
            self.name = name
            s = readsav(data_dir+name+'_result.dat')
            self.rv = (np.asarray(s.rv) - np.mean(s.rv)) * 1.e3
            self.sig = np.asarray(s.sig) * 1.e3
            self.mask = np.ones(len(self.rv), dtype=bool)
            self.bis_err = np.asarray(s.sig)  # change these
            self.fwhm_err = np.asarray(s.sig) # change these
            for a in ['date','shk','bis','fwhm','shk_err']:
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
            
        def plot_periodogram(self, logpmin=-1, logpmax=4, y=None, dy=None, ax=None):
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
            
        def bin(self, t=1./3./24.):
            # bin observations by intervals of t days (default 20 minutes)                        
            if any(s.date != np.sort(s.date)):
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
                        for (x,dx) in [(self.rv, self.sig), (self.fwhm, self.fwhm_err), 
                                        (self.bis, self.bis_err), (self.shk, self.shk_err)]:
                            if np.sum(m) > 1: # average over good elements
                                x[i] = np.average(x[ind[m]], weights=dx[ind[m]])
                                dx[i] = np.std(x[ind[m]])/np.sqrt(len(ind[m])-1.)  # is this correct??
                            else: # just keep the good element
                                x[i] = x[ind[m]]
                                dx[i] = dx[ind[m]]
                        self.mask[i] = True
                    delete = np.append(delete, ind[1:])
                    i += len(ind)
                else:
                    i += 1
            delete = np.asarray(delete, dtype=int) # idk why it was float before??
            for attr in ['date', 'rv', 'sig', 'mask', 'bis', 'bis_err', 'fwhm', 'fwhm_err', 'shk', 'shk_err']:
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
            
        def subtract_linear(self):
            # subtract a linear trend
            par0 = calc_line(self.date[0],self.date[-1],self.rv[0],self.rv[-1]) # first guess parameters
            soln = leastsq(linear_resid, par0, args=(self.date[self.mask], self.rv[self.mask], self.sig[self.mask]))
            par = soln[0]
            self.rv -= linear(par, self.date)
            self.linpar = par
            
        def subtract_activity(self, marker_name, errors=True):
            # fit & subtract a linear relation between RV & activity marker marker_name
            valid = ['shk', 'bis', 'fwhm']
            if marker_name not in valid:
                print "Activity index {0} not recognized. Acceptable values: {1}".format(marker_name, valid)
            marker = getattr(self, marker_name)
            marker_err = getattr(self, marker_name+'_err')
            # fit:
            i0, i1 = np.argmax(marker), np.argmin(marker)
            par0 = calc_line(marker[i0],self.rv[i0],marker[i1],self.rv[i1]) # first guess parameters
            par = lf2d.bestfit(self.rv[self.mask], self.sig[self.mask], marker[self.mask], marker_err[self.mask], par0=par0)
            # subtract:
            self.rv -= linear(par, marker)
            setattr(self, marker_name+'par', par)
            if errors:
                par_err = lf2d.emcee_error(self.rv[self.mask], self.sig[self.mask], marker[self.mask], marker_err[self.mask], par0=par)
                # TODO: inflate errors on sig
                setattr(self, marker_name+'par_err', par_err)
            
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

def linear_resid(par,x,y,yerr):
    model = linear(par,x)
    return (y - model)/yerr     
  

    
if __name__ == "__main__":
    
    s = Star('HIP22263')
    s.mask_bad()
    #s.plot_rv()
    s.bin()
    #s.subtract_linear()
    #s.plot_rv(ax=plt.gca())
    
    #s.plot_periodogram()
    
    plt.errorbar(s.shk, s.rv, yerr=s.sig, xerr=s.shk_err, fmt='o')
    s.subtract_activity('shk', errors=False)
    plt.errorbar(s.shk, s.rv, yerr=s.sig, xerr=s.shk_err, fmt='o')
    
    
    