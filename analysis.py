import numpy as np
from scipy.io.idl import readsav
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming.pyPeriod import TimeSeries, Gls
import shk

class Star:
        def __init__(self, name, data_dir="/Users/mbedell/Documents/Research/HARPSTwins/Results/"):
            self.name = name
            s = readsav(data_dir+name+'_result.dat')
            self.rv = (s.rv - np.mean(s.rv)) * 1.e3
            self.sig = s.sig * 1.e3
            for a in ['date','shk','bis','fwhm','shk_err']:
                setattr(self, a, getattr(s,a))
                
        def __repr__(self):
            return "Star object with name = {0}, n_data = {1}".format(self.name, len(self.rv))
            
        def plot_rv(self, ax=None):
            # RV time-series plot
            if ax is None:
                fig,ax = plt.subplots(1,1)
            ax.errorbar(self.date, self.rv, self.sig, fmt='o')
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
            lc = TimeSeries(self.date, y, dy)
            lp = Gls(lc, ofac=10, hifac=1, freq=freq_grid)
            if ax is None:
                fig,ax = plt.subplots(1,1)
            ax.plot(period_grid, lp.power)
            ax.axhline(lp.powerLevel(0.1))
            ax.axhline(lp.powerLevel(0.01))
            ax.set_xscale('log')
            ax.set_xlabel('Period (d)')
            ax.set_ylabel('Power')            
            
        def subtract_linear(self):
            # subtract a linear trend
            par0 = calc_line(self.date[0],self.date[-1],self.rv[0],self.rv[-1]) # first guess parameters
            soln = leastsq(linear_resid, par0, args=(self.date, self.rv, self.sig))
            par = soln[0]
            self.rv -= linear(par, self.date)
            self.linpar = par
            
        def subtract_shk(self):
            # subtract a linear relation with S_HK
            i0, i1 = np.argmax(self.shk), np.argmin(self.shk)
            par0 = calc_line(self.shk[i0],self.shk[i1],self.rv[i0],self.rv[i1]) # first guess parameters
            soln = leastsq(linear_resid, par0, args=(self.shk, self.rv, self.sig))
            par = soln[0]
            self.rv -= linear(par, self.shk)
            self.shkpar = par
            
        def subtract_fwhm(self):
            # subtract a linear relation with FWHM
            i0, i1 = np.argmax(self.fwhm), np.argmin(self.fwhm)
            par0 = calc_line(self.fwhm[i0],self.fwhm[i1],self.rv[i0],self.rv[i1]) # first guess parameters
            soln = leastsq(linear_resid, par0, args=(self.fwhm, self.rv, self.sig))
            par = soln[0]
            self.rv -= linear(par, self.fwhm)
            self.fwhmpar = par
            
        def subtract_bis(self):
            # subtract a linear relation with BIS
            i0, i1 = np.argmax(self.bis), np.argmin(self.bis)
            par0 = calc_line(self.bis[i0],self.bis[i1],self.rv[i0],self.rv[i1]) # first guess parameters
            soln = leastsq(linear_resid, par0, args=(self.fwhm, self.rv, self.sig))
            par = soln[0]
            self.rv -= linear(par, self.bis)
            self.bispar = par
            
def calc_line(x0, x1, y0, y1):
    """" 
    from two data points (x0, y0), (x1, y1), calculate line parameters
    returns: [intercept, slope] as np array
    """
    slope = (y1 - y0)/(x1 - x0)
    intercept = y0 - slope*x0
    return np.asarray([intercept, slope])
   
def linear(par,x):
    return par[0] + par[1]*x

def linear_resid(par,x,y,yerr):
    model = linear(par,x)
    return (y - model)/yerr     
  

    
if __name__ == "__main__":
    
    s = Star('HIP14501')
    s.plot_rv()
    s.subtract_linear()
    s.plot_rv(ax=plt.gca())
    
    s.plot_periodogram()
    
    