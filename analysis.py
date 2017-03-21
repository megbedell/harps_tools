import numpy as np
from scipy.io.idl import readsav
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
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
            
        def subtract_linear(self):
            # subtract a linear trend
            date,rv,sig = self.date,self.rv,self.sig
            slope0 = (rv[-1] - rv[0])/(date[-1] - date[0])
            int0 = rv[0] - slope0*date[0]
            par0 = np.asarray([int0, slope0])
            soln = leastsq(linear_resid, par0, args=(date, rv, sig))
            par = soln[0]
            self.rv -= linear(par, date)
            self.linpar = par
            
        def subtract_shk(s):
            # subtract a linear relation with S_HK
            par0 = np.asarray([0.0, 1.0])
            soln = leastsq(linear_resid, par0, args=(self.shk, self.rv, self.sig))
            par = soln[0]
            self.rv -= linear(par, self.shk)
            self.shkpar = par
            
   
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
    
    