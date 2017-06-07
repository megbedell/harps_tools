import numpy as np
from numpy import log, exp, pi, sqrt, sin, cos, tan, arctan

def calc_ma(T0, t, period):
    # calculate mean anomaly
    days = t - T0
    phase = days/period % 1.0
    ma = phase * 2.0 * pi
    return ma
    
def calc_ea(ma, ecc):
    # calculate eccentric anomaly from mean anomaly, eccentricity
    tolerance = 1e-3
    ea = np.copy(ma)
    while True:
        diff = ea - ecc * sin(ea) - ma
        ea -= diff / (1. - ecc * cos(ea))
        if abs(diff).all() <= tolerance:
            break
    return ea
 
    
def calc_rvs(t,par):
    '''
    Calculate RV(t) given par
    par: [period, K, ecc, omega, M0, offset]
    where omega is the argument of periastron
    and Tp is time at periastron
    '''
    P,K,ecc,omega,tp,offset = par
    
    # enforce boundaries on parameters:
    #if (P < 0.0 or K < 0.0 or ecc < 0.0 or ecc > 0.999 or omega < 0. or omega > 2.*pi or M0 < 0. or M0 > 2.*pi):
    #    return np.zeros_like(t)
    #P = max([0.0, P])
    #K = max([0.0, K])
    #ecc = min([max([0.0, ecc]), 0.99])
    #omega = min([max([-pi, omega]), pi])
    #M0 = min([max([-pi, M0]), pi])
    
    ma = 2. * pi / P * (t - tp)  # mean anomaly
    ea = calc_ea(ma, ecc)  # eccentric anomaly

    f = 2.0 * np.arctan2(sqrt(1+ecc)*sin(ea/2.0), sqrt(1-ecc)*cos(ea/2.0)) # true anomaly
    rvs = - K * (cos(omega + f) + ecc*cos(omega))
    return rvs + offset


def calc_ipower(et,t,y,dy,period):
    ecc = et[0]
    T0 = et[1]
    ma = calc_ma(T0, t, period)  # mean anomaly
    ea = calc_ea(ma, ecc)  # eccentric anomaly from Kepler's equation
    v = arctan(2.0 * sqrt((1+ecc)/(1-ecc)) * tan(ea/2.0))  # true anomaly
   
    # now solve the equation RV = c + a*cos(v) + b*sin(v):
    w = 1.0/dy**2 / (1.0/dy**2).sum()  # normalized weights
    wmean = (w*y).sum()  # weighted mean of y
    YY = (w * (y-wmean)**2).sum()
    YC = (w*(y-wmean)*cos(v)).sum()
    YS = (w*(y-wmean)*sin(v)).sum()
    CC = (w*cos(v)**2).sum() - (w**2*cos(v)**2).sum()
    SS = (w*sin(v)**2).sum() - (w**2*sin(v)**2).sum()
    CS = (w*cos(v)*sin(v)).sum() - (w **2*cos(v)*sin(v)).sum()
    D = CC * SS - CS**2
    
    return ((YY * D)) * (SS*YC**2 + CC*YS**2 - 2.0*CS*YC*YS).sum()   # returns INVERSE POWER for minimization purposes
    
def calc_npower(et,t,y,dy,period):
    ecc = et[0]
    T0 = et[1]
    ma = calc_ma(T0, t, period)  # mean anomaly
    ea = calc_ea(ma, ecc)  # eccentric anomaly from Kepler's equation
    v = arctan(2.0 * sqrt((1+ecc)/(1-ecc)) * tan(ea/2.0))  # true anomaly
   
    # now solve the equation RV = c + a*cos(v) + b*sin(v):
    w = 1.0/dy**2 / (1.0/dy**2).sum()  # normalized weights
    wmean = (w*y).sum()  # weighted mean of y
    YY = (w * (y-wmean)**2).sum()
    YC = (w*(y-wmean)*cos(v)).sum()
    YS = (w*(y-wmean)*sin(v)).sum()
    CC = (w*cos(v)**2).sum() - (w**2*cos(v)**2).sum()
    SS = (w*sin(v)**2).sum() - (w**2*sin(v)**2).sum()
    CS = (w*cos(v)*sin(v)).sum() - (w **2*cos(v)*sin(v)).sum()
    D = CC * SS - CS**2
    
    return -1.0/((YY * D)) * (SS*YC**2 + CC*YS**2 - 2.0*CS*YC*YS).sum()   # returns NEGATIVE POWER for minimization purposes
