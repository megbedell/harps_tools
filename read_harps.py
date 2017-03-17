import numpy as np
from astropy.io import fits
import glob
import shk

def read_spec(spec_file):
    '''Read a HARPS 1D spectrum file from the ESO pipeline

    Parameters
    ----------
    spec_file : string
    name of the fits file with the data (s1d format)
    
    Returns
    -------
    wave : np.ndarray
    wavelength (in Angstroms)
    flux : np.ndarray
    flux value
    '''
    sp = fits.open(spec_file)
    header = sp[0].header
    n_wave = header['NAXIS1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    index = np.arange(n_wave, dtype=np.float64)
    wave = crval1 + index*cdelt1
    flux = sp[0].data
    return wave, flux
    
def read_spec_2d(spec_file, blaze=False, flat=False):
    '''Read a HARPS 2D spectrum file from the ESO pipeline

    Parameters
    ----------
    spec_file : string
    name of the fits file with the data (e2ds format)
    blaze : boolean
    if True, then divide out the blaze function from flux
    flat : boolean
    if True, then divide out the flatfield from flux
    
    Returns
    -------
    wave : np.ndarray (shape n_orders x 4096)
    wavelength (in Angstroms)
    flux : np.ndarray (shape n_orders x 4096)
    flux value 
    '''
    path = spec_file[0:str.rfind(spec_file,'/')+1]
    sp = fits.open(spec_file)
    header = sp[0].header
    flux = sp[0].data
    wave_file = header['HIERARCH ESO DRS CAL TH FILE']
    try:
        ww = fits.open(path+wave_file)
        wave = ww[0].data
    except:
        print "Wavelength solution file {0} not found!".format(wave_file)
        return
    if blaze:
        blaze_file = header['HIERARCH ESO DRS BLAZE FILE']
        bl = fits.open(path+blaze_file)
        blaze = bl[0].data
        flux /= blaze
    if flat:
        flat_file = header['HIERARCH ESO DRS CAL FLAT FILE']
        fl = fits.open(path+flat_file)
        flat = fl[0].data
        flux /= flat
    return wave, flux
    

def read_ccfs(filename):
    '''Read a HARPS CCF file from the ESO pipeline

    Parameters
    ----------
    filename : string
    name of the fits file with the data

    Returns
    -------
    velocity : np.ndarray
    velocity (in km/s)
    ccf : np.ndarray
    ccf value
    rv : float
    the pipeline-delivered RV **without drift correction applied** (in km/s)
    '''
    sp = fits.open(filename)
    header = sp[0].header

    # get the relevant header info
    n_vels = header['NAXIS1']
    n_orders = header['NAXIS2']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    
    # construct an (n_orders, n_vels) array of velocities
    index = np.arange(n_vels)
    velocity_one = crval1 + index*cdelt1
    velocity = np.tile(velocity_one, (n_orders,1))

    # get the ccf values
    ccf = sp[0].data
    
    # get the header RV
    rv = header['HIERARCH ESO DRS CCF RV']
        
    return velocity, ccf, rv    

def read_wavepar(filename):
    '''Parse wavelength solution on a HARPS file from the ESO pipeline

    Parameters
    ----------
    filename : string
    name of the fits file with the data (can be ccf, e2ds, s1d)

    Returns
    -------
    wavepar : np.ndarray
    a 72 x 4 array of wavelength solution parameters
    '''
    sp = fits.open(filename)
    header = sp[0].header
    
    #n_orders = header['NAXIS2']
    n_orders = 72
    wavepar = np.arange(4*n_orders, dtype=np.float).reshape(n_orders,4)
    for i in np.nditer(wavepar, op_flags=['readwrite']):
        i[...] = header['HIERARCH ESO DRS CAL TH COEFF LL{0}'.format(str(int(i)))]
    return wavepar


def read_snr(filename):
    '''Parse SNR from header of a HARPS file from the ESO pipeline

    Parameters
    ----------
    filename : string
    name of the fits file with the data (can be ccf, e2ds, s1d)

    Returns
    -------
    snr : np.ndarray
    SNR values taken near the center of each order
    '''
    sp = fits.open(filename)
    header = sp[0].header
    
    #n_orders = header['NAXIS2']
    n_orders = 72
    snr = np.arange(n_orders, dtype=np.float)
    for i in np.nditer(snr, op_flags=['readwrite']):
        i[...] = header['HIERARCH ESO DRS SPE EXT SN{0}'.format(str(int(i)))]
    return snr


def headers(filelist, outfile='data.csv'):
    # reads relevant info from fits headers & saves as a csv
    # input is a list of all bis files
    date, obj, bjd, rv, e_rv = [], [], [], [], []
    exptime, airm, drift = [], [], []
    bis, fwhm, s_hk, e_s_hk = [], [], [], []
    for f in filelist:
        b = fits.open(f)
        header = b[0].header
        date = np.append(date, header['DATE'])
        obj = np.append(obj, header['OBJECT'])
        bjd = np.append(bjd, header['HIERARCH ESO DRS BJD'])
        rv = np.append(rv, header['HIERARCH ESO DRS CCF RVC'])
        e_rv = np.append(e_rv, header['HIERARCH ESO DRS DVRMS'])
        exptime = np.append(exptime, header['EXPTIME'])
        airm = np.append(airm, header['HIERARCH ESO TEL AIRM START'])
        drift = np.append(drift, header['HIERARCH ESO DRS DRIFT SPE RV'])
        bis = np.append(bis, header['HIERARCH ESO DRS BIS SPAN'])
        fwhm = np.append(fwhm, header['HIERARCH ESO DRS CCF FWHM'])
        # load up the spectrum & measure activity index:
        i = str.find(f, 'bis')
        f2 = f[:i]+'s1d_A.fits'
        sp = fits.open(f2)
        header = sp[0].header
        n_wave = header['NAXIS1']
        crval1 = header['CRVAL1']
        cdelt1 = header['CDELT1']
        index = np.arange(n_wave, dtype=np.float64)
        wave = crval1 + index*cdelt1
        flux = sp[0].data
        s_hk_i, e_s_hk_i = shk.calc_shk(wave, flux, rv[-1])
        s_hk = np.append(s_hk, s_hk_i)
        e_s_hk = np.append(e_s_hk, e_s_hk_i)
    np.savetxt(outfile, np.transpose([date, obj, bjd, rv, e_rv, exptime, airm, drift, bis, fwhm, s_hk, e_s_hk]),
            header='date, obj, bjd, rv, e_rv, exptime, airm, drift, bis, fwhm, s_hk, e_s_hk', delimiter=',',
            fmt='%s')
        
        
        

if __name__ == "__main__":
    #data_dir = '/Users/mbedell/Documents/Research/GJ876/'
    #filelist = glob.glob(data_dir+'HARPS*_bis_*_A.fits')
    #headers(filelist, outfile=data_dir+'gj876.csv')
    
    data_dir = '/Users/mbedell/Documents/Research/HARPSTwins/Data/Reduced/'
    filelist = np.append(glob.glob(data_dir+'*/HARPS*_bis_*_A.fits'), glob.glob(data_dir+'archive/*/HARPS*_bis_*_A.fits'))
    filelist = np.append(filelist, glob.glob(data_dir+'18Sco/*/HARPS*_bis_*_A.fits'))
    print "{0} files found. Reading them now...".format(len(filelist))
    outfile = '/Users/mbedell/Documents/Research/HARPSTwins/Results/all.csv'
    headers(filelist, outfile=outfile)
    print "Results saved to: {0}".format(outfile)