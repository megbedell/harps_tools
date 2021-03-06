import numpy as np
from astropy.io import fits
import glob
import shk
from starchive import identifiers

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
    try:
        wave_file = header['HIERARCH ESO DRS CAL TH FILE']
    except KeyError: # HARPS-N
        wave_file = header['HIERARCH TNG DRS CAL TH FILE']
    wave_file = str.replace(wave_file, 'e2ds', 'wave') # just in case of header mistake..
                                                       # ex. HARPS.2013-03-13T09:20:00.346_ccf_M2_A.fits
    try:
        ww = fits.open(path+wave_file)
        wave = ww[0].data
    except:
        print("Wavelength solution file {0} not found!".format(wave_file))
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
    '''Parse wavelength solution on a HARPS file from the ESO pipeline

    Parameters
    ----------
    filename : string
    name of the fits file with the data (can be ccf, e2ds, s1d)

    Returns
    -------
    snr : np.ndarray
    M-order length array of SNR values taken near the center of each order
    '''
    sp = fits.open(filename)
    header = sp[0].header
    
    n_orders = 72
    snr = np.arange(n_orders, dtype=np.float)
    for i in np.nditer(snr, op_flags=['readwrite']):
        try:
            i[...] = header['HIERARCH ESO DRS SPE EXT SN{0}'.format(str(int(i)))]
        except: # if HARPS-N
            i[...] = header['HIERARCH TNG DRS SPE EXT SN{0}'.format(str(int(i)))]
    return snr


def headers(filelist, outfile='data.csv', offset=(0.0,0.0), shk=True):
    """
    Reads relevant RV-related information from FITS headers.
    
    Parameters
    ----------
    filelist : list of str
        N-length list of file names for BIS files. 
        (Note: CCF files will do for most of this info,
        but the BIS and FWHM parameters will fail.)
    outfile : str (default: 'data.csv')
        Name of output table.
    offset : (float, float) tuple (default: (0.0, 0.0))
        Offset of RVs post July 2015 HARPS upgrade relative to pre-upgrade, in (offset [m/s], 
        uncertainty [m/s]) form. Will be automatically subtracted off from the returned RVs.
    shk : bool (default : True)
        If True, load up the spectrum and measure the S_HK activity index.
        Note that this relies upon the S1D file being in the same folder as the BIS file
        (which it will be by default if you untarred all the ESO archive data products)
    
    Returns
    -------
    data : list of arrays
        List of various quantities that were saved in `outfile`. Each entry
        is an N-length array.
    """
    date, obj, bjd, rv, e_rv = [], [], [], [], []
    exptime, airm, drift, progid = [], [], [], []
    bis, fwhm, s_hk, e_s_hk, snr = [], [], [], [], []
    for f in filelist:
        b = fits.open(f)
        header = b[0].header
        date = np.append(date, header['DATE'])
        obj = np.append(obj, header['OBJECT'])
        bjd = np.append(bjd, header['HIERARCH ESO DRS BJD'])
        rv = np.append(rv, header['HIERARCH ESO DRS CCF RVC'] * 1e3) # m/s
        e_rv = np.append(e_rv, header['HIERARCH ESO DRS DVRMS'])
        exptime = np.append(exptime, header['EXPTIME'])
        airm = np.append(airm, header['HIERARCH ESO TEL AIRM START'])
        drift = np.append(drift, header['HIERARCH ESO DRS DRIFT SPE RV'])
        bis = np.append(bis, header['HIERARCH ESO DRS BIS SPAN'])
        fwhm = np.append(fwhm, header['HIERARCH ESO DRS CCF FWHM'])
        progid = np.append(progid, header['HIERARCH ESO OBS PROG ID'])
        snr = np.append(snr, header['HIERARCH ESO DRS SPE EXT SN65'])
        # adjust for instrument upgrade:
        if np.float(bjd[-1]) > 2457218.5:
            rv[-1] = rv[-1] - offset[0]
            e_rv[-1] = np.sqrt(e_rv[-1]**2 + offset[1]**2)
        if shk:
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
            s_hk_i, e_s_hk_i = shk.calc_shk(wave, flux, rv[-1] / 1e3)
            s_hk = np.append(s_hk, s_hk_i)
            e_s_hk = np.append(e_s_hk, e_s_hk_i)
        else:
            s_hk = np.zeros_like(rv) + np.nan
            e_s_hk = np.zeros_like(rv)
    np.savetxt(outfile, np.transpose([filelist, date, obj, bjd, rv, e_rv, exptime, progid, airm, drift, bis, fwhm, s_hk, e_s_hk, snr]),
            header='filename, date, obj, bjd, rv, e_rv, exptime, progid, airm, drift, bis, fwhm, s_hk, e_s_hk, snr', delimiter=',',
            fmt='%s')
    return [filelist, date, obj, bjd, rv, e_rv, exptime, progid, airm, drift, bis, fwhm, s_hk, e_s_hk, snr]
                
def single_star_data(hipname, header_data, conv=None):
    # takes a star name in format 'HIPx', applies cuts to header_data to select only that star
    [filelist, date, obj, bjd, rv, e_rv, exptime, progid, airm, drift, bis, fwhm, s_hk, e_s_hk] = header_data
    # selection mask
    if conv is None:
        conv = identifiers.Converter()
    hipno = int(hipname[3:])
    hdno = conv.hiptohd(hip_number)
    hdname = 'HD' + str(hdno)
    mask = [str.replace(o, '-', '') in [hipname, hdname] for o in obj]
    if hipname == 'HIP79672': # special case
        mask[obj == '18_Sco'] = True
        mask[obj == '18Sco'] = True
    elif hipname == 'HIP36512':
        mask[obj == 'HD59711A'] = True
    elif hipname == 'HD221287':
        mask[obj == 'HDTE221287'] = True
    elif hipname == 'HIP85042':
        mask[obj == 'HD157347_std'] = True
    elif hipname == 'HIP7585':
        mask[obj == 'HIP7505'] = True
    elif hipname == 'HIP36515':
        mask[obj == 'HD059967'] = True
    # select
    star_data = [0. for i in range(len(header_data))] # initialize
    for i,x in enumerate(header_data):
        star_data[i] = x[mask]
    return star_data
        
def write_systemic(bjd, rv, e_rv, starname='star', starmass=1.0, systemic_dir='/Applications/Systemic/datafiles/'):
    # takes lists/arrays of data for a star, writes out Systemic files
    vels_file = systemic_dir+starname+'.vels'
    with open(vels_file, 'w') as f:
        for i in range(len(bjd)):
            f.write('{0:16.8f}    {1:10.3f}    {2:5.3f}    \n'.format(bjd[i], rv[i], e_rv[i]))
    with open(systemic_dir+starname+'.sys', 'w') as f:
        f.write('Data  {{\n	RV[] "{vels_file}"\n}}\n"{starname}" {{    \n  Mass  {starmass} \n}}'.format())            
    
        
        
        

if __name__ == "__main__":
    
    #data_dir = '/Users/mbedell/Documents/Research/HARPSTwins/Data/Reduced/'
    #filelist = np.append(glob.glob(data_dir+'*/HARPS*_bis_*_A.fits'), glob.glob(data_dir+'archive/*/HARPS*_bis_*_A.fits'))
    #filelist = np.append(filelist, glob.glob(data_dir+'18Sco_archive/*/HARPS*_bis_*_A.fits'))
    #print "{0} files found. Reading them now...".format(len(filelist))
    #outfile = '/Users/mbedell/Documents/Research/HARPSTwins/Results/all.csv'
    #data = headers(filelist, outfile=outfile, , offset=(15.4,0.4))
    #print "Results saved to: {0}".format(outfile)
    
    data_dir = '/Users/mbedell/Documents/Research/pi-men/'
    filelist = glob.glob(data_dir+'HARPS*_bis_*_A.fits')
    print("{0} files found. Reading them now...".format(len(filelist)))
    outfile = '/Users/mbedell/Documents/Research/pi-men/harps_rvs.csv'
    data = headers(filelist, outfile=outfile, offset=(0.0,0.0))
    print("Results saved to: {0}".format(outfile))
    
    
    
    