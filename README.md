# harps_tools
some codes I use with HARPS data.

### about this repository
most of these codes date back to my phd and will be of limited value! 
please don't judge me for the idl.

### potentially useful tools
you might find the following scripts useful:

- `unpack_eso.sh` : bash script to go from a directory full of ESO archive tarballs to all of those FITS files in a single, flat directory. just run it in the directory.
- `read_harps.py` : some python functions to load up data and read in various header keywords from HARPS data. in particular, the `read_headers()` function (example call in `__main__`), which scrapes a list of HARPS FITS files and generates a csv file with the RV timeseries and other useful ancillaries from their headers, might be useful.
- `calc_shk.py` : Calculates the Ca II H&K index (S_HK) from a spectrum.
