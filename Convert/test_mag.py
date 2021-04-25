from astropy.io import fits
import numpy as np
import magnitude as M
def calc_mag(fn,path_passband):
    # fn = "spec-6582-56274-0034.fits"
    dt = fits.open(fn)
    data = dt[1].data
    wave = 10**(data["LOGLAM"])
    flux = data["FLUX"]*1e-17
    list_bands = ["NUV", "u", "g", "r",  "i", "z", "y", "GU", "GV", "GI"]
    mags = []
    for band in list_bands:
        data_band = np.loadtxt(path_passband + band +"_throughput.txt",skiprows=1)
        mags.append(M.flux_to_mag(wave, flux, data_band[:,0], data_band[:,1]))
    return mags

path_passband = "./passband/csst/"
fn = "spec-6582-56274-0034.fits"
mags = calc_mag(fn, path_passband)
print(mags)