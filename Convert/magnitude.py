"""
 this is the code for magnitude conversion between different systems and also the calculation of the magnitude
 from flux
"""
from astropy.io import fits
from astropy.table import Table
import numpy as np
from scipy import interpolate
import os

def flux_to_mag(obs_wave, obs_flux, filter_wave, filter_efficiency):
    """
    obs_flux: flux of the spectrum in unit erg/s/cm**2/A
    filter_wave: the wavelength in unit A
    obs_waveL: the wavelength of the observed spectra in unit A
    """
    wave_conversion = interpolate.interp1d(filter_wave, filter_efficiency, kind="linear")
    obs_efficiency = wave_conversion(obs_wave)

    # calculate the effective wavelength
    mean_flux = np.sum(obs_efficiency*obs_flux)/np.sum(obs_efficiency)
    mean_wave = np.sum(filter_wave*filter_efficiency)/np.sum(filter_efficiency)
    AB_mag_jy = mean_flux / ((2.997925 * (10 ** 18)) / ((mean_wave * 1) ** 2)) / 10 ** (-23)
    AB_mag = -2.5 * (np.log10(AB_mag_jy) - 23) - 48.6
    return AB_mag

def calc_mag_for_one(fn,path_passband):
    list_bands = ["NUV", "u", "g", "r",  "i", "z", "y", "GU", "GV", "GI"]
    mags = []
    dt = fits.open(fn)
    data = dt[1].data
    flux = data["FLUX"]*1e-17      # unit in erg/s/cm**2/A
    wave = 10**(data["LOGLAM"])
    for band in list_bands:
        data_band = np.loadtxt(path_passband + band + "_throughput.txt", skiprows=1)
        mags.append(flux_to_mag(wave, flux, data_band[:, 0], data_band[:, 1]))
    return mags

def calc_mag_for_plate(fn,path_passband):
    import matplotlib.pyplot as plt
    list_bands = ["NUV", "u", "g", "r",  "i", "z", "y", "GU", "GV", "GI"]
    wave_0 = 3.5792
    wave_stp = 0.0001
    dt = fits.open(fn)
    data_flux = dt[0].data
    data = dt[5].data
    npixel = len(data_flux[0,:])
    print("pixels are ", npixel,'--------')
    nstar = len(data["RA"])
    wave = 10**(np.linspace(wave_0, wave_0+(npixel-1)*wave_stp, npixel))
    mags = np.zeros((nstar, len(list_bands)))
    plt.plot(wave, data_flux[9,:])
    plt.ylim(-10,50)
    plt.show()
    # for i_star in range(nstar):
    #     flux = data_flux[i_star,:]*1e-17
    #     mags_tmp = []
    #     for band in list_bands:
    #         data_band = np.loadtxt(path_passband + band + "_throughput.txt", skiprows=1)
    #         mags_tmp.append(flux_to_mag(wave, flux, data_band[:, 0], data_band[:, 1]))
    #     mags[i_star,:] = mags_tmp
    return mags

# ---------------------------------------------------------
# here below can be changed
#  this is the path for the passthrough data file
# choice for type_datafile: "star"  or  "plate"
# "star" means each data file contains the data of just one star
# "plate" mean each data file contains the data of a plate including many stars
type_datafile = "plate"
path_passband = "./passband/csst/"
# this is the path for the data files
# path_data = "/Users/htian/Documents/work/CSST_Science/ZhangYanxia/8/"
path_data = "/Users/htian/Documents/work/CSST_Science/ZhangYanxia/tt/"
filename_output = "mags.fits"
# -----------------------------------------------------------
# -----------------------------------------------------------
# obtain the data file list in the path
list_of_files = []
for root, dirs, files in os.walk(path_data):
    for file in files:
        list_of_files.append(file)

list_bands = ["NUV", "u", "g", "r",  "i", "z", "y", "GU", "GV", "GI"]
fn_mag = []
mag_NUV = []
mag_u = []
mag_g = []
mag_r = []
mag_i = []
mag_z = []
mag_y = []
mag_GU = []
mag_GV = []
mag_GI = []
# calculate the magnitudes and save to the fits file
if type_datafile=="star":
    for name in list_of_files:
        if name[-4:]=="fits":
            mags=calc_mag_for_one(path_data+name, path_passband)
            fn_mag.append(name)
            mag_NUV.append(mags[0])
            mag_u.append(mags[1])
            mag_g.append(mags[2])
            mag_r.append(mags[3])
            mag_i.append(mags[4])
            mag_z.append(mags[5])
            mag_y.append(mags[6])
            mag_GU.append(mags[7])
            mag_GV.append(mags[8])
            mag_GI.append(mags[9])
    output = Table()
    output["fn"] = np.array(fn_mag)
    output["mag_NUV"] = mag_NUV
    output["mag_u"] = mag_u
    output["mag_g"] = mag_g
    output["mag_r"] = mag_r
    output["mag_i"] = mag_i
    output["mag_z"] = mag_z
    output["mag_y"] = mag_y
    output["mag_GU"] = mag_GU
    output["mag_GV"] = mag_GV
    output["mag_GI"] = mag_GI
    output.write(filename_output, overwrite=True)
if type_datafile=="plate":
    for name in list_of_files:
        if name[-4:]=="fits":
            mags=calc_mag_for_plate(path_data+name, path_passband)
            output = Table()
            output["mag_NUV"] = mags[:,0]
            output["mag_u"] = mags[:,1]
            output["mag_g"] = mags[:,2]
            output["mag_r"] = mags[:,3]
            output["mag_i"] = mags[:,4]
            output["mag_z"] = mags[:,5]
            output["mag_y"] = mags[:,6]
            output["mag_GU"] = mags[:,7]
            output["mag_GV"] = mags[:,8]
            output["mag_GI"] = mags[:,9]
            output.write("mags_"+name, overwrite=True)

    #%%