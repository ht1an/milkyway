#%%
import sys
sys.path.append('./Dynamics')
import Orbits as O
import numpy as np
from astropy import units
#%%
from astropy.io import fits
fn = "CWY.fits"
dt = fits.open(fn)
data = dt[1].data
ra = data["ra"]
dec = data["dec"]
dist = 1/(data["parallax"]+0.027)
# dist = data["rMedPhotogeo"]/1000
# dist = data["rMedGeo"]/1000
rv = 153
mu_ra = data["pmra"]
mu_dec = data["pmdec"]
vx, vy, vz = 0.,230,0. #0, 1, 0
t = np.linspace(0,5,100001)
Ort = O.Calc_Orbit(ra[0],dec[0], dist[0], mu_ra[0], mu_dec[0], rv ,t=t)
#%%
import matplotlib.pyplot as plt
Ort.plot(d1="x",d2="y")
plt.plot(Ort.x(), Ort.y(),'k.')
plt.show()

Ort.plot(d1="x",d2="z")
plt.plot(Ort.x(), Ort.z(),'k.')
plt.show()