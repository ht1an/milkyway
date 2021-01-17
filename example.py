#%%
import sys
sys.path.append('./Dynamics')
import Orbits as O
sys.path.append('./Convert')
import MockConvert as MC
import numpy as np
from astropy import units
#%%
ppath = "/Users/htian/Desktop/"
# dist_tag = "parallax"
# dist_tag = "rMedPhotogeo"
# dist_tag = "rMedGeo"
dist_tag = None
Gabs_list = [-5.417, -5.623]
Age_list = [0.01360, 0.01260]
I_ind = 0
arrow_x = 25
arrow_y = 25
Gabs = Gabs_list[I_ind]
Age = Age_list[I_ind]
#%%
from astropy.io import fits
fn = "CWY.fits"
dt = fits.open(fn)
data = dt[1].data
ra = data["ra"]
dec = data["dec"]
if dist_tag is None:
    Gmag = data["phot_g_mean_mag"]
    ll = data["l"][0]
    bb = data["b"][0]
    print(ll,bb)
    d_mag, mabsmag = MC.Mags_to_distance(Gabs, Gmag, 2.489, 1, 100, ll, bb)
    dist = np.array([d_mag])
    print(d_mag, mabsmag, Gabs)
else:
    if dist_tag == "parallax":
        dist = 1/(data[dist_tag]+0.017)
    else:
        dist = data[dist_tag]/1000
# dist = data["rMedPhotogeo"]/1000
# dist = data["rMedGeo"]/1000
rv = 153
mu_ra = data["pmra"]
mu_dec = data["pmdec"]
vx, vy, vz = 0.,230,0. #0, 1, 0
max_age = 5
t = np.linspace(0,max_age,10001)
tb = np.linspace(0,-max_age,10001)
ind_mint = np.argmin(np.abs(tb+Age))
Ort, ts = O.Calc_Orbit(ra[0],dec[0], dist[0], mu_ra[0], mu_dec[0], rv ,t=t)
Ortb, tsb = O.Calc_Orbit(ra[0],dec[0], dist[0], mu_ra[0], mu_dec[0], rv ,t=tb)
#%%
import matplotlib.pyplot as plt
fig_xy = plt.figure(figsize=(8,6))
plt.scatter(Ort.x(ts), Ort.y(ts),s=5, c = t, cmap='jet',vmin=-max_age,vmax=max_age)
plt.scatter(Ortb.x(tsb), Ortb.y(tsb),s=5, c = tb, cmap='jet',vmin=-max_age,vmax=max_age)
plt.colorbar()
# plt.vlines(x=Ortb.x(tsb[ind_mint]),ymin=-50,ymax=50,color='gray')
# plt.hlines(y=Ortb.y(tsb[ind_mint]),xmin=-50,xmax=50,color='gray')
plt.plot(Ort.x(), Ort.y(),'kp',markersize=5)
if dist_tag is None:
    plt.arrow(Ortb.x(tsb[ind_mint])-arrow_x, Ortb.y(tsb[ind_mint])-arrow_y,arrow_x,arrow_y,color="gray",
              head_width=10,length_includes_head=True,width=3,
              label=f"(x, y)=({Ortb.x(tsb[ind_mint]):.2f},{Ortb.y(tsb[ind_mint]):.2f})")
    plt.legend(fontsize=15)
else:
    plt.title(f"Distance from {dist_tag}")
plt.xlabel("X(kpc)", fontsize=15)
plt.ylabel("Y(kpc)", fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
if dist_tag is None:
    plt.savefig(ppath + f"XY_Gabs{Gabs:.2f}.pdf")
else:
    plt.savefig(ppath+f"XY_{dist_tag}.pdf")
plt.show()

fig_xz = plt.figure(figsize=(8,6))
plt.scatter(Ort.x(ts), Ort.z(ts),s=5, c = t, cmap='jet',vmin=-max_age,vmax=max_age)
plt.scatter(Ortb.x(tsb), Ortb.z(tsb),s=5, c = tb, cmap='jet',vmin=-max_age,vmax=max_age)
plt.colorbar()
# plt.vlines(x=Ortb.x(tsb[ind_mint]),ymin=-50,ymax=50,color='gray')
# plt.hlines(y=Ortb.z(tsb[ind_mint]),xmin=-50,xmax=50,color='gray')
plt.plot(Ort.x(), Ort.z(),'kp',markersize=5)
if dist_tag is None:
    plt.arrow(Ortb.x(tsb[ind_mint])-arrow_x, Ortb.z(tsb[ind_mint])-arrow_y,arrow_x,arrow_y,color="gray",
              head_width=10,length_includes_head=True,width=3,
              label=f"(x, z)=({Ortb.x(tsb[ind_mint]):.2f},{Ortb.z(tsb[ind_mint]):.2f})")
    plt.legend(fontsize=15)
else:
    plt.title(f"Distance from {dist_tag}")
plt.xlabel("X(kpc)", fontsize=15)
plt.ylabel("Z(kpc)", fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
if dist_tag is None:
    plt.savefig(ppath + f"XZ_Gabs{Gabs:.2f}.pdf")
else:
    plt.savefig(ppath+f"XZ_{dist_tag}.pdf")
plt.show()


print(Ortb.z(tsb))