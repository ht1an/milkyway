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
max_age = 0.02    # in unit gyr
arrow_x = 10
arrow_y = 10
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
    print(ll,bb,Gmag)
    d_mag, mabsmag, ebv = MC.Mags_to_distance(Gabs, Gmag, 2.489, 1, 100,10001, ll, bb)
    dist = np.array([d_mag])
    print(f"dist: {d_mag}, absmag:{mabsmag}, Gabs:{Gabs}, ebv:{ebv}",)
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
t = np.linspace(0,max_age,10001)
tb = np.linspace(0,-max_age,10001)
ind_mint = np.argmin(np.abs(tb+Age))
Ort, ts = O.Calc_Orbit(ra[0],dec[0], dist[0], mu_ra[0], mu_dec[0], rv ,t=t,PV_Sun=[11.1, 12.24, 7.25,8,220])
Ortb, tsb = O.Calc_Orbit(ra[0],dec[0], dist[0], mu_ra[0], mu_dec[0], rv ,t=tb,PV_Sun=[11.1, 12.24, 7.25,8,220])
#%%
import matplotlib.pyplot as plt
fig_xy = plt.figure(figsize=(8,6))
plt.scatter(Ort.x(ts), Ort.y(ts),s=5, c = t*1000, cmap='jet',vmin=-max_age*1000,vmax=max_age*1000)
plt.scatter(Ortb.x(tsb), Ortb.y(tsb),s=5, c = tb*1000, cmap='jet',vmin=-max_age*1000,vmax=max_age*1000)
plt.colorbar()
# plt.vlines(x=Ortb.x(tsb[ind_mint]),ymin=-50,ymax=50,color='gray')
# plt.hlines(y=Ortb.y(tsb[ind_mint]),xmin=-50,xmax=50,color='gray')
plt.plot(Ort.x(), Ort.y(),'kp',markersize=5)
if dist_tag is None:
    plt.arrow(Ortb.x(tsb)[ind_mint]-arrow_x, Ortb.y(tsb)[ind_mint]-arrow_y,arrow_x,arrow_y,color="gray",
              head_width=4,length_includes_head=True,width=1 ,
              label=f"(x, y)=({Ortb.x(tsb)[ind_mint]:.2f},{Ortb.y(tsb)[ind_mint]:.2f})")
    plt.legend(fontsize=15)
else:
    plt.title(f"Distance from {dist_tag}")
plt.xlabel("X(kpc)", fontsize=15)
plt.ylabel("Y(kpc)", fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim((0,20))
plt.ylim((-20,20))
if dist_tag is None:
    plt.savefig(ppath + f"XY_Gabs{Gabs:.2f}.pdf")
else:
    plt.savefig(ppath+f"XY_{dist_tag}.pdf")
plt.show()

# fig_xz = plt.figure(figsize=(8,6))
# plt.fill_between([-100,100],[5,5],[-5,-5],color="gray",alpha=0.2)
# plt.scatter(Ort.x(ts), Ort.z(ts),s=5, c = t*1000, cmap='jet',vmin=-max_age*1000,vmax=max_age*1000)
# plt.scatter(Ortb.x(tsb), Ortb.z(tsb),s=5, c = tb*1000, cmap='jet',vmin=-max_age*1000,vmax=max_age*1000)
# plt.colorbar()
# # plt.vlines(x=Ortb.x(tsb[ind_mint]),ymin=-50,ymax=50,color='gray')
# # plt.hlines(y=Ortb.z(tsb[ind_mint]),xmin=-50,xmax=50,color='gray')
# plt.plot(Ort.x(), Ort.z(),'kp',markersize=5)
# if dist_tag is None:
#     plt.arrow(Ortb.x(tsb)[ind_mint]-arrow_x, Ortb.z(tsb)[ind_mint]-arrow_y,arrow_x,arrow_y,color="gray",
#               head_width=4,length_includes_head=True,width=1,
#               label=f"(x, z)=({Ortb.x(tsb)[ind_mint]:.2f},{Ortb.z(tsb)[ind_mint]:.2f})")
#     plt.legend(fontsize=15)
# else:
#     plt.title(f"Distance from {dist_tag}")
# plt.xlabel("X(kpc)", fontsize=15)
# plt.ylabel("Z(kpc)", fontsize=15)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.xlim((0,20))
# plt.ylim((-20,20))
# if dist_tag is None:
#     plt.savefig(ppath + f"XZ_Gabs{Gabs:.2f}.pdf")
# else:
#     plt.savefig(ppath+f"XZ_{dist_tag}.pdf")
# plt.show()
# #
# #
# # fig_lb = plt.figure(figsize=(8,6))
# # plt.scatter(Ort.ll(ts), Ort.bb(ts),s=5, c = t*1000, cmap='jet',vmin=-max_age*1000,vmax=max_age*1000)
# # plt.scatter(Ortb.ll(tsb), Ortb.bb(tsb),s=5, c = tb*1000, cmap='jet',vmin=-max_age*1000,vmax=max_age*1000)
# # plt.colorbar()
# # plt.plot(Ort.ll(), Ort.bb(),'kp',markersize=5)
# # plt.xlabel("l(degree)",fontsize=15)
# # plt.ylabel("b(degree)",fontsize=15)
# # plt.xticks(fontsize=14)
# # plt.yticks(fontsize=14)
# # if dist_tag is None:
# #     plt.savefig(ppath + f"llbb_Gabs{Gabs:.2f}.pdf")
# # else:
# #     plt.savefig(ppath+f"llbb_{dist_tag}.pdf")
# # plt.show()
# #
# # fig_pm = plt.figure(figsize=(8,6))
# # plt.scatter(Ort.pmll(ts), Ort.pmbb(ts),s=5, c = t*1000, cmap='jet',vmin=-max_age*1000,vmax=max_age*1000)
# # plt.scatter(Ortb.pmll(tsb), Ortb.pmbb(tsb),s=5, c = tb*1000, cmap='jet',vmin=-max_age*1000,vmax=max_age*1000)
# # plt.colorbar()
# # plt.plot(Ort.pmll(), Ort.pmbb(),'kp',markersize=5)
# # plt.xlabel("$\\mu_l$(mas yr$^{-1}$)",fontsize=15)
# # plt.ylabel("$\\mu_b$(mas yr$^{-1}$)",fontsize=15)
# # plt.xticks(fontsize=14)
# # plt.yticks(fontsize=14)
# # if dist_tag is None:
# #     plt.savefig(ppath + f"pm_Gabs{Gabs:.2f}.pdf")
# # else:
# #     plt.savefig(ppath+f"pm_{dist_tag}.pdf")
# # plt.show()
# #
# # # fig_lb = plt.figure(figsize=(8,6))
# # # plt.quiver(Ort.ll(ts), Ort.bb(ts),Ort.pmll(ts), Ort.pmbb(ts),scale=50)
# # # plt.quiver(Ortb.ll(tsb), Ortb.bb(tsb),Ortb.pmll(tsb), Ortb.pmbb(tsb),scale=50)
# # # plt.scatter(Ort.ll(ts), Ort.bb(ts),s=5, c = t*1000, cmap='jet',vmin=-max_age*1000,vmax=max_age*1000)
# # # plt.scatter(Ortb.ll(tsb), Ortb.bb(tsb),s=5, c = tb*1000, cmap='jet',vmin=-max_age*1000,vmax=max_age*1000)
# # # plt.colorbar()
# # # plt.plot(Ort.ll(), Ort.bb(),'kp',markersize=5)
# # # plt.xlabel("l(degree)",fontsize=15)
# # # plt.ylabel("b(degree)",fontsize=15)
# # # plt.xticks(fontsize=14)
# # # plt.yticks(fontsize=14)
# # # if dist_tag is None:
# # #     plt.savefig(ppath + f"llbbpm_Gabs{Gabs:.2f}.pdf")
# # # else:
# # #     plt.savefig(ppath+f"llbbpm_{dist_tag}.pdf")
# # # plt.show()
# #
# # print(f"x:{Ortb.x(0)},\n",
# #       f"y:{Ortb.y(0)},\n",
# #       f"z:{Ortb.z(0)},\n",
# #       f"vx:{Ortb.vx(0)},\n",
# #       f"vy:{Ortb.vy(0)},\n",
# #       f"vz:{Ortb.vz(0)},\n",
# #       "v_tot:",np.sqrt(Ortb.vx(0)**2 + Ortb.vy(0)**2 + Ortb.vz(0)**2))