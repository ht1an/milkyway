"""
this is the code for fiting the density profiles including:
the thin/thick disk:
the halo:
"""
import numpy as np

def calc_halo_q_mc(R, Z, ns=10001):
    # the parameter and calculation are from the results of
    # Xu et al MNRAS 473, 1244–1257 (2018)
    sr_max = np.sqrt(R**2+Z**2)   # this is the maximum radius
    # generate the mock equivalent radius for ecllips, sr = sqrt(R**2+(Z/q)**2)
    sr = np.linspace(0, sr_max, ns) 
    p1 = 0.887
    p2 = 5.15
    p3 = -7.418
    p4 = 315.79
    sq = (p1 * sr ** 2 + p2 * sr) / (sr ** 2 + p3 * sr + p4) # the equation in page 1251
    msr = np.sqrt(R ** 2 + Z ** 2 / (sq ** 2))               # the calculated equivalent radius
    diff = np.abs(sr - msr)
    ind_min = np.argmin(diff)
    return sq[ind_min], sr[ind_min], diff[ind_min]



def density_halo(R0, Z0, q_halo=None, ns = None):
    n = len(R0)
    if q_halo is None:
        q_halo = np.zeros_like(R0)
        if ns is None:
            for i in range(n):
                q_halo[i] = calc_halo_q_mc(R0[i], Z0[i])
        else:
            for i in range(n):
                q_halo[i] = calc_halo_q_mc(R0[i], Z0[i], ns=ns)
    return (R0 / (np.sqrt(R0 ** 2 + (Z0 / q_halo) ** 2))) ** 5

def density_disk(Z0, H0, nu0=1):
    return nu0 / (np.cosh((Z0) / (2 * H0))) ** 2

def model_DDH(x, R0, y, prange):
    #     obs data
    Z_obs = y[:, 0]         # the height to the disk
    rho_obs = y[:, 1]       # density corresponding to Z
    # parameters
    P_tot = x[0]            #    total nu
    #     thin disk part
    f_thin = x[1]          #    fraction of the thin disk
    H_thin = x[2]          #    scale height of the thin disk
    #     thick disk part
    f_thick = x[3]         #    fraction of the thick disk
    H_thick = x[4]         #    scale height of the thick disk
    # calculate q for halo
    q_halo = np.zeros_like(Z_obs)
    for i in range(len(Z_obs)):
        q_halo[i], _, _ = calc_halo_q_mc(R0, Z_obs[i])

    if (P_tot < prange[0, 0]) or (P_tot > prange[0, 1]) or \
            (f_thin < prange[1, 0]) or (f_thin > prange[1, 1]) or \
            (H_thin < prange[2, 0]) or (H_thin > prange[2, 1]) or \
            (f_thin < prange[3, 0]) or (f_thin > prange[3, 1]) or \
            (H_thick < prange[4, 0]) or (H_thick > prange[4, 1]) or (f_thin + f_thick > 1):
        return -np.inf
    else:
        P_thin_model = 10 ** P_tot * f_thin * density_disk(H_thin, Z_obs)
        P_thick_model = 10 ** P_tot * f_thick * density_disk(H_thick, Z_obs)
        P_halo_model = 10 ** P_tot * (1-f_thin-f_thick)*density_halo(R0,Z_obs,q_halo)
        #     the chi squre
        lnp = np.sum(-0.5 * (np.log(P_thin_model + P_thick_model + P_halo_model) - np.log(rho_obs)) ** 2)
        return lnp


def model_DH(x, R0, y, prange):
    #     obs data
    Z_obs = y[:, 0]         # the height to the disk
    rho_obs = y[:, 1]       # density corresponding to Z
    # parameters
    P_tot = x[0]            #    total nu
    #     disk part
    f_disk = x[1]          #    fraction of disk
    H_disk = x[2]          #    scale height of disk
    # calculate q for halo
    q_halo = np.zeros_like(Z_obs)
    for i in range(len(Z_obs)):
        q_halo[i], _, _ = calc_halo_q_mc(R0, Z_obs[i])

    if (P_tot < prange[0, 0]) or (P_tot > prange[0, 1]) or \
            (f_disk < prange[1, 0]) or (f_disk > prange[1, 1]) or \
            (H_disk < prange[2, 0]) or (H_disk > prange[2, 1]) :
        return -np.inf
    else:
        P_disk_model = 10 ** P_tot * f_disk * density_disk(H_disk, Z_obs)
        P_halo_model = 10 ** P_tot * (1-f_disk)*density_halo(R0,Z_obs,q_halo)
        #     the chi squre
        lnp = np.sum(-0.5 * (np.log(P_disk_model + P_halo_model) - np.log(rho_obs)) ** 2)
        return lnp


def model_DD(x, R0, y, prange):
    #     obs data
    Z_obs = y[:, 0]  # the height to the disk
    rho_obs = y[:, 1]  # density corresponding to Z
    # parameters
    P_tot = x[0]  # total nu
    #     thin disk part
    f_thin = x[1]  # fraction of the thin disk
    H_thin = x[2]  # scale height of the thin disk
    #     thick disk part
    f_thick = 1 - f_thin  # fraction of the thick disk
    H_thick = x[3]  # scale height of the thick disk

    if (P_tot < prange[0, 0]) or (P_tot > prange[0, 1]) or \
            (f_thin < prange[1, 0]) or (f_thin > prange[1, 1]) or \
            (H_thin < prange[2, 0]) or (H_thin > prange[2, 1]) or \
            (H_thick < prange[3, 0]) or (H_thick > prange[3, 1]) or (f_thin + f_thick > 1):
        return -np.inf
    else:
        P_thin_model = 10 ** P_tot * f_thin * density_disk(H_thin, Z_obs)
        P_thick_model = 10 ** P_tot * f_thick * density_disk(H_thick, Z_obs)
        #     the chi squre
        lnp = np.sum(-0.5 * (np.log(P_thin_model + P_thick_model) - np.log(rho_obs)) ** 2)
        return lnp


def fit_density(R0, y, p0, N, N_burn_in, func, nwalkers, Prange):
    import emcee
    ndim = p0.shape[1]
    sampler = emcee.EnsembleSampler(nwalkers,ndim, func, args=[R0, y, Prange])
    pos, prob, state = sampler.run_mcmc(p0, N_burn_in)
    sampler.reset()
    sampler.run_mcmc(pos, N)
    samples = sampler.chain[:, N_burn_in:, :].reshape((-1, ndim))
    popt = np.median(samples, axis=0)
    pcov = np.zeros((ndim, ndim))
    for i in range(ndim):
        for j in range(ndim):
            pcov[i, j] = (np.sum((samples[:, i] - popt[i]) *
                                 (samples[:, j] - popt[j]))) / len(samples)
    return popt, pcov, samples
