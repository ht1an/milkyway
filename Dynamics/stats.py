"""
this is used to calculate the statistics
"""
import numpy as np


def beta_anisotropic(v_theta, v_phi, v_r):
    std_theta = np.nanstd(v_theta)
    std_phi = np.nanstd(v_phi)
    std_r = np.nanstd(v_r)
    return 1 - 0.5 * (std_theta**2 + std_phi**2)/(std_r**2)