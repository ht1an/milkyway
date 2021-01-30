"""
this includes the correction of the quantities such as the selection or the distance, rv
"""
def correct_system_offset(true_value, obs_value, fmt="median"):
    import numpy as np
    format = ["median", "mean"]
    if fmt in format:
        if fmt == format[0]:
            V_sig = np.median(true_value - obs_value, (16, 50, 84))
            return V_sig
        else:
            Vmean = np.mean(true_value - obs_value)
            std = np.std(true_value - obs_value)
            return np.array([Vmean-std, Vmean, Vmean+std])





