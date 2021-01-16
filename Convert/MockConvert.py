
def TrueToObs(X_true, X_err):
    import numpy as np
    return X_true + np.random.normal(0,1,len(X_true)) * X_err