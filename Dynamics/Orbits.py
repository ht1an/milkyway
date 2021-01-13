
# class SingleStar
def Calc_Orbit(px, py, pz, vx, vy, vz,t=None, PV_Sun=None, Potential=None):
    import numpy as np
    from astropy import units
    # from galpy.actionAngle import actionAngleStaeckel
    from galpy.orbit import Orbit
    # import galpy.util.bovy_conversion as bc
    from galpy.potential import MWPotential2014

    if PV_Sun == None:
        U_sun, V_sun, W_sun = 11.1, 12.24, 7.25
        X_sun = 8.3  # Reid et al 2014
        V_LSR = 232
    else:
        U_sun, V_sun, W_sun = PV_Sun[0:3]
        X_sun = PV_Sun[3]
        V_LSR = PV_Sun[4]
    if Potential == None:
        pot = MWPotential2014
    else:
        pot = Potential

    o = Orbit([px/X_sun, py/X_sun, pz/X_sun, vx/V_LSR, vy/V_LSR, vz/V_LSR],ro=X_sun,vo=V_LSR)
    if t == None:
        ts = np.linspace(0,10,1001)*units.Gyr
    else:
        ts = t*units.Gyr
    o.integrate(ts, pot=pot)
    return o

px, py, pz = 10, 0, 0
vx, vy, vz = 230, 0, 0
Ort = Calc_Orbit(px,py, pz, vx, vy, vz )
import matplotlib.pyplot as plt
plt.plot(Ort.x, Ort.y,'k.')
plt.show()