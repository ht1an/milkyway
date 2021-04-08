

def calc_ebv(ra, dec, dist, radec=True, kpc=True, degree=True):
    from dustmaps.bayestar import BayestarQuery
    bayestar = BayestarQuery(version='bayestar2019')
    from astropy.coordinates import SkyCoord
    import astropy.units as units
    import numpy as np
    if not kpc:
        dist = dist/1000
    if radec :
        import galpy.util.bovy_coords as gub
        llbb = gub.radec_to_lb(ra, dec, degree=degree)
        ll = llbb[:, 0]
        bb = llbb[:, 1]
    else:
        ll = np.copy(ra)
        bb = np.copy(dec)
    coords = SkyCoord(ll * units.deg, bb * units.deg, \
                          distance=dist * units.kpc, frame='galactic')
    ebv = bayestar(coords, mode='median')
    return ebv