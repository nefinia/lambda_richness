import numpy as np
from math import *


def degtohms(ra, dec):
    ra_h = int(ra / 15)
    ra_m = int((ra / 15. - ra_h) * 60.)

    ra_s = ((ra / 15. - ra_h) * 60. - ra_m) * 60.
    dec_d = int(dec)
    dec_m = abs(int((dec - dec_d) * 60.))
    dec_s = abs((abs((dec - dec_d) * 60.) - dec_m) * 60.)

    ra_hms = [ra_h, ra_m, ra_s]
    dec_dms = [dec_d, dec_m, dec_s]

    print 'ra', ra, '->', ra_hms, 'dec', dec, '->', dec_dms
    return ra_hms, dec_dms


def gaussian(z, u, s): return np.exp(-(z - u) ** 2 / (2. * s ** 2)) / (sqrt(2 * np.pi) * s)


def pdfw(x, w=1, p=1.e-4):
    """Weighted probability density function.
       Created by Sofia G. Gallego, May 29 2015

       Usage
       =====
       input:
         x: value points
         w: weight of the values
         p: calculation precision

       output:
         x where the maximum of the weighted pdf is.
       """

    x = np.array(x)
    if type(w) == type(1):
        w = np.zeros(len(x)) + 1.
    else:
        w = np.array(w)
    steps = (max(x) * 1.1 - min(x) * .9) * (np.arange(int(1 / p))) * float(p) + min(x) * .9
    values = []

    corr = 1  # abs(steps[0]-steps[1])/(max(w)-min(w))/sqrt(p)

    for i in steps:
        values.append(sum(gaussian(i, x, corr / w)))

    mv = np.where(values == max(values))[0]
    pdfw = steps[mv][0]

    print 'found max pdf =', max(values), 'in x =', pdfw
    return pdfw  # ,steps,values,corr
