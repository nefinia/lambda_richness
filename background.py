#!/usr/bin/env python
preamble = """
####################################################
#
# Background galaxies estimator (for Rykoff+ 2012 lambda estimator)
# Autor: Sofia G. Gallego, August 28 2014 
# 
####################################################
"""

import numpy as np
from math import *


class bkg_data:
    # Create class from previously estimated background
    sigma_g = []
    color = []
    mag = []
    ncolorbins = 0
    colorbinsize = 0
    colorzero = 0.0
    nmagbins = 0
    magbinsize = 0
    magzero = 0.0

    def __init__(self, ncolorbins, colorbinsize, colorzero, nmagbins, magbinsize, magzero, data):
        self.ncolorbins = ncolorbins
        self.colorbinsize = colorbinsize
        self.colorzero = colorzero
        self.magzero = magzero
        self.nmagbins = nmagbins
        self.magbinsize = magbinsize
        self.color = np.arange(ncolorbins) * .1 + colorzero
        self.mag = np.arange(nmagbins) * magbinsize + magzero
        self.sigma_g = [[] for i in np.arange(nmagbins)]
        for i in range(nmagbins): self.sigma_g[i] = data[i]
        self.sigma_g = np.array(self.sigma_g)
        self.sigma_e = np.sqrt(self.sigma_g) / self.sigma_g


class background_cat:
    sigma_g = []
    sigma_e = []
    color = []
    mag = []
    ncolorbins = 0
    colorbinsize = 0
    colorzero = 0.0
    nmagbins = 0
    magbinsize = 0
    magzero = 0.0
    h_mag = []
    h_color = []

    def __init__(self, ncolorbins, colorbinsize, colorzero, nmagbins, magbinsize, magzero, color, mag, area,
                 smooth=True):
        self.ncolorbins = ncolorbins
        self.colorbinsize = colorbinsize
        self.colorzero = colorzero
        self.nmagbins = nmagbins
        self.magbinsize = magbinsize
        self.magzero = magzero
        self.color = np.arange(ncolorbins) * colorbinsize + colorzero
        self.mag = np.arange(nmagbins) * magbinsize + magzero
        self.sigma_g = [[] for i in np.arange(nmagbins)]

        colormin = colorzero
        colormax = colorzero + ncolorbins * colorbinsize
        magmin = magzero
        magmax = magzero + nmagbins * magbinsize

        if smooth:
            good = np.where((mag >= magmin) & (mag <= magmax) & (color >= colormin) & (color <= colormax))[0]
            print 'total', len(mag), 'good', len(good)
            mag = mag[good]
            color = color[good]
            weight = np.zeros(len(mag)) + 1.
            x = (color - colormin) * ncolorbins / (colormax - colormin)
            y = (mag - magmin) * nmagbins / (magmax - magmin)
            h2d = cic(weight, y, nmagbins, x, ncolorbins)

        else:
            h2d = \
                np.histogram2d(mag, color, bins=(nmagbins, ncolorbins), range=[[magmin, magmax], [colormin, colormax]])[
                    0]

        self.h_mag = np.histogram(mag, bins=nmagbins, range=[magmin, magmax])
        self.h_color = np.histogram(color, bins=ncolorbins, range=[colormin, colormax])

        self.sigma_g = h2d / (area * magbinsize * colorbinsize)
        self.sigma_g = np.array(self.sigma_g)
        self.sigma_e = np.sqrt(self.sigma_g) / self.sigma_g


def bkgplot(bkgs, colors=['g-r', 'r-i', 'i-z'], ns=[131, 132, 133], figsize=(18, 10), log=False, show=True, save=True,
            outname='bkg.png', area=False, h2d=None, ramax=None, ramin=None, decmin=None, decmax=None):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.ion()

    if area:
        extent = [ramax, ramin, decmin, decmax]
        plt.figure(figsize=(10, 10))
        plt.xlabel('RA (degrees)')
        plt.ylabel('DEC (degrees)')
        plt.imshow(h2d, extent=extent, interpolation='nearest')
        plt.colorbar(shrink=.9)

    else:
        plt.figure(figsize)
        for bkg, c, n in zip(bkgs, colors, ns):
            f = plt.subplot(n, axisbg='midnightblue', ylabel='i', xlabel=c)
            extent = [bkg.color[0], bkg.color[-1], bkg.mag[-1], bkg.mag[0]]
            print extent
            print 'max value', np.amax(bkg.sigma_g)
            if log:
                aa = f.imshow(bkg.sigma_g, norm=LogNorm(vmin=1, vmax=np.amax(bkg.sigma_g)), interpolation='nearest',
                              extent=extent)
            else:
                aa = f.imshow(bkg.sigma_g, interpolation='nearest', extent=extent)
            plt.xticks(np.arange(bkg.color[0], bkg.color[-1] + 1., 1.0))
        plt.colorbar(aa, label='Number counts')

    if show: plt.show()
    if save: plt.savefig(outname)
    plt.close()


def cic(value, x, nx, y=None, ny=1, wraparound=False):
    """ Interpolate an irregularly sampled field using Cloud in Cell
    method.

    This function interpolates an irregularly sampled field to a
    regular grid using Cloud In Cell (nearest grid point gets weight
    1-dngp, point on other side gets weight dngp, where dngp is the
    distance to the nearest grid point in units of the cell size).
    
    Inputs
    ------
    value: array, shape (N,)
        Sample weights (field values). For a temperature field this
        would be the temperature and the keyword average should be
        True. For a density field this could be either the particle
        mass (average should be False) or the density (average should
        be True).
    x: array, shape (N,)
        X coordinates of field samples, unit indices: [0,NX>.
    nx: int
        Number of grid points in X-direction.
    y: array, shape (N,), optional
        Y coordinates of field samples, unit indices: [0,NY>.
    ny: int, optional
        Number of grid points in Y-direction.
    wraparound: bool (False)
        If True, then values past the first or last grid point can
        wrap around and contribute to the grid point on the opposite
        side (see the Notes section below).

    Returns
    -------
    dens: ndarray, shape (nx, ny)
        The grid point values.

    Notes
    -----
    Example of default allocation of nearest grid points: nx = 4, * = gridpoint.

      0   1   2   3     Index of gridpoints
      *   *   *   *     Grid points
    |---|---|---|---|   Range allocated to gridpoints ([0.0,1.0> -> 0, etc.)
    0   1   2   3   4   posx

    Example of ngp allocation for wraparound=True: nx = 4, * = gridpoint.

      0   1   2   3        Index of gridpoints
      *   *   *   *        Grid points
    |---|---|---|---|--    Range allocated to gridpoints ([0.5,1.5> -> 1, etc.)
      0   1   2   3   4=0  posx


    References
    ----------
    R.W. Hockney and J.W. Eastwood, Computer Simulations Using Particles
        (New York: McGraw-Hill, 1981).

    Modification History
    --------------------
    IDL code written by Joop Schaye, Feb 1999.
    Avoid integer overflow for large dimensions P.Riley/W.Landsman Dec. 1999
    Translated to Python by Neil Crighton, July 2009.
    
    Examples
    --------
    >>> nx = 20
    >>> ny = 10
    >>> posx = np.random.rand(size=1000)
    >>> posy = np.random.rand(size=1000)
    >>> value = posx**2 + posy**2
    >>> field = cic(value, posx*nx, nx, posy*ny, ny)
    # plot surface
    """

    def findweights(pos, ngrid):
        """ Calculate CIC weights.
        
        Coordinates of nearest grid point (ngp) to each value. """

        if wraparound:
            # grid points at integer values
            ngp = np.fix(pos + 0.5)
        else:
            # grid points are at half-integer values, starting at 0.5,
            # ending at len(grid) - 0.5
            ngp = np.fix(pos) + 0.5

        # Distance from sample to ngp.
        distngp = ngp - pos

        # weight for higher (right, w2) and lower (left, w1) ngp
        weight2 = np.abs(distngp)
        weight1 = 1.0 - weight2

        # indices of the nearest grid points
        if wraparound:
            ind1 = ngp
        else:
            ind1 = ngp - 0.5
        ind1 = ind1.astype(int)

        # print 'ind',ind1,'max min ind',max(ind1),min(ind1)

        ind2 = ind1 - 1

        # Correct points where ngp < pos (ngp to the left).
        ind2[distngp < 0] += 2

        # Note that ind2 can be both -1 and ngrid at this point,
        # regardless of wraparound. This is because distngp can be
        # exactly zero.
        bad = (ind2 == -1)
        ind2[bad] = ngrid - 1
        if not wraparound:
            weight2[bad] = 0.
        bad = (ind2 == ngrid)
        ind2[bad] = 0
        bad1 = (ind1 == ngrid)
        ind1[bad1] = ngrid - 1

        if not wraparound:
            weight2[bad] = 0.

        if wraparound:
            ind1[ind1 == ngrid] = 0

        return dict(weight=weight1, ind=ind1), dict(weight=weight2, ind=ind2)

    def update_field_vals(field, value, a, b=None, debug=False):
        """ This updates the field array (and the totweight array if
        average is True).

        The elements to update and their values are inferred from
        a,b,c and value.
        """
        print 'Updating field vals'

        # weight per coordinate
        weights = a['weight'] * b['weight']
        # Don't modify the input value array, just rebind the name.
        value = weights * value
        indices = []

        for i in range(len(value)):
            field[a['ind'][i]][b['ind'][i]] += value[i]

        if debug: print i, weights[i], value[i], field[a['ind'][i]][b['ind'][i]]

    # pdb.set_trace()
    nx, ny = (int(i) for i in (nx, ny))
    nxny = nx * ny
    value = np.asarray(value)

    print 'Resampling %i values to a %i by %i grid' % (
        len(value), nx, ny)

    # normalise data such that grid points are at integer positions.
    # x = (x - x.min()) / x.ptp() * nx
    # if y!=None: y = (y - y.min()) / y.ptp() * ny
    # if z!=None: z = (z - z.min()) / z.ptp() * nz

    x1, x2 = findweights(np.asarray(x), nx)
    ind = []
    ind.append([x1, x2])
    if y is not None:
        y1, y2 = findweights(np.asarray(y), ny)
        ind.append([y1, y2])

    # float32 to save memory for big arrays (e.g. 256**3)
    field = np.zeros(shape=(nx, ny))  # field = np.zeros(shape=(nx,ny,nz)).squeeze()
    print 'shape field', field.shape

    update_field_vals(field, value, x1, y1)
    update_field_vals(field, value, x2, y1)
    update_field_vals(field, value, x1, y2)
    update_field_vals(field, value, x2, y2)

    return field
