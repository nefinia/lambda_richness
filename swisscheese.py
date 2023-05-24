#!/usr/bin/env python
preamble = """
####################################################
#
# Catalog with no clustered galaxies generator
# Autor: Sofia C. Gallego, September 2014 
# 
####################################################
"""

import numpy as np
import pyfits as pf
from math import *
import os


def set_clusters(det_type='orca', det_file='', folderout=''):
    "Use detected clusters as holes for the cheese"

    ra = []
    dec = []
    z = []
    extent = []
    N_gal = []

    if det_type == 'orca':
        print 'ORCA input'
        N_assoc = []

        f1 = '%sclusters_info.dat' % folderout

        if os.path.isfile(f1):
            'Info file already exists', f1
            mat = np.loadtxt(f1)
        else:
            from dataIO import load
            det = load(det_file)
            print 'Done loading detection file'
            cls = det.clusters

            for cl in cls:
                ra.append(np.degrees(cl.x))
                dec.append(np.degrees(cl.y))
                z.append(cl.z)
                extent.append(np.degrees(cl.extent))
                N_gal.append(cl.ngalaxies)

            mat = np.array([ra, dec, z, extent, N_gal])
            f1 = '%sclusters_info.dat' % folderout
            print f1, mat.T
            np.savetxt(f1, mat.T, fmt=('%f', '%f', '%f', '%f', '%d'))

    else:
        print 'Please specify an existing detection type'
        return

    return mat


def get_cheese(ra_c, dec_c, extent, N_gal, ra, dec, g, r, i, zz, eG, eR, eI, eZ, clid=-1, size=60, cat_type='SDSS',
               folderout='', full=False):
    # ra, dec and size should be in degrees, size in arcmin

    if full:
        ra_ins = ra
        dec_ins = dec
        incls = range(len(ra_c))

    else:
        if clid == -1:
            print 'Please provide cluster id'
            raise SystemExit(0)

        print 'Number of galaxies in cluster%d =' % clid, N_gal[clid]

        in_cluster = np.where((ra - ra_c[clid]) ** 2 + (dec - dec_c[clid]) ** 2 <= extent[clid] ** 2)[0]
        mat = np.array([ra[in_cluster], dec[in_cluster], g[in_cluster], r[in_cluster], i[in_cluster], zz[in_cluster],
                        eG[in_cluster], eR[in_cluster], eI[in_cluster], eZ[in_cluster]])
        np.savetxt('%s/cluster%d_all-gals-inside' % (folderout, clid), mat.T,
                   fmt=('%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f'))

        ramax = ra_c[clid] + size / 120.
        ramin = ra_c[clid] - size / 120.
        decmax = dec_c[clid] + size / 120.
        decmin = dec_c[clid] - size / 120.

        incls = []
        incls.append(clid)

        others = np.where(np.array(range(len(ra_c))) != clid)[0]
        print others

        for j in others:
            if is_cluster_inside(ra_c[j], dec_c[j], extent[j], ramin, ramax, decmin, decmax):
                incls.append(j)

        print 'Number of clusters inside area', len(incls)

        print 'Catalog limits', min(ra), max(ra), min(dec), max(dec)
        print 'Cheese limits', ramin, ramax, decmin, decmax

        inside = np.where((ra < ramax) & (dec < decmax) & (ra > ramin) & (dec > decmin))[0]

        ra_ins = ra[inside]
        dec_ins = dec[inside]

    ra_inc = ra_ins
    dec_inc = dec_ins
    ra_exc = []
    dec_exc = []

    for k in incls:
        inc = np.where((ra_inc - ra_c[k]) ** 2 + (dec_inc - dec_c[k]) ** 2 > extent[k] ** 2)[0]
        ra_inc = ra_inc[inc]
        dec_inc = dec_inc[inc]

        if not full:
            exc = np.where((ra_ins - ra_c[k]) ** 2 + (dec_ins - dec_c[k]) ** 2 <= extent[k] ** 2)[0]
            for j in exc:
                ra_exc.append(ra_ins[j])
                dec_exc.append(dec_ins[j])

    inc = np.where(ra_inc)[0]
    print len(inc)

    mat = np.array([ra_inc, dec_inc, g[inc], r[inc], i[inc], zz[inc], eG[inc], eR[inc], eI[inc], eZ[inc]])

    if full:
        np.savetxt('%s/%s-full_cheese' % (folderout, cat_type), mat.T,
                   fmt=('%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f'))
        print 'Number of background galaxies outside clusters in the full area', len(
            inc)  # ,'number of galaxies removed',len(exc)

    else:
        np.savetxt('%s/cluster%d_cheese' % (folderout, clid), mat.T,
                   fmt=('%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f'))
        print 'Number of galaxies outside clusters in area', len(inside), 'number of galaxies removed', len(ra_exc)

    return mat


def is_cluster_inside(ra, dec, rad, ramin, ramax, decmin, decmax):
    """A cluster is inside the area when the intersection of the cluster boundary (assuming a circle)
       with the area boundaries exists"""

    left = is_intersection(1, -2 * dec, dec ** 2 + (ramin - ra) ** 2 - rad ** 2, decmin, decmax)
    right = is_intersection(1, -2 * dec, dec ** 2 + (ramax - ra) ** 2 - rad ** 2, decmin, decmax)
    top = is_intersection(1, -2 * ra, ra ** 2 + (decmin - dec) ** 2 - rad ** 2, ramin, ramax)
    bottom = is_intersection(1, -2 * ra, ra ** 2 + (decmax - dec) ** 2 - rad ** 2, ramin, ramax)
    center = False
    if ra >= ramin and ra <= ramax and dec >= decmin and dec <= decmax: center = True

    if (left or right or top or bottom or center):
        return True

    else:
        return False


def is_intersection(a, b, c, minlim, maxlim):
    "Solving quadratic equation constricted to image boundaries"

    d = b ** 2 - 4 * a * c

    if d < 0:
        return False

    else:
        sqrd = sqrt(d)
        sol1 = (-b + sqrd) / (2 * a)
        sol2 = (-b + sqrd) / (2 * a)

        if ((minlim < sol1 < maxlim) | (minlim < sol2 < maxlim)):
            return True
        else:
            return False


if 1:

    # cat_type='CFHTLS'#'CFHTLS','SDSS','Conito',etc...
    cat_type = 'CFHTLS'  # -sextractor'
    full = True
    det_type = 'orca'

    if cat_type == 'CFHTLS':
        fits = 'input/cfhtls-w1_magcuts.fits'
        folderout = 'input/'
        det_file = 'input/test_run.dat'
        data = pf.getdata(fits, 1)
        ra = data['alpha']
        dec = data['delta']
        g = data['G']
        r = data['R']
        i = data['I']
        z = data['Z']
        eG = data['eG']
        eR = data['eR']
        eI = data['eI']
        eZ = data['eZ']

    mat = set_clusters(det_type, det_file, folderout)
    print 'Done setting clusters'

    mat = np.array(mat)
    ra_c = mat[0]
    dec_c = mat[1]
    # z=mat[2]
    extent = mat[3]
    N_gal = mat[4]

    clid = 0

    size = 60

    if full:
        data_cheese = get_cheese(ra_c, dec_c, extent, N_gal, ra, dec, g, r, i, z, eG, eR, eI, eZ, cat_type=cat_type,
                                 folderout=folderout, full=full)
        print 'End'


    else:
        for clid in range(len(mat[0])):
            print clid
            data_cheese = get_cheese(ra_c, dec_c, extent, N_gal, ra, dec, g, r, i, z, eG, eR, eI, eZ, clid, size=size,
                                     cat_type=cat_type, folderout=folderout)
            # tbhdu=pf.new_table(data_cheese) #problem!!!! makes a sdss fits file :S
            # tbhdu.writeto('/data3/scgalleg/orca/background/cfhtls_w1-cluster%d-bkg_gals.fits'%clid)
            print 'Done getting swiss cheese for cluster%d' % clid  # , len(data_cheese)

        print 'End'
