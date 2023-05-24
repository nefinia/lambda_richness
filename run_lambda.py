#!/usr/bin/env python
from dataIO import load
from tools import rdarg2
from tools_sofi import degtohms, pdfw
import numpy as np
import os, sys
from math import *
from lambda_richness import lambda_richness, lambda_angdist
from background import background_cat, bkg_data
from pyfits import getdata

catalog = 'CFHTLS-full'

outname = 'output/lambdas-%s' % catalog
outname2 = outname + '-galaxies'

overwrite = True
if os.path.isfile(outname) and not overwrite: raise SystemExit('Output file already exists %s' % outname)
f1 = open(outname, 'w')
f2 = open(outname2, 'w')

cat = 'input/test_run.dat'
det = load(cat)
print 'Done loading detection'
cls = det.clusters

src = 'input/cfhtls-w1_magcuts.fits'
data = getdata(src, 1)
ra = np.array(data['alpha'])
dec = np.array(data['delta'])
zphot = np.array(data['photoz'])
zspec = np.array(data['zspec'])

G = np.array(data['G'])
R = np.array(data['R'])
I = np.array(data['I'])
Z = np.array(data['Z'])
eG = np.array(data['eG'])
eR = np.array(data['eR'])
eI = np.array(data['eI'])
eZ = np.array(data['eZ'])

lambda_list = []
zs = []
ras = []
decs = []
Ngals = []
extents = []

i = 0
f1.write(
    '#cluster_id cluster_dna name ra[deg] dec[deg] bcg_ra bcg_dec z_orca z_spec_lambda lambda lambda_error Ngals_orca Ngals_lambda Ngals_spec basis width slope intercept extent[arcmin] extent[Mpc]\n')
f2.write(
    '#cluster_id cluster_dna cluster_ra cluster_dec bcg_ra bcg_dec cluster_z_orca cluster_z_spec_lambda lambda ra dec zphot zspec radius prob weight G R I Z eG eR eI eZ\n')

#Area of the field. Very important for the background estimation!!
#This is the final unmasked area of CFHTLS W1 patch
#If swisscheese.py is applied the number should be modified!
area = 58.47
ncolorbins = 40
magzero = 14
colorbinsize = 0.1
nmagbins = 100
magbinsize = 0.1
clabels = ['gr', 'ri', 'iz']
index = {'gr': 0,  'ri': 1, 'iz': 2}
colors = [G - R, R - I, I - Z]
c2 = [R, I, Z]
color_errs = [np.sqrt(eG ** 2 + eR ** 2), np.sqrt(eR ** 2 + eI ** 2), np.sqrt(eI ** 2 + eZ ** 2)]
mags = [I, I, I]
colorzero = [-1, -1.5, -1.7]
bkg_full = []

for j in range(3):
    file = 'input/bkg_%s.dat' % clabels[j]
    print 'Loading background %s' % file
    if os.path.isfile(file):
        bkg = bkg_data(ncolorbins, colorbinsize, colorzero[j], nmagbins, magbinsize, magzero, np.loadtxt(file))
    else:
        bkg = background_cat(ncolorbins, colorbinsize, colorzero[j], nmagbins, magbinsize, magzero, colors[j],
                                   mags[j], area, smooth=True)
        np.savetxt(file, bkg.sigma_g)
    bkg_full.append(bkg)

for cl in cls:
    print '\nCluster', i

    z_orca = float(cl.z)
    z_cl = z_orca
    zspecl = -1
    extent = np.degrees(cl.extent)

    slope, width, intercept = cl.cmr(fit=True, colour=cl.basis)
    width = sqrt(width ** 2 + .05 ** 2)  # intrinsic color dispersion .05

    ra_bcg = np.degrees(cl.bcg.x)
    dec_bcg = np.degrees(cl.bcg.y)
    Ngals_orca = len(cl.galaxies)

    ang = lambda_angdist(z_cl)
    rad = ang * np.sqrt((np.radians(ra)-cl.bcg.x)**2+(np.radians(dec)-cl.bcg.y)**2)

    k = index[cl.basis]
    bkg = bkg_full[k]
    radfilter = np.where(rad <= cl.extent * ang)[0]

    color = colors[k][radfilter]
    color_err = color_errs[k][radfilter]
    d = color - (slope * (c2[k][radfilter] - 20.) + intercept)

    mag = I[radfilter]
    rad = rad[radfilter]

    d_err = np.sqrt(color_err*color_err+width*width)
    d_wt = (1. / (np.sqrt(2.*pi)*d_err))*np.exp(-(d*d)/(2.*d_err*d_err))
    lamb, pwt, rlambda, lambda_error, redseq = lambda_richness(z_cl, color, color_err, mag, rad, bkg, d_wt,
                                                               alpha=-1., _print=False)

    if len(d_wt) == 0:
        print 'No gals??'
        lamb = -1
        pwt = None

    if pwt == None:
        Ngals_used = 0
        zspecl = -1
        pwt = np.zeros(len(rad))
        pwt = [pwt, pwt]

    else:
        Ngals_used = len(np.where(pwt[1] > .3)[0])

        if len(zspec) > 0:
            zspecm = zspec[radfilter]

            # usar pdf!!!!!
            good = np.where((pwt[0] > .3) & (zspecm > 0.) & (zspecm < 1.5))[0]  # (zspecm>z_cl-.05)&(zspecm<z_cl+.05)&
            if len(good) > 1:
                # zspecl=np.median(zspecm[good])
                zspecl = pdfw(zspecm[good], pwt[0][good], p=1e-3)
                print 'Spectral redshift found!', zspecl, 'number of spectra used', len(good)
            else: zspecl = -1
        else:
            zspecl = -1
            good = []

    print 'cluster', i, 'lambda', lamb
    print 'total gals selected', len(color), 'gals with p>.3', Ngals_used, 'orca gals', Ngals_orca
    print 'orca redshift', z_orca, 'spec redshift', zspecl

    a = degtohms(np.degrees(cl.x), np.degrees(cl.y))
    if a[1][0] >= 0:
        sign = '+'
    else:
        sign = '-'
    # Code name for each cluster
    name = 'GMP_J%02d%02d%02d%s%02d%02d%02.1f' % (a[0][0], a[0][1], a[0][2], sign, abs(a[1][0]), a[1][1], a[1][2])

    f1.write('%d %d %s %f %f %f %f %f %f %f %f %d %d %d %s %f %f %f %f %f\n' % (
    i, cl.dna, name, np.degrees(cl.x), np.degrees(cl.y), ra_bcg, dec_bcg, z_orca, zspecl, lamb, lambda_error,
    Ngals_orca, Ngals_used, len(good), cl.basis, width, slope, intercept, extent * 60., ang * cl.extent))

    ss0 = '%d %d %f %f %f %f %f %f %f' % (
    i, cl.dna, np.degrees(cl.x), np.degrees(cl.y), ra_bcg, dec_bcg, z_orca, zspecl, lamb)
    ss = ''
    members = np.where(pwt[0] > .3)[0]

    if len(zspec) > 0:
        zss = zspec[radfilter]
    else:
        zss = np.zeros(len(radfilter))

    for mm in members:
        ss += '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n' % (
        ss0, ra[radfilter][mm], dec[radfilter][mm], zphot[radfilter][mm], zss[mm], rad[mm], \
        pwt[0][mm], pwt[1][mm], G[radfilter][mm], R[radfilter][mm], I[radfilter][mm], Z[radfilter][mm],
        eG[radfilter][mm], eR[radfilter][mm], eI[radfilter][mm], eZ[radfilter][mm])

    if len(members) > 0:
        f2.write(ss)
    else:
        print 'No bcg?? no galaxies with p>.3???'

    ss = ''
    ss0 += ' %d' % len(good)
    i += 1

f1.close()
f2.close()