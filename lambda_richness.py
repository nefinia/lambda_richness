#!/usr/bin/env python
notes = """\
####################################################
#
# Python adaptation of the IDL code by Eli S. Rykoff
# See Rykoff et al 2012 for references
# Autor: Sofia C. Gallego, August 5 2014 
# 
####################################################

 
###############################################################
#
# NAME:
#   LAMBDA_RICHNESS
# PURPOSE:
#   Optimized red-sequence redshift estimator for galaxy clusters
# EXPLANATION:
#   Calculate lambda richness from Rykoff et al. (2011).  Lambda is an
#   optimized richness estimator, designed to minimize the scatter in
#   LX at fixed richness with extensive tests on the maxBCG cluster
#   catalog (Koester et al. 2007).  The probabilistic methodology is
#   very robust to various perturbations, including changes in the
#   color used to calculate the red sequence, background shifts, and
#   zero-point uncertainties.
#   
#   The specific implementation here is appropriate for SDSS data and
#   galaxy clusters at 0.1<z<0.35.  If modifying, please be sure to
#   update the red sequence model, estimation of m*, and background structure.
#
#   The aperture used to measure each cluster is optimized based on the
#   richness of the cluster, such that the cutoff radius is given by:
#     rc = r0*(lambda/100)^beta
#
# CALLING SEQUENCE:
#   lambda = lambda_richness(z,color,color_err,imag,r,r0,beta,lstarcut,h0,bkg,
#                            rlambda,lambda_err,wtvals,pvals)
#
# INPUTS:
#   z - the redshift of the cluster
#   color - the g-r color (dereddened MODEL_MAG) for each galaxy within the
#          region of interest (at least 1.5 h^-1 Mpc from the center)
#   color_err - the g-r color error for each galaxy
#   imag - the i-band magnitude (deredded CMODEL_MAG) for each galaxy
#   r - the radius from the center (h^-1 Mpc) for each galaxy
#
# OPTIONAL INPUTS:
#   r0 - the constant in the radial scaling (h^-1 Mpc) [default 1.0]
#   beta - the exponent in the radial scaling [default 0.2]
#   lstarcut - the fraction of L* to use as a cut [default 0.2]
#   h0 - value for H0 used to define R0, r [default 100.0]
#   bkg - the background structure, consisting of the following fields:
#         .color - the color color bins [findgen(40)*colorbinsize-1.0]
#         .imag - the imag magnitude bins [findgen(100)*imagbinsize+12.0]
#         .colorbinsize - the size of the color bins [0.1 mag]
#         .imagbinsize - the size of the mag bins [0.1 mag]
#         .sigma_g - background density in N/deg.sq/mag/mag
#
# RETURN VALUE:
#   lambda - the lambda richness for the cluster
#
# OPTIONAL OUTPUTS:
#   bkg - background matrix in color-magnitude space 
#   rlambda - the cutoff radius used: rlambda = r0*(lambda/100)^beta
#   lambda_err - error in lambda estimated from Eqn. 3 in Rykoff et al. 2011
#   wtvals - weights for individual galaxies brighter than the luminosity
#             cutoff and within rlambda.  lambda = total(wtvals)
#             identically 0 for all galaxies outside the luminosity/radial
#             range
#   pvals - probabilities for all galaxies input.  For galaxies brighter than
#             the luminosity cutoff and within rlambda wtvals=pvals.
#             Outside the range, we can still calculate a probability for the
#             galaxy, with the caveat that the sum of probabilities will not
#             be properly normalized.
#
# CALLED ROUTINES:
#   LAMBDA_NFW_WEIGHT, LAMBDA_LUM_FUNC, LAMBDA_ANGDIST, LAMBDA_BCOUNTS,
#    LAMBDA_WEIGHTS, LAMBDA_BISECTION_ZERO, LAMBDA_SDSS_BACKGROUND
#
# REFERENCES:
#   Rykoff, E. S. et al 2012
# 
# REVISION HISTORY
#   Author: Eli S. Rykoff, LBNL, March 2011
#-
######################################################################
"""

import numpy as np
from math import *
from scipy import integrate


##############################################################
# MAIN ROUTINE
##############################################################


def lambda_richness(z, color, color_err, mag, r, bkg=None, d_wt=None, wtvals=None, pvals=None, r0=1., beta=.2,
                    lstarcut=.2, h0=100., alpha=-1., _print=False, R08=False, R12=False, R14=False, optoutput=True):
    "MAIN ROUTINE"

    ## z is redshift
    ## color, color_err are MODEL colors and errors
    ## mag is CMODEL magnitude (typically i magnitude)
    ## r is radius from center in h^-1 Mpc

    if len(color) == 0:
        print 'Empty data!'
        return -1, None

    if _print:
        print 'Bkg maximum', np.amax(bkg.sigma_g)

    if bkg == None:
        print 'Please specify background'
        raise SystemExit(0)

    color = np.array(color)
    color_err = np.array(color_err)
    mag = np.array(mag)
    r = np.array(r)

    if type(z) not in [type(0.5), type(1)]:
        print 'ERROR: can only specify one redshift z.'
        return -1, None

    if len(color) != len(color_err) or len(color) != len(mag) or len(color) != len(r):
        print 'ERROR: number of elements of color, color_err, mag, r must be identical.'
        return -1, None

    ## First, the luminosity function

    ## we need mstar...

    ## Rozo+ 2008 mstar, this is valid for 0.05<z<0.35
    if R08 or R12:
        mstar = 12.27 + 62.36 * z - 289.79 * z * z + 729.69 * pow(z, 3.) - 709.42 * pow(z, 4.)  # eq 12 Rykoff+ 2012

    # Rykoff+ 2014 mstar, valid from 0.05<z<0.7, appropiate for SDSS DR8
    elif R14:
        if z <= .5:
            mstar = 22.44 + 3.36 * log(z) + .273 * log(z) ** 2 - .0618 * log(z) ** 3 - .0227 * log(z) ** 4
        else:
            mstar = 22.94 + 3.08 * log(z) - 11.22 * log(z) ** 2 - 27.11 * log(z) ** 3 - 18.02 * log(z) ** 4

    else:
        if z <= .5:
            mstar = 22.44 + 3.36 * log(z) + .273 * log(z) ** 2 - .0618 * log(z) ** 3 - .0227 * log(z) ** 4
        else:
            mstar = 22.94 + 3.08 * log(z) - 11.22 * log(z) ** 2 - 27.11 * log(z) ** 3 - 18.02 * log(z) ** 4

    ## Calculate the limiting mag

    limmag = mstar - 2.5 * log10(lstarcut)

    if _print:
        print 'mstar', mstar, 'limmag', limmag

    lum_wt = lambda_lum_func(mag, mstar, limmag, alpha=alpha)

    if _print:
        print 'lum_wt', lum_wt


    sigma_int = .05

    ## Now calculate the color filter G(c)
    bright = np.where(mag <= limmag)[0]

    if R08:
        # Color filter as in Rozo+ 2008
        d = color - 0.625 + 3.149 * z

    elif R12 or R14 or len(bright) == 0:
        # Color filter as in Rykoff+ 2012
        ## define the red sequence tilt and intercept as a function of redshift (valid just for g-r color!)
        m = -0.008969 - 0.0701 * z #m17
        b = 0.58907 + 3.2982 * z #b17

        ## calculate G(c) filter
        d = color - (m * (mag - 17.0) + b)

        redseqprop = [m, b]

    else:
        # Fit values with the red sequence of each cluster
        m, b = np.polyfit(mag[bright], color[bright], 1)
        print 'polifyt', m, b
        d = color - m * mag - b

    if d_wt is None:
        d_err = np.sqrt(color_err * color_err + sigma_int * sigma_int)
        d_wt = (1. / (np.sqrt(2. * pi) * d_err)) * np.exp(-(d * d) / (2. * d_err * d_err))

    if _print:
        print 'd_wt', d_wt

    ## Third, calculate the NFW filter

    nfw_wt = lambda_nfw_weight(r)

    if _print:
        print 'nfw_wt', nfw_wt

    ## We can now calculate the filter u:
    ucounts = d_wt * nfw_wt * lum_wt

    if _print:
        print 'ucounts', ucounts

    ## make sure we don't have any infinities
    ucounts[np.isinf(ucounts)] = 0.

    ## and now we need the background
    bcounts = lambda_bcounts(bkg, z, r, color, mag)

    if _print:
        print 'bcounts', bcounts

    ## finally, the weight function for the luminosity cutoff
    ## 1.0 if it's above the cut, 0.0 below
    # theta_i[np.where(mag < limmag)]=1.
    theta_i = np.zeros(len(mag))
    br = np.where(mag < limmag)[0]

    if len(br) > 0:
        theta_i[br] = 1.

    if _print:
        print 'theta_i', theta_i

    lamb, pwt = lambda_bisection_zero(r0, beta, ucounts, bcounts, r, theta_i)
    lambda_err = -1

    if wtvals is not None:
        wtvals = pwt[1]

    if lamb > 0:
        rlambda = r0 * pow(lamb / 100., beta)

        ## and calculate the error on lambda
        if wtvals is not None:
            bar_p = sum(wtvals * wtvals) / sum(wtvals)
            lambda_err = sqrt(lamb * (1. - bar_p))
            print 'lambda error', lambda_err

    else:
        print 'Que paso?! lambda error!!'
        return -1, None, 0, -1, [0, 0]

    if not R08 and not R12 and not R14:
        redseqprop = [b, m - 20 * b]
    else:
        redseqprop = [0, 0]

    if optoutput: return lamb, pwt, rlambda, lambda_err, redseqprop
    else: return lamb


def lambda_nfw_weight(r, rscale=0.15):
    "Routine to compute the NFW filter weight (see Eqns 6-8 in Rykoff et al. 2012)"
    # from operator import countOf

    x = r / rscale

    nx = len(x)
    sigx = np.zeros(len(x))

    ## core is 100 kpc
    corer = 0.1
    corex = corer / rscale

    low = np.where(x < corex)[0]
    nlow = len(low)
    mid = np.where((x >= corex) & (x < 1.0))[0]
    nmid = len(mid)
    high = np.where((x > 1.0) & (x < 10.0 / rscale))[0]
    nhigh = len(high)
    other = np.where((x > 0.999) & (x < 1.001))[0]
    nother = len(other)

    if nlow > 0:
        arg = np.sqrt((1. - corex) / (1. + corex))
        pre = 2. / np.sqrt(1. - corex ** 2.)
        front = 1. / (corex ** 2. - 1)
        sigx[low] = front * (1. - pre * 0.5 * log((1. + arg) / (1. - arg)))

    if nmid > 0:
        arg = np.sqrt((1. - x[mid]) / (1. + x[mid]))
        pre = 2. / (np.sqrt(1. - x[mid] ** 2))
        front = 1. / (x[mid] ** 2 - 1)
        sigx[mid] = front * (1. - pre * 0.5 * np.log((1 + arg) / (1. - arg)))

    if nhigh > 0:
        arg = np.sqrt((x[high] - 1.) / (x[high] + 1.))
        pre = 2. / np.sqrt(x[high] ** 2 - 1.)
        front = 1. / (x[high] ** 2 - 1.)
        sigx[high] = front * (1. - pre * np.arctan(arg))

    if nother > 0:
        xlo = 0.999
        xhi = 1.001
        arglo = np.sqrt((1. - xlo) / (1. + xlo))
        prelo = 2. / np.sqrt(1. - xlo ** 2)
        frontlo = 1. / (xlo ** 2 - 1)
        testlo = frontlo * (1. - prelo * 0.5 * log((1 + arglo) / (1. - arglo)))
        arghi = np.sqrt((xhi - 1.) / (xhi + 1.))
        prehi = 2. / (sqrt(xhi ** 2 - 1.))
        fronthi = 1. / (xhi ** 2 - 1.)
        testhi = fronthi * (1. - prehi * np.arctan(arghi))
        sigx[other] = (testlo + testhi) / 2.0

    wt = 2. * pi * r * sigx

    ## we need to make sure that things aren't blowing up at 0:
    smallr = np.where(r < 1e-6)[0]
    if len(smallr) > 0:
        wt[smallr] = 2. * pi * (1e-6) * sigx[smallr]

    # print 'wt nfw',wt

    return wt


def lambda_lum_func(mag, mstar, limmag, alpha=-1.):
    "Routine to compute luminosity filter weight (see Eqn 10 in Rykoff et al. 2012)"

    ## integrate
    n = integrate.quad(schechter, 10, limmag, args=(alpha, mstar))

    wts = schechter(mag, alpha, mstar) / n[0]

    return wts


def schechter(m, alpha, mstar):
    "Schechter luminosity function"

    s1 = 10. ** (-0.4 * (m - mstar) * (alpha + 1))
    s2 = -10. ** (-0.4 * (m - mstar))

    return s1 * np.exp(s2)


def lambda_oneovere(z, omega_l=0.7):
    "Routine to compute 1/E(z) for angular diameter distance calculation"

    omega_m = 1. - omega_l

    return 1. / sqrt(omega_m * pow(1 + z, 3.) + omega_l)


def lambda_angdist(z, h0=100, omega_l=0.7):
    "Routine to compute angular diameter distance"

    c = 2.99792458e5
    dh = c / h0

    dm = integrate.quad(lambda_oneovere, 0, z, args=omega_l)

    return dh * dm[0] / (1 + z)


def lambda_bcounts(bkg, z, r, color, mag, h0=100):
    "Routine to compute b(x) = 2*pi*r*Sigma_g(m,c) (see S3.4 in Rykoff et al. 2011)"

    ncolor = bkg.ncolorbins
    nmag = bkg.nmagbins

    colorindex = np.array((color - bkg.color[0]) * ncolor / ((bkg.color[-1] + bkg.colorbinsize) - bkg.color[0]),
                          dtype='int')
    magindex = np.array((mag - bkg.mag[0]) * nmag / ((bkg.mag[-1] + bkg.magbinsize) - bkg.mag[0]), dtype='int')
    colorpos = colorindex * bkg.colorbinsize + bkg.colorzero
    magpos = magindex * bkg.magbinsize + bkg.magzero

    ## Check for overruns
    badcolor = np.where((colorindex >= ncolor) | (colorindex < 0))[0]
    badmag = np.where((magindex >= nmag) | (magindex < 0))[0]

    colorindex[badcolor] = 0
    magindex[badmag] = 0

    sigma_g = np.array([bkg.sigma_g[i][j] for i, j in zip(magindex, colorindex)])

    ## and set all bad elememts to infinity: these galaxies should not be given any
    ## weight

    sigma_g[badcolor] = float('inf')
    sigma_g[badmag] = float('inf')

    # I'm not using inside cluster's galaxies so the bin could be zero --> no need to change value
    # badcombination = np.where(sigma_g == 0.0)[0]
    # if len(badcombination)>0:
    #    print 'Bkg bin zero in %d galaxies!'%len(badcombination)
    # sigma_g[badcombination] = float('inf')


    ## get the conversion factor to go from deg. to Mpc
    angdist = lambda_angdist(z, h0)

    c = (angdist * 1.0) * (pi / 180.)

    ## calculate b_i
    bcounts = 2. * pi * r * (sigma_g / c / c)

    return bcounts


def lambda_weights(inlambda, r0, beta, ucounts, bcounts, r, theta_i):
    "Routine to calculate p_i = lambda*u/(lambda*u+b) (see Eqn. 1 in Rykoff et al. 2011)"

    ## figure out the radial cut
    rc = r0 * pow(inlambda / 100., beta)

    ## calculate the normalization for the NFW filter for this radial cut
    ##  (see Eqns. 8-9 in Rykoff et al. 2011)
    nfwnorm = exp(
        1.65169 - 0.547850 * log(rc) + 0.138202 * log(rc) * log(rc) - 0.0719021 * pow(log(rc), 3.) - 0.0158241 * pow(
            log(rc), 4) - 0.000854985 * pow(log(rc), 5.))

    ## we only want to weight those that are inside...
    ## our weight function is theta
    theta_r = np.zeros(len(r))
    ins = np.where(r <= rc)[0]
    nins = len(ins)
    if (nins > 0):
        theta_r[ins] = 1.0  ## 0.0 otherwise

    ## store both "p" the probabilities for all galaxies in the region, and "wt"
    ## the weights used to calculate lambda
    nuc = np.zeros(len(ucounts))[0]
    pwt = [nuc, nuc]

    pwt[0] = (inlambda * ucounts * nfwnorm) / (inlambda * ucounts * nfwnorm + bcounts)

    ## fix any NaNs (if both ucounts and bcounts are 0.0)

    pwt[1] = pwt[0] * theta_r * theta_i

    return pwt


def lambda_bisection_zero(r0, beta, ucounts, bcounts, r, theta_i, tol=1e-3, wtvals=None, pvals=None):
    "Routine to calculate the zero in the function lambda - Sum(p)=0 (see Eqn. 2 in Rykoff et al. 2011)"

    ## simple bisector to find the zero of the function

    lamhi = 800.0
    lamlo = 0.5
    outlo = -1
    while abs(lamhi - lamlo) > 2 * tol:
        mid = (lamhi + lamlo) / 2.
        if outlo < 0:
            pwt = lambda_weights(lamlo, r0, beta, ucounts, bcounts, r, theta_i)
            # print 'type(pwt[1])',type(pwt[1])
            outlo = sum(pwt[1])
        pwt = lambda_weights(mid, r0, beta, ucounts, bcounts, r, theta_i)
        # print pwt
        outmid = sum(pwt[1])

        if outlo < 1.0:
            outlo = 0.9  ## provides stability at low
        if (outlo - lamlo) * (outmid - mid) > 0:
            ## we need the upper bracket
            lamlo = mid
            outlo = -1
        else:
            lamhi = mid  ## we need the lower bracket

    lamb = (lamlo + lamhi) / 2.0
    pwt = lambda_weights(lamb, r0, beta, ucounts, bcounts, r, theta_i)
    pvals = pwt[0]
    wtvals = pwt[1]

    if lamb < 1.0:
        ## failure to find any galaxies
        lamb = -1

    return lamb, pwt
