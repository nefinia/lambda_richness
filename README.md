# lambda_richness

Lambda richness estimator

minimal inputs:
- redshift of the cluster
- color (initially g-r)
- color error
- magnitude (initially i-band but in principle can change)
- distance from the center of the cluster for each galaxy
- background structure

minimal output:
-  the lambda optical richness of the cluster


Main diferences from the Rykoff+12 lambda richness code:

- Translated to python (instead of idl)
- Update mstar estimation (updated to Rykoff+14)
- Red sequence parameters are fitted for each cluster (slope, width, intercept). Therefore the color distribution (d_wt in the code) 
  will depend on the current cluster, instead of a model which is a function of redshift (see eq. 32 of Rykoff+14)
- Background is estimated from the galaxy catalog (in the webpage there was a matrix with values valid only for SDSS). 
  The algorithm applies a cloud-in-cell technique for smoothing (background.py).
- Extra cluster redshift estimator based on the probabilities of the galaxies to be part of the cluster (pwt in the code) and the spectral redshifts if available.
- OPTIONAL: swisscheese.py will upgrade the background estimation by removing all galaxies inside a cluster (haven't really tested lately)

For a test on the code execute run_lambda.py

The inputs (in folder input) are:
- orca detection file: test_run.dat
- galaxy catalog: cfhtls-w1_magcuts.fits
- background matrix for each color set: bkg_gr.dat, bkg_ri.dat, bkg_iz.dat

The outputs are:
- A catalog of clusters with their lambda value and other properties
- A catalog of galaxies with probability of being part of each cluster higher than 0.3
