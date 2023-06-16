# Lambda Richness Estimator

The Lambda Richness Estimator is a Python implementation of a code for estimating the optical richness of galaxy clusters based on various input parameters. It provides a measure of the cluster's richness, which is a proxy for its mass and can be used to study galaxy cluster populations.

## Minimal Inputs

To run the Lambda Richness Estimator, the following minimal inputs are required:

- Redshift of the cluster
- Color (initially g-r)
- Color error
- Magnitude (initially i-band, but can be changed)
- Distance from the center of the cluster for each galaxy
- Background structure information

## Minimal Output

The Lambda Richness Estimator provides the following minimal output:

- The lambda optical richness of the cluster

## Main Differences from the Rykoff+12 Lambda Richness Code

Compared to the original Rykoff+12 lambda richness code, the Lambda Richness Estimator offers several improvements:

1. Translated to Python from IDL for easier use and compatibility.
2. Updated estimation of the stellar mass (mstar) based on Rykoff+14.
3. Red sequence parameters (slope, width, intercept) are fitted for each cluster individually, resulting in a color distribution (d_wt in the code) specific to the current cluster, rather than a fixed model dependent on redshift (see eq. 32 of Rykoff+14).
4. Background estimation is performed using the galaxy catalog itself, rather than relying on a pre-defined matrix specific to SDSS. The algorithm applies a cloud-in-cell technique for smoothing (background.py).
5. An additional cluster redshift estimator is provided based on the probabilities of galaxies being part of the cluster (pwt in the code) and spectral redshifts if available.
6. OPTIONAL: The swisscheese.py module can be used to improve the background estimation by removing all galaxies inside a cluster (please note that this feature hasn't been extensively tested).

## Getting Started

To run a test on the code, execute the run_lambda.py script. The input files required are located in the "input" folder, including:

- Orca detection file: test_run.dat
- Galaxy catalog: cfhtls-w1_magcuts.fits
- Background matrix for each color set: bkg_gr.dat, bkg_ri.dat, bkg_iz.dat

The outputs include:

- A catalog of clusters with their lambda value and other properties
- A catalog of galaxies with a probability of being part of each cluster higher than 0.3

Feel free to explore and customize the code to suit your specific needs and datasets.

## Contributors

This Lambda Richness Estimator code was developed by Sofia G. Gallego and is based on the work of Rykoff et al. (2012) and Rykoff et al. (2014).
