# M4 (NGC 6121) Proper Motion Membership Catalog

This repository contains the M4 membership catalog produced by () from Gaia DR2 proper 
motion measurements, as well as the code used to create the catalog.  The contents of this repo 
are distributed under the BSD 3-Clause License (see LICENSE file for more information).

## Dependencies

The following Python packages are used by the main code (pm_membership.py):
[astropy](http://www.astropy.org/), [matplotlib](https://matplotlib.org/), 
[numpy](http://www.numpy.org/),
[pickle](https://docs.python.org/3/library/pickle.html), 
[scikit-learn](http://scikit-learn.org/stable/)

## Repo Contents

* pm_membership.py -- the code used to generate the membership catalog, as well as produce 
 a plot showing the fit results. The code takes two arguments to run: 1. the name of the 
input VOTable file following the Gaia data format, and 2. the name of the output pickle file
that stores the fitted Gaussian mixture model.

* gaia_query.py -- this code can be used to query the Gaia data archive for the appropriate 
data.  The call in the code currently reproduces the Gaia data I used to make the fitted model.

* gmm_fitted.pkl -- the fitted model.  This is a Python pickle of an instance of an 
[sklearn.mixture.GaussianMixture](http://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html)
object.  This pickle was generated using Python 2, but can be read in to Python 3 with the 
`encoding='latin'` option to your `pickle.load()` call.

* catalog.txt -- the catalog in plain text.  See header for column information.

## How to Use



## Acknowledgment and Attribution

The official publication for this catalog is (fill in when published).  
This repo also has a DOI, ().
Any use of this catalog or the associated code should cite the Wallace et al. article, 
and any use of the data or code in the repo should also cite the repo DOI.