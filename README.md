# M4 (NGC 6121) Proper Motion Membership Catalog

This repository contains the M4 membership catalog produced by () from Gaia DR2 proper 
motion measurements, as well as the code used to create the catalog.  The contents of this repo 
are distributed under the BSD 3-Clause License (see [LICENSE](LICENSE.md) 
file for more information).

## Dependencies

The following Python packages are used by the main code (pm_membership.py):
[astropy](http://www.astropy.org/), [matplotlib](https://matplotlib.org/), 
[numpy](http://www.numpy.org/),
[pickle](https://docs.python.org/3/library/pickle.html), 
[scikit-learn](http://scikit-learn.org/stable/).

Additionally, the code to query the Gaia archive (gaia_query.py)
uses astropy and [astroquery](https://astroquery.readthedocs.io/en/latest/).

## Repo Contents

* [pm_membership.py](pm_membership.py) -- the code used to generate the 
membership catalog, as well as produce a plot showing the fit results. The 
code takes two arguments to run: 1. the name of the input VOTable file 
following the Gaia data format, and 2. the name of the output pickle file
that stores the fitted Gaussian mixture model.

* [gaia_query.py](gaia_query.py) -- this code can be used to query the Gaia 
data archive for the appropriate data.  The call in the code currently 
reproduces the Gaia data I used to make the fitted model.

* gmm_fitted*.pkl -- the fitted models.  These are Python pickles of instances 
of an 
[sklearn.mixture.GaussianMixture](http://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html)
object.  The pickles were generated using Python 2, but can be read in to 
Python 3 with the `encoding='latin'` option to your `pickle.load()` call.
The first component of the output to the `predict_proba()` method is the
cluster membership probability (the second component is the field 
membership probability).

    * [gmm_fitted_maglessthan_19.0.pkl](gmm_fitted_maglessthan_19.0.pkl) -- the 
fitted model for objects with Gaia *G* magnitudes less (brighter) than 19.  

    * [gmm_fitted_maggreaterthan_19.0.pkl](gmm_fitted_maggreaterthan_19.0.pkl) -- 
the fitted model for objects with Gaia *G* magnitudes greater (dimmer) than 19

* [membership_catalog.txt](membership_catalog.txt) -- the catalog in plain 
text.  Contains Gaia DR2 ID numbers, our calculated membership probabilities, 
and some other information (*G* magnitudes, RA, dec, proper motion values and 
errors on the latter three). See header for column and unit information.

## How to Use

**If all you want is the catalog published as-is,**
[membership_catalog.txt](membership_catalog.txt) has what you want in plain 
text, and you don't need to worry about anything else in this repo.

**If all you want is the fitted model,** then 
[gmm_fitted_maglessthan_19.0.pkl](gmm_fitted_maglessthan_19.0.pkl) and 
[gmm_fitted_maggreaterthan_19.0.pkl](gmm_fitted_maggreaterthan_19.0.pkl)
are what you want.  These are the fitted `sklearn.mixture.GaussianMixture` 
objects from `scikit-learn`.  Read them in as you would any pickle (they
were generated using Python 2, so reading into Python 3 will require the
`encoding='latin'` option.

**If you want to understand how the catalog was generated, or use
the same code to generate your own modified catalog,** you will want to look
at [gaia_query.py](gaia_query.py), which shows the ADQL query to the Gaia
data archive for the input data for the membership calculation, and 
[pm_membership.py](pm_membership.py), which has code to fit the model, 
generate the membership catalog from the fitted model, and plot up
a visualization of the fitted model.

The code can be ran as-is as following to regenerate the membership catalog
as published:

```
python gaia_query.py foo.xml # Downloads a portion of the Gaia data archive
python pm_membership.py foo.xml gmm_fitted 19
```

The first line gets the relevant data from the Gaia data archive and saves
it to the file `foo.xml`.  The second line then reads the data from `foo.xml`,
fits and saves the models to pickle files prepended `gmm_fitted`, and with
the model fits split between objects brighter than and dimmer than *G*=19.

See [pm_membership.py](pm_membership.py) for more detailed information on the
code itself.

## Caveats

Initialization

Sample size effects

magnitude effect


## Use, Acknowledgment, and Attribution

The contents of this repo are distributed under the BSD 3-Clause License 
(see [LICENSE](LICENSE.md) file for more information).

The official publication for this catalog is (fill in when published).  
This repo also has a DOI, ().
Any use of this catalog or the associated code should cite the Wallace et al. article, 
and any use of the data or code in the repo should also cite the repo DOI.