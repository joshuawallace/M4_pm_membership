[![DOI](https://zenodo.org/badge/154380587.svg)](https://zenodo.org/badge/latestdoi/154380587)


# M4 (NGC 6121) Proper Motion Membership Catalog

This repository contains the M4 membership catalog produced by 
[Wallace 2018](http://iopscience.iop.org/article/10.3847/2515-5172/aaf1a2/meta) 
from *Gaia* DR2 proper 
motion measurements, as well as the code used to create the catalog.  The contents of this repo 
are distributed under the BSD 3-Clause License (see [LICENSE](LICENSE.md) 
file for more information).

I would like to acknowledge the assistance of my academic advisors, 
Joel Hartman and Gaspar Bakos, as well as significant advice from 
[Waqas Bhatti](https://github.com/waqasbhatti).


## Dependencies

The following Python packages are used by the main code (pm_membership.py):
[astropy](http://www.astropy.org/), [matplotlib](https://matplotlib.org/), 
[numpy](http://www.numpy.org/),
[pickle](https://docs.python.org/3/library/pickle.html), 
[scikit-learn](http://scikit-learn.org/stable/).

Additionally, the code to query the *Gaia* archive (gaia_query.py)
uses [astropy] (http://www.astropy.org/)
 and [astroquery](https://astroquery.readthedocs.io/en/latest/).

## Repo Contents

* [pm_membership.py](pm_membership.py) -- the code used to generate the 
membership catalog, as well as produce a plot showing the fit results. The 
code takes two arguments to run: 1. the name of the input VOTable file 
following the *Gaia* data format, and 2. the name of the output pickle file
that stores the fitted Gaussian mixture model.

* [gaia_query.py](gaia_query.py) -- this code can be used to query the *Gaia* 
data archive for the appropriate data.  The call in the code currently 
reproduces the *Gaia* data I used to make the fitted model.

* gmm_fitted*.pkl -- the fitted models.  These are Python pickles of instances 
of an 
[sklearn.mixture.GaussianMixture](http://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html)
object.  The pickles were generated using Python 2, but can be read in to 
Python 3 with the `encoding='latin'` option to your `pickle.load()` call.
The first component of the output to the `predict_proba()` method is the
cluster membership probability (the second component is the field 
membership probability).

    * [gmm_fitted_maglessthan_19.0.pkl](gmm_fitted_maglessthan_19.0.pkl) -- the 
fitted model for objects with *Gaia* *G* magnitudes less (brighter) than 19.  

    * [gmm_fitted_maggreaterthan_19.0.pkl](gmm_fitted_maggreaterthan_19.0.pkl) -- 
the fitted model for objects with *Gaia* *G* magnitudes greater (dimmer) than 19

* [membership_catalog.txt](membership_catalog.txt) -- the catalog in plain 
text.  Contains *Gaia* DR2 ID numbers, our calculated membership probabilities, 
and some other information (*G* magnitudes, RA, dec, proper motion values and 
errors on the latter three). See header for column and unit information. 
**Please note that only objects brighter than *G*=19 are currently included
in this catalog file.**

## How to Use

**If all you want is the published catalog as-is,**
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
at [gaia_query.py](gaia_query.py), which shows the ADQL query to the *Gaia*
data archive for the input data for the membership calculation, and 
[pm_membership.py](pm_membership.py), which has code to fit the model, 
generate the membership catalog from the fitted model, and plot up
a visualization of the fitted model.

The code can be ran as-is as following to regenerate the membership catalog
as published:

```
# Download a portion of the *Gaia* data archive
python gaia_query.py foo.xml 

# Fit the model and outputs the fitted model, membership catalog, and some plots
python pm_membership.py foo.xml gmm_fitted 19
```

The first line gets the relevant data from the *Gaia* data archive and saves
it to the file `foo.xml`.  The second line then reads the data from `foo.xml`,
fits and saves the models to pickle files prepended `gmm_fitted`, and with
the model fits split between objects brighter than and dimmer than *G*=19.

See [pm_membership.py](pm_membership.py) for more detailed information on the
code itself.

## Caveats

* **Completeness** - In crowded regions such as globular clusters, the 
completeness of the *Gaia* DR2 source catalog does not match that obtained 
in other regions with a similar number of observations.  Near the cluster 
center, completeness of source detection is much less that of the survey
as a whole, reaching only ~50% for objects as bright as *G*=16.
For more details, please check out 
["Gaia Data Release 2: Catalogue Validation" (Arenou et al. 2018)](https://arxiv.org/pdf/1804.09375.pdf), Section 3.
In particular, Figure 6 shows completeness for M4 (NGC 6121), and Table B.1 
provides a quantification of this figure.


* **Change in fit as a function of magnitude** - Since errors in the proper
motion measurements increase with magnitude, and since the fraction of stars
that are cluster members is a function of magnitude, different values for
`mag_cutoff` (the third argument given to pm_membership.py) will produce
different fits.  A more complete proper motion membership catalog 
calculation would likely bin objects by magnitude and have a separate 
mixture model fit for each magnitude bin (more than the maximum of two
bins the code now provides for), but that is not implemented here.  
This change in fit as a function of magnitude is the primary reason
that the current catalog was chosen to cut off at *G*=19, as the 
fainter objects require some more care and consideration than the catalog
currently offers.

* **Change in fit depending on distance from cluster center** - Since
membership probability is additionally dependent on distance from 
cluster center---a given star near the center of the cluster is
much more likely to be a cluster member than a field star, without
knowing any other prior information, given the density of cluster
stars in the area, with this probability decreasing with increasing 
distance from the cluster center---how far from the cluster center
one goes to choose objects to include in the fit also affects
the fitted model.  We use objects within 30 arcminutes of the cluster
center in the code in this repo. We also tried 6 and 60 arcminutes and
found that the final fits did not change much between these three
data sets.  A more complete proper motion membership catalog would
perhaps bin objects into annuli centered on the cluster center
and having a separate mixture model fit for each annulus,
but that is not implemented here.



* **Initialization** - A Gaussian mixture model is an unsupervised machine 
learning model, meaning it finds structure in the data without providing 
classification labels for the detected structure.  Labels must be determined
afterward.  In the setup as currently exists in this repo, it is the *first*
component of the fitted model that corresponds to cluster membership 
probability.  The initial means for the fitting have been set such that this
should be the case for any M4-centered data set that is used for the fitting,
but we cannot make any absolute guarantee on this.  Please always double
check which component of the fit corresponds to the cluster membership
probability if adapting the code for your own purposes.  If the `initial_means`
or `random_state` arguments to the [`read_catalog_in_and_fit()`](https://github.com/joshuawallace/M4_pm_membership/blob/master/pm_membership.py#L46)
function in
pm_membership.py are altered, please take extra care on this front.  If you
use this code to create a proper motion membership catalog for a cluster
other than M4, you will probably need to change the value for `initial_means`
if you want to have any level of certainty as to which component corresponds
to the cluster membership probability.



## Use, Acknowledgment, and Attribution

The contents of this repo are distributed under the BSD 3-Clause License 
(see [LICENSE](LICENSE.md) file for more information).

The official publication for this catalog is 
["M4 Membership Catalog from *Gaia* Proper Motions"](http://iopscience.iop.org/article/10.3847/2515-5172/aaf1a2/meta),
authored by Joshua Wallace in 2018 and published in *Research Notes of the American 
Astronomical Society*.  
This repo also has a DOI at Zenodo.  The following Zenodo link
[10.5281/zenodo.1488301](https://doi.org/10.5281/zenodo.1488301) 
should resolve to the most recent version, while this link 
[10.5281/zenodo.1488302](https://zenodo.org/record/1488302#.W_w6g5y1thE) will resolve 
to version 1.0, the current officially released version.
**Any use of this catalog or the associated code should cite the Wallace 2018 
article (Bibtex entry below), 
and any use of the data or code in the repo should also cite the repo DOI.**



Bibtex entry for the Wallace 2018 article:
```
@article{2515-5172-2-4-213,
  author={Joshua J. Wallace},
  title={M4 Membership Catalog from Gaia Proper Motions},
  journal={Research Notes of the AAS},
  volume={2},
  number={4},
  pages={213},
  url={http://stacks.iop.org/2515-5172/2/i=4/a=213},
  year={2018},
  abstract={}
}
```