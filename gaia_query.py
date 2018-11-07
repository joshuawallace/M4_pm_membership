#!/usr/bin/env python

"""
gaia_query.py, code to pull the data from the Gaia data archive.

Author: Joshua Wallace, advised by Joel Hartman and Gaspar Bakos.
Also special thanks to Semyeong Oh for helping me with this and 
overall helping me understand the Gaia data archive.

The code as written here will extract the data that were used to 
create the fitted model (gmm_fitted.pkl) as presently in the
directory.

Usage: gaia_query.py  output_catalog_filename

     output_catalog_filename: the filename for the VOTable output
          of this Astroquery call
"""

from astroquery.gaia import Gaia
import sys
from astropy.io.votable import writeto


if len(sys.argv) != 2:
    print ("Needed command line arguments: output_catalog_filename ")
    raise RuntimeError("Wrong number of command line arguments.  Got " +\
                           str(len(sys.argv)-1) + " instead of 1")

output_filename = sys.argv[1]

# Create the query string, following ADQL format
# The coordinates are from http://adsabs.harvard.edu/abs/2010AJ....140.1830G
query_string =  "SELECT source_id, phot_g_mean_mag, ra, dec, ra_error, "
query_string += "dec_error, pmra, pmdec, pmra_error, pmdec_error"
query_string += "FROM gaiadr2.gaia_source "
query_string += "WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',245.89675000000003,-26.52575,0.5))=1"

# Launch the query and get reuslts
table = Gaia.launch_job(query_string).get_results()


# Write the table to output
writeto(table,output_filename)
