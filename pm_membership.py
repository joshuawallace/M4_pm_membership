#!/usr/bin/env python

"""
pm_membership.py, code to fit and use a model fit to proper motions to 
determine cluster membership for M4.

Author: Joshua Wallace, advised by Joel Hartman and Gaspar Bakos

If you use this code, please cite (fill in information after publication)

While this code is focused on M4, it can be modified to work on other
clusters or moving groups, and may work with no other alternation than 
pointing the code to the correct input file.

The code is segmented into functions for reading in and fitting the catalog
and reading in and using the fitted model.  If you would like to use the
model I already fitted, you will only need the latter functionality.

Usage: pm_membership.py input_catalog_filename output_pickle_filename

     input_catalog_filename: the filename for the input catalog, assumed to
          be a VOTable from the Gaia data archive or formatted similarly

     output_pickle_filename: the filename for the pickle file that will store
          the fitted Gaussian mixture model for later use

"""

from __future__ import absolute_import, division, print_function
import sys, pickle
import numpy as np
from astropy.io.votable import parse_single_table
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture as GM


def read_catalog_in_and_fit(input_filename, output_filename, mag_cutoff=19.3,
                            initial_means=[[-12.5,-19], [-3,-4]],
                            random_state=546):
    """
    input_filename - the input catalog, assume to be a VOTable and assumed
         to have similar format as the Gaia DR2 source catalog
    output_filename - the name for the output pickle file, which will store
         the fitted model
    mag_cutoff - the dimmest (maximum) G-band magnitude for stars to use 
         in the model fit
    initial_means - the initial guesses for the means of the two Gaussian
         components.  Random initializations (the default) lead to random
         assignments: sometimes the first component is what gets fit to 
         the cluster proper motions, other times its the second component.
         For consistency, the default initial_means provided here is intended
         so that the first component fits the cluster proper motions, which
         have a mean value around the default value shown above.
    random_state - here we have provided a fixed seed value for the random
         number generator for consistency across runs.  Changing this value
         should have little to no effect on the final fit.

    Returns the fitted model and the proper motions
    """

    # Read in the catalog
    table = parse_single_table(input_filename)
    catalog_data = table.array

    # Set up the Gaussian mixture model
    gm_model = GM(n_components=2,max_iter=300,means_init=initial_means,
                  random_state=random_state)
    
    # Mask out catalog entries without measured proper motions and dimmer
    # than the cutoff magnitude
    pm_mask = (catalog_data['phot_g_mean_mag'] < mag_cutoff) & \
        (~np.isnan(catalog_data['pmra'])) & (~np.isnan(catalog_data['pmdec']))
    
    # Extract the data to fit and fit it
    data_to_fit = [item for item in zip(catalog_data['pmra'][pm_mask],
                                        catalog_data['pmdec'][pm_mask])]
    gm_model.fit(data_to_fit)

    # Dump the fitted model to a pickle file, and also return it
    with open(output_filename,"wb") as f:
        pickle.dump(gm_model,f)
    return (gm_model, (catalog_data['pmra'][pm_mask],
                       catalog_data['pmdec'][pm_mask]))


def read_model_in_and_use(model_filename,proper_motions):
    """
    model_filename - the name for the output pickle file, which stores
         the fitted model
    proper_motions - list/array of proper motions to be fed into the model,
         of format [ [pm_ra_1, pm_dec_1], [pm_ra_2, pm_dec_2], ...]
    """

    # Read in 
    with open(model_filename,"rb") as f:
        fitted_model = pickle.load(f)
    
    # Output the membership probability from the proper motions
    prob = fitted_model.predict_proba(proper_motions)

    # And return just the first component's membership probability
    return np.array([item[0] for item in prob])


def plot_up(proper_motions,membership_probabilities,
            plot_filename="gmm_fit.png"):
    """
    proper_motions - the proper motions [pm_ra_list, pm_dec_list]
    membership_probabilities - the membership probabilities associated
         with the proper motions above
    plot_filename - name of the file for the output plot
    """
    
    # Plot up the points and make colorbar
    sc = plt.scatter(proper_motions[0],proper_motions[1],
                     c=membership_probabilities,cmap='brg',
                     s=1,alpha=.6)
    cbar = plt.colorbar(sc,label='Star cluster membership probability')

    # Add labels and save
    plt.xlabel("Proper motion in RA (mas/yr)")
    plt.ylabel("Proper motion in dec (mas/yr)")
    plt.xlim(-22,10)
    plt.ylim(-27,5)
    plt.tight_layout()
    plt.savefig(plot_filename,dpi=300)
    plt.close()
    


if __name__ == "__main__":
    # First, read in the desired input and output filenames
    if len(sys.argv) != 3:
        print ("Needed command line arguments: input_catalog_filename " +\
                   "output_pickle_filename")
        raise RuntimeError("Wrong number of command line arguments.  Got " +\
                               str(len(sys.argv)-1) + " instead of 2")

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]


    # Next, read in the catalog and fit membership model
    fitted_model, proper_motions = read_catalog_in_and_fit(input_filename, 
                                                           output_filename)
    # We can use the fitted_model directly:
    sample_proper_motions = [ [-12.5,-19],
                              [-3,-4],
                              [-13,-20],
                              [-10,-20],
                              [100,100]]
    probability_output = fitted_model.predict_proba(sample_proper_motions)
    print([item[0] for item in probability_output]) # Just the first component

    # Or we can use the function above:
    function_probability_output = read_model_in_and_use(output_filename,
                                                        sample_proper_motions)
    print(function_probability_output)

    # Let's plot it up
    proper_motion_tofit = [item for item in zip(proper_motions[0],
                                                proper_motions[1])]
    full_prob_output = read_model_in_and_use(output_filename,
                                             proper_motion_tofit)
    plot_up(proper_motions,full_prob_output)
