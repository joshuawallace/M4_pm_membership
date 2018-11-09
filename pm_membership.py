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

This code has some redunancies/unnecessary steps for a full end-to-end
calculation, but these were included to provide e.g. examples of how 
to read in the fitted model from a pickle file and use the contents.

Usage: pm_membership.py  input_catalog_filename  output_pickle_filename_prefix
       mag_cutoff

     input_catalog_filename: the filename for the input catalog, assumed to
          be a VOTable from the Gaia data archive or formatted similarly

     output_pickle_filename_prefix: the prefix for the filename for the pickle 
          file that will store the fitted Gaussian mixture model for later use

     mag_cutoff: the magnitude cutoff, where the model fitting gets split. 
          separate fits are made for objects brighter and dimmer (in Gaia G
          magnitude) than this value.  To not split the fitting, use 
          'None' or 'none'.
"""

from __future__ import absolute_import, division, print_function
import sys, pickle
import numpy as np
from astropy.io.votable import parse_single_table
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture as GM


def read_catalog_in_and_fit(input_filename, output_filename_prefix, 
                            mag_cutoff=None,
                            initial_means=[[-12.5,-19], [-3,-4]],
                            random_state=546):
    """
    Read in the Gaia data and fit the Gaussian mixture model

    input_filename - the input catalog, assume to be a VOTable and assumed
         to have similar format as the Gaia DR2 source catalog
    output_filename_prefix - the prefix of the name for the output pickle file, 
         which will store the fitted model
    mag_cutoff - the magnitude to split the model over.  Separate fits will
         be made for objects brighter and dimmer than this magnitude.
         A value of None will use all stars for a fit, no splitting
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

    Returns the fitted model and the read-in catalog data
    """

    # Read in the catalog
    table = parse_single_table(input_filename)
    catalog_data = table.array

    # Set up the Gaussian mixture model
    if mag_cutoff is not None:
        gm_model_bright = GM(n_components=2,max_iter=300,means_init=initial_means,
                             random_state=random_state)
        gm_model_dim = GM(n_components=2,max_iter=300,means_init=initial_means,
                          random_state=random_state)
    else:
        gm_model = GM(n_components=2,max_iter=300,means_init=initial_means,
                      random_state=random_state)
    
    # Mask out catalog entries without measured proper motions and dimmer
    # than the cutoff magnitude
    pm_mask_fitting = (~np.isnan(catalog_data['pmra'])) &\
        (~np.isnan(catalog_data['pmdec']))
    if mag_cutoff is not None:
        pm_mask_fitting_bright = pm_mask_fitting &\
            (catalog_data['phot_g_mean_mag'] < mag_cutoff)
        pm_mask_fitting_dim = pm_mask_fitting &\
            (catalog_data['phot_g_mean_mag'] >= mag_cutoff)
    
    # Extract the data to fit and fit it
    if mag_cutoff is not None:
        data_to_fit_bright = [item for item in zip(catalog_data['pmra'][pm_mask_fitting_bright],
                                                   catalog_data['pmdec'][pm_mask_fitting_bright])]
        data_to_fit_dim = [item for item in zip(catalog_data['pmra'][pm_mask_fitting_dim],
                                                   catalog_data['pmdec'][pm_mask_fitting_dim])] 
        gm_model_bright.fit(data_to_fit_bright)
        gm_model_dim.fit(data_to_fit_dim)
    else:
        data_to_fit = [item for item in zip(catalog_data['pmra'][pm_mask_fitting],
                                            catalog_data['pmdec'][pm_mask_fitting])]
        gm_model.fit(data_to_fit)


    # Dump the fitted model to a pickle file, and also return it
    if mag_cutoff is not None:
        with open(output_filename_prefix + "_maglessthan_" +\
                      str(mag_cutoff) + ".pkl","wb") as f:
            pickle.dump(gm_model_bright,f)
        with open(output_filename_prefix + "_maggreaterthan_" +\
                      str(mag_cutoff) + ".pkl","wb") as f:
            pickle.dump(gm_model_dim,f)
        return ((gm_model_bright,gm_model_dim), catalog_data)
    else:
        with open(output_filename_prefix + ".pkl","wb") as f:
            pickle.dump(gm_model,f)
        return (gm_model, catalog_data)


def read_model_in_and_use(model_filename_prefix, catalog_data, mag_cutoff=None):
    """
    Read the model in from the given pickle file and use it

    model_filename_prefix - the prefix for the name for the output 
         pickle file, which stores the fitted model
    catalog_data - the read-in data for the objects, assumed to be in
         similar format as the Gaia DR2 data.  See the 
         read_catalog_in_and_fit() function for an example for how to
         read in the catalog
    mag_cutoff - the magnitude the model was splitted over.  A value of 
         None means all stars were used in a single fit, no splitting
    """

    # Read in the fitted model
    if mag_cutoff is not None:
        with open(output_filename_prefix + "_maglessthan_" +\
                      str(mag_cutoff) + ".pkl","rb") as f:
            fitted_model_bright = pickle.load(f)
        with open(output_filename_prefix + "_maggreaterthan_" +\
                      str(mag_cutoff) + ".pkl","rb") as f:
            fitted_model_dim = pickle.load(f)
    else:
        with open(model_filename_prefix + ".pkl","rb") as f:
            fitted_model = pickle.load(f)
    
    # From the catalog data, get the proper motions and calculate probability
    # and create a dictionary matching ID to membership probability
    prob_dict = {}
    pm_mask_calculating = (~np.isnan(catalog_data['pmra'])) &\
        (~np.isnan(catalog_data['pmdec']))
    if mag_cutoff is not None:
        bright_mask = pm_mask_calculating & (catalog_data['phot_g_mean_mag'] < mag_cutoff)
        dim_mask = pm_mask_calculating & (catalog_data['phot_g_mean_mag'] >= mag_cutoff)

        # First, the bright
        data_to_calc_bright = [item for item in zip(catalog_data['pmra'][bright_mask],
                                                    catalog_data['pmdec'][bright_mask])]
        prob_bright = fitted_model_bright.predict_proba(data_to_calc_bright)
        prob_bright_first_comp = [item[0] for item in prob_bright] # Just the first component
        IDs_masked_bright = catalog_data['source_id'][bright_mask]
        for i in range(len(IDs_masked_bright)):
            prob_dict[IDs_masked_bright[i]] = prob_bright_first_comp[i]

        # Now, the dim
        data_to_calc_dim = [item for item in zip(catalog_data['pmra'][dim_mask],
                                                    catalog_data['pmdec'][dim_mask])]
        prob_dim = fitted_model_dim.predict_proba(data_to_calc_dim)
        prob_dim_first_comp = [item[0] for item in prob_dim] # Just the first component
        IDs_masked_dim = catalog_data['source_id'][dim_mask]
        for i in range(len(IDs_masked_dim)):
            prob_dict[IDs_masked_dim[i]] = prob_dim_first_comp[i]
    else:
        data_to_calc = [item for item in zip(catalog_data['pmra'][pm_mask_calculating],
                                             catalog_data['pmdec'][pm_mask_calculating])]
        prob = fitted_model.predict_proba(data_to_calc)
        prob_first_comp = [item[0] for item in prob] # Just the first component
        IDs_masked = catalog_data['source_id'][pm_mask_calculating]
        for i in range(len(IDs_masked)):
            prob_dict[IDs_masked[i]] = prob_first_comp[i]
               
    # And return the probability dict
    return prob_dict


def plot_up(catalog_data,membership_probabilities,
            mag_cutoff=None,
            plot_filename="gmm_fit.png"):
    """
    Show the membership probabilities graphically

    catalog_data - the read-in data for the objects, assumed to be in
         similar format as the Gaia DR2 data.  See the 
         read_catalog_in_and_fit() function for an example for how to
         read in the catalog
    membership_probabilities - a dictionary relating membership 
         probabilities to source IDs
    mag_cutoff - the magnitude the model was splitted over.  A value of 
         None means all stars were used in a single fit, no splitting,
         or to plot things up in one plot regardless of whether there
         was actually a mag_cutoff
    plot_filename - name of the file for the output plot
    """

    # Extract all the relevant information, assuming that the keys to 
    # membership_probabilities dictionary is the most complete list of what
    # we want to plot up, and then plot it up
    if mag_cutoff is None:
        pmra = catalog_data['pmra'][np.isin(catalog_data['source_id'],
                                            membership_probabilities.keys())]
        pmdec = catalog_data['pmdec'][np.isin(catalog_data['source_id'],
                                            membership_probabilities.keys())]
        source_ids_ordered = catalog_data['source_id'][np.isin(catalog_data['source_id'],
                                            membership_probabilities.keys())]
        probs = [membership_probabilities[ID] for ID in source_ids_ordered]


        # Plot up the points and make colorbar
        sc = plt.scatter(pmra,pmdec,c=probs,cmap='brg',s=1,alpha=.5)
        cbar = plt.colorbar(sc,label='Star cluster membership probability')

        # Add labels and save
        plt.xlabel("Proper motion in RA (mas/yr)")
        plt.ylabel("Proper motion in dec (mas/yr)")
        plt.xlim(-22,10)
        plt.ylim(-27,5)
        plt.tight_layout()
        plt.savefig(plot_filename,dpi=300)
        plt.close()

    else: # Split up plots by the magnitude cutoff
        IDs_bright = catalog_data['source_id'][catalog_data['phot_g_mean_mag'] < mag_cutoff]
        IDs_dim = catalog_data['source_id'][catalog_data['phot_g_mean_mag'] >= mag_cutoff]

        pmra_bright = catalog_data['pmra'][(np.isin(catalog_data['source_id'],
                                                   membership_probabilities.keys())) &
                                           (np.isin(catalog_data['source_id'],
                                                   IDs_bright))]
        pmra_dim = catalog_data['pmra'][(np.isin(catalog_data['source_id'],
                                                   membership_probabilities.keys())) &
                                           (np.isin(catalog_data['source_id'],
                                                   IDs_dim))]
        pmdec_bright = catalog_data['pmdec'][(np.isin(catalog_data['source_id'],
                                                   membership_probabilities.keys())) &
                                           (np.isin(catalog_data['source_id'],
                                                   IDs_bright))]
        pmdec_dim = catalog_data['pmdec'][(np.isin(catalog_data['source_id'],
                                                   membership_probabilities.keys())) &
                                           (np.isin(catalog_data['source_id'],
                                                   IDs_dim))]
        source_ids_bright_ordered = catalog_data['source_id'][(np.isin(catalog_data['source_id'],
                                                    membership_probabilities.keys())) &
                                                              (np.isin(catalog_data['source_id'],
                                                                       IDs_bright))]
        source_ids_dim_ordered = catalog_data['source_id'][(np.isin(catalog_data['source_id'],
                                                    membership_probabilities.keys())) &
                                                              (np.isin(catalog_data['source_id'],
                                                                       IDs_dim))]
        probs_bright = [membership_probabilities[ID] for ID in source_ids_bright_ordered]
        probs_dim = [membership_probabilities[ID] for ID in source_ids_dim_ordered]

        ## Bright
        # Plot up the points and make colorbar
        sc_bright = plt.scatter(pmra_bright,pmdec_bright,c=probs_bright,cmap='brg',s=1,alpha=.5)
        cbar_bright = plt.colorbar(sc_bright,label='Star cluster membership probability')

        # Add labels and save
        plt.xlabel("Proper motion in RA (mas/yr)")
        plt.ylabel("Proper motion in dec (mas/yr)")
        plt.xlim(-22,10)
        plt.ylim(-27,5)
        plt.tight_layout()
        last_period = plot_filename.rfind('.')
        bright_filename = plot_filename[:last_period] + "_brighter" + plot_filename[last_period:]
        plt.savefig(bright_filename,dpi=300)
        plt.close()

        ## Dim
        # Plot up the points and make colorbar
        sc_dim = plt.scatter(pmra_dim,pmdec_dim,c=probs_dim,cmap='brg',s=1,alpha=.5)
        cbar_dim = plt.colorbar(sc_dim,label='Star cluster membership probability')

        # Add labels and save
        plt.xlabel("Proper motion in RA (mas/yr)")
        plt.ylabel("Proper motion in dec (mas/yr)")
        plt.xlim(-22,10)
        plt.ylim(-27,5)
        plt.tight_layout()
        last_period = plot_filename.rfind('.')
        dim_filename = plot_filename[:last_period] + "_dimmer" + plot_filename[last_period:]
        plt.savefig(dim_filename,dpi=300)
        plt.close()
    

def produce_catalog(catalog_data,membership_probabilities,
                    output_catalog_filename="membership_catalog.txt"):
    """
    After calculating the probability memberships, write the 
    final output membership catalog.

    catalog_data - the read-in data for the objects, assumed to be in
         similar format as the Gaia DR2 data.  See the 
         read_catalog_in_and_fit() function for an example for how to
         read in the catalog
    membership_probabilities - a dictionary relating membership 
         probabilities to source IDs
    """

    # Get all the information together
    source_ids_ordered = catalog_data['source_id'][np.isin(catalog_data['source_id'],
                                        membership_probabilities.keys())]
    G_mag = catalog_data['phot_g_mean_mag'][np.isin(catalog_data['source_id'],
                                        membership_probabilities.keys())]
    ra = catalog_data['ra'][np.isin(catalog_data['source_id'],
                                        membership_probabilities.keys())]
    dec = catalog_data['dec'][np.isin(catalog_data['source_id'],
                                        membership_probabilities.keys())]
    ra_err = catalog_data['ra_error'][np.isin(catalog_data['source_id'],
                                        membership_probabilities.keys())]
    dec_err = catalog_data['dec_error'][np.isin(catalog_data['source_id'],
                                        membership_probabilities.keys())]
    pmra = catalog_data['pmra'][np.isin(catalog_data['source_id'],
                                        membership_probabilities.keys())]
    pmdec = catalog_data['pmdec'][np.isin(catalog_data['source_id'],
                                        membership_probabilities.keys())]
    pmra_err = catalog_data['pmra_error'][np.isin(catalog_data['source_id'],
                                        membership_probabilities.keys())]
    pmdec_err = catalog_data['pmdec_error'][np.isin(catalog_data['source_id'],
                                        membership_probabilities.keys())]
    probs = [membership_probabilities[ID] for ID in source_ids_ordered]

    # Write it out
    with open(output_catalog_filename,"w") as f:
        f.write("# Gaia_DR2_ID        Gaia_G_mag       memb_prob"+\
                    "              RA(deg)              RA_err(deg)"+\
                    "              dec(deg)             dec_err(deg)"+\
                    "         pm_RA(mas/yr)     pm_RA_err(mas/yr)"+\
                    "        pm_dec(mas/yr)    pm_dec_err(mas/yr)\n")
        formatter = "%-19s  %10s  %14.6e  %19s  %23.17e  %20s  %23.17e  "+\
            "%20s  %20s  %20s  %20s\n"
        for i in range(len(source_ids_ordered)):
            f.write(formatter % (source_ids_ordered[i], G_mag[i],
                                 probs[i], ra[i], ra_err[i]/(1000*3600), 
                                 dec[i], dec_err[i]/(1000*3600), pmra[i], 
                                 pmra_err[i], pmdec[i], pmdec_err[i]))


if __name__ == "__main__":
    # First, read in the desired input and output filenames
    if len(sys.argv) != 4:
        print ("Needed command line arguments: input_catalog_filename   " +\
                   "output_pickle_filename_prefix   mag_cutoff")
        raise RuntimeError("Wrong number of command line arguments.  Got " +\
                               str(len(sys.argv)-1) + " instead of 3")

    input_filename = sys.argv[1]
    output_filename_prefix = sys.argv[2]
    if sys.argv[3] not in ['None', 'none']:
        mag_cutoff = float(sys.argv[3])
    else:
        mag_cutoff = None

    # Next, read in the catalog and fit membership model
    fitted_model, catalog_data = read_catalog_in_and_fit(input_filename, 
                                                         output_filename_prefix,
                                                         mag_cutoff=mag_cutoff)

    # We can use the fitted_model directly:
    sample_proper_motions = [ [-12.5,-19],
                              [-3,-4],
                              [-13,-20],
                              [-10,-20],
                              [100,100]]
    if isinstance(fitted_model, tuple):
        probability_output = fitted_model[0].predict_proba(sample_proper_motions)
    else:
        probability_output = fitted_model.predict_proba(sample_proper_motions)
    print([item[0] for item in probability_output]) # Just the first component

    # Or we can use the function above:
    full_prob_output = read_model_in_and_use(output_filename_prefix,
                                             catalog_data,mag_cutoff=mag_cutoff)

    # Let's plot it up
    plot_up(catalog_data,full_prob_output,mag_cutoff=mag_cutoff)

    # And now write out the plain text catalog
    produce_catalog(catalog_data,full_prob_output)
