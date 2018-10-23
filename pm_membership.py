#!/usr/bin/env python

"""
pm_membership.py, code to fit and use a model fit to proper motions to 
determine cluster membership for M4.

Author: Joshua Wallace, advised by Joel Hartman and Gaspar Bakos

If you use this code, please cite (fill in information after publication)

While this code is focused on M4, it can be modified to work on other
clusters or moving groups, and may work with no other alternation than 
pointing the code to the correct input file.
"""

from __future__ import absolute_import, division, print_function
import sys, pickle
from astropy.io.votable import parse_single_table
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture as GM


# First, read in the desired input and output file names
