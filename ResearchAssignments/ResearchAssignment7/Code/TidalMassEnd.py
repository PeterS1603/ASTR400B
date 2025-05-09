#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 22:58:24 2025

@author: petershea
"""

# Imports

#from Plotting_Snaps import RotateFrame
#from CenterOfMass2 import CenterOfMass
#from MassProfile import MassProfile
#import matplotlib.pyplot as plt
#from Tidal_features import read_indexed_stars
#from Tidal_featires import get_selected_stars_txt

########## Begin Code ##########
'''
This will be used to determine where the majority of the mass from these features 
is at a given point in the simulation allowing a plot to be created detailing 
the existance of these tidal features. Further, once most of the mass has returned 
to be inside this critical radius the feature can be assumed to no longer exist 
and such lifetime of these features can be estimated
'''
### Plot 2 - selected particles within ~99.5%? mass radius of galaxy ###
def MassPercentHist(galaxy,snap):
    '''Function to plot percent of each tidal features disk particles within a 
    key radius of the COM of the host galaxy
       
    Parameters
    ----------
    galaxy : string
        name of the galaxy in question: 'MW' or 'M31'
    snap : int
        snap number indicating simulation snapshot of interest

    Returns
    -------
    None
    
    '''
    
    # read in selected particles for each snap
    
    # create mass profile for merger remnant
    
    # determine radius in which ~99.5%? of merger mass resides
    
    # determine particles within and outside of this radius
    
    # determines mass ratio for tidal features created at each snap
    
    # plots as a histogram for each snap
    
    