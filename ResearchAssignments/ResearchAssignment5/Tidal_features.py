#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 14:42:26 2025

@author: petershea
"""
import os
import numpy as np
import fnmatch

from Plotting_Snaps import Density_plot_stars

"""
This code is meant to take text files of indexed stars created in star_select_visualizer.py 
and plot the data over a face on density map of each galaxy. It plots each particle
identified to be a part of a tidal feature in the same or earlier snapshots of the simulation.
"""

def read_indexed_stars(filename):
    '''function to read in selected particles from file
    
    Parameters
    ----------
    filename : string
        name of the file containing indexed particles identified in tidal features

    Returns
    -------
    selected_indices : list
        list disk particle indicies
    '''
    selected_indices = np.loadtxt(filename, dtype=int)
    
    return selected_indices
        

def get_selected_stars_txt(galaxy, snap):
    '''Function to get paths to text files containing indexed particles
    
    Parameters
    ----------
    galaxy : string
        name of the galaxy in question: 'MW' or 'M31'
    snap : int
        snap number indicating simulation snapshot of interest

    Returns
    -------
    txt_files : list
        list of path strings totext files containing indexed particles
    '''
    # Format the expected pattern
    directory = f"{galaxy}_Phase_Plots/"
    pattern = f"{galaxy}_{snap}_SelectedParticles_tidal*.txt"

    # List all files and filter with fnmatch
    txt_files = [f for f in os.listdir(f"{galaxy}_Phase_Plots") if fnmatch.fnmatch(f, pattern)]
    
    # Append directory to filename to get path
    for i in range(len(txt_files)):
        txt_files[i] = directory + txt_files[i]
    
    return txt_files

# defines arrays to loop over for all galaxy and snaps investigated
galaxies = np.array(["MW","M31"])
snaps = np.array([275,280,285,335,340,345,410,415,420,425,430,435,470,475,480,801])

# defines dicts to store particle indexes and which snaps were identified to produce these features
all_events = {} # particle indexes
key_snaps = {}  # snaps with tidal features identified

for galaxy in galaxies: # loops over galaxies
    snaps_important = [] # list to store snaps where tidal features were identified
    for snap in snaps: # loops over all snaps
        selected = [] # list to store particle indices 
        for file in get_selected_stars_txt(galaxy, snap): # loops over all files for a given snap
            particles = read_indexed_stars(file) # reads in the data
            for i in particles: # loops over indexed particles
                if i not in selected: # only appends a given particle once for a given snap
                    selected.append(i) 
                    
        if len(selected) != 0: # only appends snap and array of indicies if tidal feature identified
            all_events[f"{galaxy}_{snap}"] = np.array(selected)
            snaps_important.append(snap)

    key_snaps[galaxy] = np.array(snaps_important) # stores snaps in dict

for key in key_snaps.keys(): # loops over galaxies
    snaps = key_snaps[key] # defines important snaps for given galaxy
    for i in range(len(snaps)-1,-1,-1): # iterates backwards through snaps creating density plot and scattering selected particles
        ax, rn = Density_plot_stars(key, snaps[i])
        for i in range(i,-1,-1):
            mask = all_events[f"{key}_{snaps[i]}"]
            ax.scatter(rn[mask,0],rn[mask,1],marker='.',s=1,label=f'Snap: {snaps[i]}')
        ax.legend()
    
    ax, rn = Density_plot_stars(key, 801) # plots the end simulation state
    for i in range(len(snaps)-1,-1,-1):
        mask = all_events[f"{key}_{snaps[i]}"]
        ax.scatter(rn[mask,0],rn[mask,1],marker='.',s=1,label=f'Snap: {snaps[i]}')
    ax.legend()
        
        
    
    