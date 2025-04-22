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

def read_indexed_stars(filename):
    selected_indices = np.loadtxt(filename, dtype=int)
    
    return selected_indices
        

def get_selected_stars_txt(galaxy, snap):
    # Format the expected pattern
    directory = f"{galaxy}_Phase_Plots/"
    pattern = f"{galaxy}_{snap}_SelectedParticles_tidal*.txt"

    # List all files and filter with fnmatch
    txt_files = [f for f in os.listdir(f"{galaxy}_Phase_Plots") if fnmatch.fnmatch(f, pattern)]
    
    # Append directory to filename to get path
    for i in range(len(txt_files)):
        txt_files[i] = directory + txt_files[i]
    
    return txt_files

galaxies = np.array(["MW","M31"])
snaps = np.array([275,280,285,335,340,345,410,415,420,425,430,435,470,475,480,801])


all_events = {}
key_snaps = {}

for galaxy in galaxies:
    snaps_important = []
    for snap in snaps:
        selected = []
        for file in get_selected_stars_txt(galaxy, snap):
            particles = read_indexed_stars(file)
            for i in particles:
                if i not in selected:
                    selected.append(i)
                    
        if len(selected) != 0:
            all_events[f"{galaxy}_{snap}"] = np.array(selected)
            snaps_important.append(snap)

    key_snaps[galaxy] = np.array(snaps_important)

for key in key_snaps.keys():
    snaps = key_snaps[key]
    for i in range(len(snaps)-1,-1,-1):
        ax, rn = Density_plot_stars(key, snaps[i])
        for i in range(i,-1,-1):
            mask = all_events[f"{key}_{snaps[i]}"]
            ax.scatter(rn[mask,0],rn[mask,1],marker='.',s=1,label=f'Snap: {snaps[i]}')
        ax.legend()
    
    ax, rn = Density_plot_stars(key, 801)
    for i in range(len(snaps)-1,-1,-1):
        mask = all_events[f"{key}_{snaps[i]}"]
        ax.scatter(rn[mask,0],rn[mask,1],marker='.',s=1,label=f'Snap: {snaps[i]}')
    ax.legend()
        
        
    
    