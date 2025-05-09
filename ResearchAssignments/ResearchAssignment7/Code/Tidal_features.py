#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 14:42:26 2025

@author: petershea
"""
import os
import numpy as np
import fnmatch
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from Plotting_Snaps import Density_plot_stars
from MassProfile import MassProfile
from CenterOfMass2 import CenterOfMass
from CenterOfMass2 import COMdefine

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
galaxies = np.array(["MW","M31"]) # "MW"
snaps = np.array([275, 280, 285, 335, 340, 345, 410, 415, 420, 425, 430, 435, 470, 475, 480, 801])

# Initialize dicts
all_events = {}  # particle indices
key_snaps = {}   # snaps with tidal features identified

for galaxy in galaxies:
    snaps_important = np.array([], dtype=int)  # Start as empty array for snaps with identified tidal features
    for snap in snaps:
        selected = np.array([], dtype=int)  # Empty array for particle indices
        for file in get_selected_stars_txt(galaxy, snap):
            particles = read_indexed_stars(file)  # This returns a numpy array
            selected = np.union1d(selected, particles)  # Keep only unique entries
        
        if selected.size != 0: # Only adds events if there was a tidal feature detected in the snap
            all_events[f"{galaxy}_{snap}"] = selected
            snaps_important = np.append(snaps_important, snap)

    key_snaps[galaxy] = snaps_important 

# Important variables for COM code
delta = 0.1

volDec_M33 = 4
volDec = 2



for galaxy in galaxies:
    print(galaxy)
    mpersnap = np.empty((0, len(key_snaps[galaxy]))) # array to store the mass fraction of each tidal feature at each snap
    for snap in snaps:
        print(snap)
        MP = MassProfile(galaxy, snap) # creates mass profile 
        
        # Defines the mass percent whcih an enclosing radius will be found for
        if galaxy == "MW":
            percent = 0.95
        elif galaxy == "M31":
            percent = 0.99
        
        r_crit = MP.Rad_MassEnclosed(2, percent) # defines a critical radius using Rad_MassEnclosed function
        
        # create center of mass object
        COMD = CenterOfMass(f"{galaxy}/{galaxy}_{snap}.txt",2)
        
        # Compute COM using disk particles
        if galaxy != 'M33':
            COMP = COMD.COM_P(volDec,delta)
        # Uses volDec_M33 = 4 for M33
        else:
            COMP = COMD.COM_P(volDec_M33,delta)

        # Determine positions of disk particles relative to COM 
        xD = COMD.x - COMP[0].value 
        yD = COMD.y - COMP[1].value 
        zD = COMD.z - COMP[2].value 
        
        m = COMD.m
        
        # total magnitude to be used with jacobi radius to determine if particle is beyond
        rtot = np.sqrt(xD**2 + yD**2 + zD**2)
        
        mpers = np.empty(len(key_snaps[galaxy])) # for each tidal feature a mass fraction determined and stored in this array
        
        for idx, snap_key in enumerate(key_snaps[galaxy]):
            print(f'Tidal_{snap_key}')
            mask = all_events[f'{galaxy}_{snap_key}'] # creates a mask using the stored particle indices
            rtot_masked = rtot[mask] # masks all except particles beyond the critical radius
            frac = np.count_nonzero(rtot_masked >= r_crit) / len(rtot_masked) # determines the mass fraction
            if snap_key <= snap: # only stores it in an array if the snap is after the tidal feature is identified
                mpers[idx] = frac
            else: mpers[idx] = np.nan

        # Stack new row
        mpersnap = np.vstack((mpersnap, mpers))
    
    # ploting the tidal structure evolution over time 
    fig, ax = plt.subplots(dpi=1000)
    
    for i in range(len(mpersnap[0])-1,-1,-1):
        ax.plot(snaps,mpersnap[:,i], label=f'Tidal Feature: {key_snaps[galaxy][i]}')
        
    ax.set_xlabel('Simulation Snapshot')
    ax.set_ylabel('Fraction of Tidal Mass')
    ax.set_title(f'Time Evolution of {galaxy} Tidal Features')
    ax.legend()

    plt.show()
        




for key in key_snaps.keys(): # loops over galaxies
    snaps = key_snaps[key] # defines important snaps for given galaxy
    for i in range(len(snaps)-1,-1,-1): # iterates backwards through snaps creating density plot and scattering selected particles
        ax, rn = Density_plot_stars(key, snaps[i]) # plots the face on column density 
        MP = MassProfile(key, snaps[i]) # creates a mass profile
        for j in range(i,-1,-1): # scatter plots all priorly identified tidal particles
            mask = all_events[f"{key}_{snaps[j]}"]
            ax.scatter(rn[mask,0],rn[mask,1],marker='.',s=3,label=f'Snap: {snaps[j]}')
        
        # determines percent of galaxy mass enclosed by the desired radius
        if key == "MW":
            percent = 0.95
        elif key == "M31":
            percent = 0.99
        
        R_crit = MP.Rad_MassEnclosed(2, percent) # determines the critical radius 
        
        circle = patches.Circle((0,0), R_crit, edgecolor='k', facecolor='none') #plots a circle 
        ax.add_patch(circle)
        
        ax.legend()
        
        plt.show()
        
    # does the same but for the last snapshot 801
    ax, rn = Density_plot_stars(key, 801) # plots the end simulation state
    for i in range(len(snaps)-1,-1,-1):
        mask = all_events[f"{key}_{snaps[i]}"]
        ax.scatter(rn[mask,0],rn[mask,1],marker='.',s=3,label=f'Snap: {snaps[i]}')
    MP = MassProfile(key, 801)
     
    R_crit = MP.Rad_MassEnclosed(2, percent)
     
    circle = patches.Circle((0,0), R_crit, edgecolor='k', facecolor='none')
    ax.add_patch(circle)
     
    ax.legend()

    plt.show()
   
    
    