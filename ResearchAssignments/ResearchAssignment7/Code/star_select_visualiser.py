#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 11:30:29 2025

@author: petershea
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.widgets import RectangleSelector


import numpy as np

from ReadFile import Read
from CenterOfMass2 import CenterOfMass
from MassProfile import MassProfile
from OrbitCOM import readOrbit
from OrbitCOM import relative
from SnapIDs import find_close_encounter_IDs
from Plotting_Snaps import RotateFrame
from Plotting_Snaps import Density_plot_Auto


"""
This code is a more functional version of the Plotting_Snaps code lacking labels 
and titles as ploting in %matplotlib qt does not handle multiple subplots well
and is required for interactive plotting. Further, it also shows both phase and 
density plots together and will scatter selected particles on the other plot to
see which particles have been selected and how the plots relate.
"""

# Important variablels for COM code
delta = 0.1

volDec_M33 = 4 # not needed but included
volDec = 2


def Star_select_visualizer_p(galaxy, snap):
    '''This function takes a single galaxy and snap number to create a face on
    phase plot and density map of the galaxy at that snap number. The plot 
    created needs to be created in a window and is interactive. Dragging a box 
    within the phase plot selects particles within it, writes their indicies to 
    a txt file, and plots all identified particles in a scatter plot ovelaying
    the face on column density 
    
    Parameters:
        galaxy: string; a string for which galaxies to create a plot for ex. 
            "MW" or "M31"
        snap: int; an integer for the snap id of the simulation state
    
    Returns: 
        None
    '''
    # Create a COM of object for given galaxy
    COMD = CenterOfMass(f"{galaxy}/{galaxy}_{snap}.txt",2)
    
    # Compute COM using disk particles
    if galaxy != 'M33':
        COMP = COMD.COM_P(volDec,delta)
        # Uses volDec_M33 = 4 for M33
    else:
        COMP = COMD.COM_P(volDec_M33,delta)
    
    COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])
    
    # Determine positions of disk particles relative to COM 
    xD = COMD.x - COMP[0].value 
    yD = COMD.y - COMP[1].value 
    zD = COMD.z - COMP[2].value 
    
    # total magnitude
    rtot = np.sqrt(xD**2 + yD**2 + zD**2)
    
    # Determine velocities of disk particles relatiev to COM motion
    vxD = COMD.vx - COMV[0].value 
    vyD = COMD.vy - COMV[1].value 
    vzD = COMD.vz - COMV[2].value 
    
    # total velocity 
    vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)
    
    # Arrays for r and v 
    r = np.array([xD,yD,zD]).T # transposed 
    v = np.array([vxD,vyD,vzD]).T
    
    # compute the rotated position and velocity vectors
    rn, vn = RotateFrame(r,v)
    
    ### Plotting
    fig, ax = plt.subplots(1,2,figsize=(22,10),dpi=300)
    
    # Plot 2D Histogram for one component of  Pos vs Vel 
    # ADD HERE
    ax[0].hist2d(rn[:,0],vn[:,1], bins=800, norm=LogNorm(), cmap='plasma')

    fig.suptitle(f'{galaxy}; Snap={snap}')
    fig.supxlabel('x [kpc]')
    ax[0].set_ylabel('yv [km/s]')
    ax[1].set_ylabel('y [kpc]')
    
    #set axis limits
    ax[0].set_xlim(-150,150)
    ax[0].tick_params(axis='x', labelsize=4)
    ax[0].tick_params(axis='y', labelsize=4)
    
    ax[1].hist2d(rn[:,0],rn[:,1], bins=800, norm=LogNorm(), cmap='plasma')

    #set axis limits
    ax[1].set_ylim(-150,150)
    ax[1].set_xlim(-150,150)
    ax[1].tick_params(axis='x', labelsize=4)
    ax[1].tick_params(axis='y', labelsize=4)
    
    ax[1].set_aspect('equal')
    
    selected_particles = np.array([])

    def on_select(eclick, erelease):
        '''This function interfaces with the interactive matplotlib plot 
        selecting particles within a box defined by clicking and dragging on 
        the plot
        
        Parameters:
            eclick: button input; pressing and holding the mouse button
            erelease: button input; releasing the held mouse button
        
        Returns:
            None
        '''
        global selected_particles  
        x_min, x_max = sorted([eclick.xdata, erelease.xdata])
        y_min, y_max = sorted([eclick.ydata, erelease.ydata])

        # Find all particles within the selected region
        mask = (rn[:, 0] >= x_min) & (rn[:, 0] <= x_max) & (vn[:, 1] >= y_min) & (vn[:, 1] <= y_max)
        selected_particles = np.where(mask)[0]  # Get the indices of selected particles
        selected_rn = rn[mask]

        # Defines filename to 
        filename = f"{galaxy}_Phase_Plots/{galaxy}_{snap}_SelectedParticles_tidal.txt"
    
        # Save the selected indices to a text file
        np.savetxt(filename, selected_particles, fmt="%d")
        
        # Prints number of selected particles
        print(f'Selected {len(selected_particles)} particles')

        while len(ax[1].collections) > 1:
            ax[1].collections[-1].remove()
            
        scatter_plot = ax[1].scatter(selected_rn[:, 0], selected_rn[:, 1], marker='.', c='k', s=1)
        
        fig.canvas.draw()
    
    rect_selector = RectangleSelector(ax[0], on_select, useblit=True, button=[1], interactive=True)
    rect_selector.set_active(True)
  
    plt.tight_layout()
    plt.show(block=True)

   
    
def Star_select_visualizer_d(galaxy, snap):
    '''This function takes a single galaxy and snap number to create a face on
    phase plot and density map of the galaxy at that snap number. The plot 
    created needs to be created in a window and is interactive. Dragging a box 
    within the density plot selects particles within it, writes their indicies to 
    a txt file, and plots all identified particles in a scatter plot ovelaying
    the phase diagram
    
    Parameters:
        galaxy: string; a string for which galaxies to create a plot for ex. 
            "MW" or "M31"
        snap: int; an integer for the snap id of the simulation state
    
    Returns: 
        None
    '''
    # Create a COM of object for given galaxy
    COMD = CenterOfMass(f"{galaxy}/{galaxy}_{snap}.txt",2)
    
    # Compute COM using disk particles
    if galaxy != 'M33':
        COMP = COMD.COM_P(volDec,delta)
        # Uses volDec_M33 = 4 for M33
    else:
        COMP = COMD.COM_P(volDec_M33,delta)
    
    COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])
    
    # Determine positions of disk particles relative to COM 
    xD = COMD.x - COMP[0].value 
    yD = COMD.y - COMP[1].value 
    zD = COMD.z - COMP[2].value 
    
    # total magnitude
    rtot = np.sqrt(xD**2 + yD**2 + zD**2)
    
    # Determine velocities of disk particles relatiev to COM motion
    vxD = COMD.vx - COMV[0].value 
    vyD = COMD.vy - COMV[1].value 
    vzD = COMD.vz - COMV[2].value 
    
    # total velocity 
    vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)
    
    # Arrays for r and v 
    r = np.array([xD,yD,zD]).T # transposed 
    v = np.array([vxD,vyD,vzD]).T
    
    # compute the rotated position and velocity vectors
    rn, vn = RotateFrame(r,v)
    
    ### Plotting
    fig, ax = plt.subplots(1,2,figsize=(22,10),dpi=300)
    
    # Plot 2D Histogram for one component of  Pos vs Vel 
    # ADD HERE
    ax[0].hist2d(rn[:,0],vn[:,1], bins=800, norm=LogNorm(), cmap='plasma')

    #set axis limits
    ax[0].set_xlim(-150,150)
    ax[0].set_xticklabels([])  # Remove x-axis tick labels
    ax[0].set_yticklabels([])
    
    ax[1].hist2d(rn[:,0],rn[:,1], bins=800, norm=LogNorm(), cmap='plasma')

    #set axis limits
    ax[1].set_ylim(-150,150)
    ax[1].set_xlim(-150,150)
    ax[1].set_xticklabels([])  # Remove x-axis tick labels
    ax[1].set_yticklabels([])
    
    ax[1].set_aspect('equal')
    
    selected_particles = []

    def on_select(eclick, erelease):
        '''This function interfaces with the interactive matplotlib plot 
        selecting particles within a box defined by clicking and dragging on 
        the plot
        
        Parameters:
            eclick: button input; pressing and holding the mouse button
            erelease: button input; releasing the held mouse button
        
        Returns:
            None
        '''
        global selected_particles  
        x_min, x_max = sorted([eclick.xdata, erelease.xdata])
        y_min, y_max = sorted([eclick.ydata, erelease.ydata])

        # Find all particles within the selected region
        mask = (rn[:, 0] >= x_min) & (rn[:, 0] <= x_max) & (rn[:, 1] >= y_min) & (rn[:, 1] <= y_max)
        selected_particles = np.where(mask)[0]  # Get the indices of selected particles
        selected_rn = rn[mask]
        selected_vn = vn[mask]

        # Defines filename to 
        filename = f"{galaxy}_Phase_Plots/{galaxy}_{snap}_SelectedParticles_tidal.txt"
    
        # Save the selected indices to a text file
        np.savetxt(filename, selected_particles, fmt="%d")
        
        # Prints number of selected particles
        print(f'Selected {len(selected_particles)} particles')
        
        # removes previous scatter plot
        while len(ax[0].collections) > 1:
            ax[0].collections[-1].remove()
        
        # plots selected particles on then other subplot
        scatter_plot = ax[0].scatter(selected_rn[:, 0], selected_vn[:, 1], marker='.', c='k', s=1)
        
        fig.canvas.draw()
    
    rect_selector = RectangleSelector(ax[1], on_select, useblit=True, button=[1], interactive=True)
    rect_selector.set_active(True)
  
    plt.show(block=True)
    
#Star_select_visualizer_p("MW", 335)