#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 12:40:07 2025

@author: petershea
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as patches


import numpy as np
import astropy.units as u

from ReadFile import Read
from CenterOfMass2 import CenterOfMass
from MassProfile import MassProfile
from OrbitCOM import readOrbit
from OrbitCOM import relative
from SnapIDs import find_close_encounter_IDs

# for contours
import scipy.optimize as so


def RotateFrame(posI,velI):
    """a function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
    
    PARAMETERS
    ----------
        posI : `array of floats`
             3D array of positions (x,y,z)
        velI : `array of floats`
             3D array of velocities (vx,vy,vz)
             
    RETURNS
    -------
        pos: `array of floats`
            rotated 3D array of positions (x,y,z) 
            such that disk is in the XY plane
        vel: `array of floats`
            rotated 3D array of velocities (vx,vy,vz) 
            such that disk angular momentum vector
            is in the +z direction 
    """
    
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    
    # normalize the angular momentum vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to
    # z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel

def Greatest_Jacobi_Radius(R_12,R_13,M2,M3,M1):
    
    # calculates jacobi radius for both otehr objects
    R_j12 = R_12 * (M2 / (2*M1))**(1/3)
    R_j13 = R_13 * (M3 / (2*M1))**(1/3)
    
    if R_j12 > R_j13:
        return R_j12
    else:
        return R_j13

delta = 0.1

volDec_M33 = 4
volDec = 2

# total masses fo MW, M31, and M33 from HW4
M_MW = 2.06e12
M_M31 = 2.06e12
M_M33 = 0.196e12

# Read in COM orbit data
MW_orbit = readOrbit('Orbits/Orbit_MW.txt')
M31_orbit = readOrbit('Orbits/Orbit_M31.txt')
M33_orbit = readOrbit('Orbits/Orbit_M33.txt') 

MW_M31_relPos = relative(MW_orbit[1],M31_orbit[1])
M33_M31_relPos = relative(M31_orbit[1],M33_orbit[1])

def Density_plot(galaxy, snap):
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
    
    # Hook for implementing Jacobi radius
    """
    R_MW_M31 = MW_M31_relPos[int(snap/5)]
    R_M33_M31 = M33_M31_relPos[int(snap/5)]
    
    R_j = Greatest_Jacobi_Radius(R_MW_M31, R_M33_M31, M_MW, M_M33, M_M31)
    print(f'Jacobi Radius = {R_j} [kpc]')
    """
    # Rotated Disk - FACE ON
    # Disk Density 
    fig, ax = plt.subplots(figsize=(11, 10))
                  # x      y
    ax.hist2d(rn[:,0],rn[:,1], bins=800, norm=LogNorm(), cmap='plasma')
    
    
    #   Hook for implementing Jacobi radius
    #circle = patches.Circle((0,0), R_j, edgecolor='k', facecolor='none')
    # Add the circle to the axes
    #ax.add_patch(circle)
    
    # Add axis labels
    ax.set_xlabel('x [kpc]', fontsize=22)
    ax.set_ylabel('y [kpc]', fontsize=22)
    ax.set_title(f'{galaxy} Face On Density; Snap={snap}', fontsize=28)
    
    #set axis limits
    ax.set_ylim(-150,150)
    ax.set_xlim(-150,150)
    
    ax.set_aspect('equal')
    
    fig.savefig(f"{galaxy}_Density_Plots/{galaxy}_{snap}.png")

galaxies = np.array(["MW","M31"])
IDs = find_close_encounter_IDs()

#for galaxy in galaxies:
#    for snap in IDs:
#        Density_plot(galaxy, int(snap))

    
    
def Phase_plot(galaxy, snap):
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
    fig = plt.figure(figsize=(11,10))
    ax = plt.subplot(111)
    
    # Plot 2D Histogram for one component of  Pos vs Vel 
    # ADD HERE
    ax.hist2d(rn[:,0],vn[:,1], bins=800, norm=LogNorm(), cmap='plasma')
    # Hook to implement Mass Profile circular velocity 
    '''
    # Mass Profile
    MP = MassProfile(galaxy,snap)
    rr = np.arange(0.01, 100.1, 0.1)*u.kpc
    Vcirc = MP.circularVelocityTotal(rr)
    
    # Overplot Circular Velocity from the MassProfile Code
    plt.plot(rr, Vcirc, color='k')
    plt.plot(-rr,-Vcirc, c='k')
    '''
    # Add axis labels
    ax.set_xlabel('x [kpc]', fontsize=22)
    ax.set_ylabel('yv [km/s]', fontsize=22)
    ax.set_title(f'{galaxy} Face On Density; Snap={snap}', fontsize=28)
    
    #set axis limits
    ax.set_xlim(-150,150)
    
    fig.savefig(f"{galaxy}_Phase_Plots/{galaxy}_{snap}.png")

galaxies = np.array(["M31"])
for galaxy in galaxies:
    for snap in IDs:
        Phase_plot(galaxy, int(snap))

