#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 20:51:55 2025

@author: petershea
"""

# Importing Libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import astropy.units as u
from astropy.constants import G
import scipy.optimize as sco

# Importing functions
from ReadFile import Read
from CenterOfMass import CenterOfMass

# defines G in the correct units
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

class MassProfile:
    
    def __init__(self, galaxy, snap):
        '''Class to derive mass profile and rotation curves for a given galaxy
        at a particular snapshot. Initializes global variables for the class
        
        Inputs:
            galaxy: string; the name of the desired galaxy (MW, M31, or M33)
            
        '''
        
        # add a string of the filenumber to the value “000”
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename="%s_"%(galaxy) + ilbl + '.txt'
        
        # Stores galaxy name as gname
        self.gname = galaxy
        
        # Reads in data from txt and stores in three variables        
        self.time, self.total, self.data = Read(self.filename)
        
        # Stores x, y, z coordinate arrays from data and assigns units 
        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc
    
    def MassEnclosed(self, ptype, radii):
        '''This function determines the mass of all particles of a particular 
        type within a given radius or array of radii
        
        Inputs:
            ptype: int; 1 - darkmatter, 2 - disk stars, 3 - bulge stars
            radii: astropy quantities [kpc]; an array of radii from the COM >0 
        
        Returns:
            masses: astropy quantities [Msun]; the array of mass values 
                containted within the given aray of radii
        '''
        
        # Creates Center of mass object for the given filename using disk stars
        COM = CenterOfMass(self.filename, 2)
        # Calculates position of the COM with error tolerance 0.1 kpc
        COM_P = COM.COM_P()
        
        # Indexes particles of a given type
        index = np.where(self.data['type'] == ptype)
        # Redefines the position of particles relative to the COM
        x_new = self.x[index] - COM_P[0]
        y_new = self.y[index] - COM_P[1]
        z_new = self.z[index] - COM_P[2]
        # Determines a distance from the COM for the particles
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)
        
        # Creates an empty array to hold mass values with the same length as 
        # the array of given radii
        masses = np.zeros(len(radii))
        
        # Iterates through radii array
        for i in range(len(radii)):
            # Stores radius as a variable
            radius = radii[i]
            # indexes particles within the given radius
            index_r = np.where(r_new < radius)
            # Uses the indexed particles to determine mass array
            m = self.data['m'][index_r]
            # Sums mass array for total mass of given ptype within radius and 
            # appends to the mass array
            masses[i] = np.sum(m)
        
        # applies proper units to the mass aray
        masses *= (1e10*u.M_sun)
        
        # Returns the mass array
        return masses
            
    def MassEnclosedTotal(self, radii):
        '''Function to determine the total mass enclosed of all particle types 
        at given radii
        
        Inputs: 
            radii: astropy quantities [kpc]; an array of radii from the COM >0 
        
        Returns:
            M_tot: astropy quantity [M_sun]; the total mass enclosed within the 
                radius
        '''
        
        # Calculates total mass for MW and M31
        if self.gname != 'M33':
            
            # Calls MassEnclosed function for each particle type and stores
            M_Halo = self.MassEnclosed(1, radii)
            M_Disk = self.MassEnclosed(2, radii)
            M_Bulge = self.MassEnclosed(3, radii)
            
            # Sums mass of all particle types
            M_tot = M_Halo + M_Disk + M_Bulge
        
        # Calculates total mass for M33 as it has no bulge
        else:
            
            # Calls MassEnclosed function for dark matter and disk stars
            M_Halo = self.MassEnclosed(1, radii)
            M_Disk = self.MassEnclosed(2, radii)
            
            # Sums mass of given prticle types
            M_tot = M_Halo + M_Disk
        
        # Returns total mass
        return M_tot
    
    def HernquistMass(self, r, a, MHalo):
        '''Function to model enclosed mass within a given radius using the 
        hernquist function 
        
        Inputs: 
            r: astropy quantity [kpc]; radius from the COM
            a: astropy quantity [kpc]; scale height of the galaxy
            MHalo: astropy quantity [Msun]; the total mass of the dark matter 
            halo of the galaxy
            
        Returns:
            HMass: astropy quantity [Msun]; Modeled enclosed masswithin radius 
                of COM
        '''
        
        # calculates the numerator of the function
        numerator = MHalo*r**2
        # calculates the denominator of the function
        denominator = (a + r)**2
        
        Hmass = numerator/denominator
        # divides and returns the value
        return Hmass
    
    def CircularVelocity(self, ptype, radii):
        '''Function to determine the circular velocity of a given particle type
        within a particular radius 
        
        Inputs:
            ptype: int; 1 - darkmatter, 2 - disk stars, 3 - bulge stars
            radii: astropy quantities [kpc]; an array of radii from the COM >0 
            
        Returns:
            v_c: astropy quantities [km/s]; The speed necessary to travel in a 
                circular orbit at the given gadii
        '''
        
        # uses the function v_c = sqrt(G*M / r) M is found by calling 
        # MassEnclosed function
        v_c = np.sqrt(G * self.MassEnclosed(ptype, radii) / radii)
        
        # Returns the circular velocity array rounded to 2 decimals
        return np.round(v_c,2)
    
    def CircularVelocityTotal(self, radii):
        '''Function to determine the Cicular velocity at a given radius due to 
        enclosed total mass 
        
        Inputs:
            radii: astropy quantities [kpc]; an array of radii from the COM >0 
        
        Returns:
            CV_tot: astropy quantity [km/s]; The speed necessary to travel in a
                circular orbit due to total enclosed mass
        '''
        
        # Calculates CV_tot for MW and M31
        if self.gname != 'M33':
            
            # Calls CircularVelociity function to determine v_c for each 
            # particle type
            CV_Halo = self.CircularVelocity(1, radii)
            CV_Disk = self.CircularVelocity(2, radii)
            CV_Bulge = self.CircularVelocity(3, radii)
            
            # Takes the magnitude of these to get total circular velocity
            CV_tot = np.sqrt(CV_Halo**2 + CV_Disk**2 + CV_Bulge**2)
        
        # calculates CV_tot for M33 as it has no bulge
        else:
            # Calls CircularVelociity function to determine v_c for each 
            # particle type
            CV_Halo = self.CircularVelocity(1, radii)
            CV_Disk = self.CircularVelocity(2, radii)
            
            # Takes the magnitude of these to get total circular velocity
            CV_tot = np.sqrt(CV_Halo**2 + CV_Disk**2)
        
        # Returns total circular velocity
        return CV_tot
    
    def HernquistVCirc(self, radii, a, MHalo):
        '''Function to model total circular velocity using mass enclosed by 
        hernquist function
        
        Inputs:
            r: astropy quantity [kpc]; radius from the COM
            a: astropy quantity [kpc]; scale height of the galaxy
            MHalo: astropy quantity [Msun]; the total mass of the dark matter 
            halo of the galaxy
        
        Returns:
            V_c: astropy quantity [km/s]; the modeled circular velocity due to 
            hernquist mass distribution
        '''
        
        # Calculates Hernquist mass array by calling function
        HMass = self.HernquistMass(radii, a, MHalo)
        
        # Calculates circular velocity using typical equation
        V_c = np.sqrt(G * HMass / radii)
        
        return V_c
        
############# Producing Mass Distributions and Rotation Curves #############

# Creates MassProfile objects for each galaxy at snap 000
MW = MassProfile('MW', 0)
M31 = MassProfile('M31', 0)
M33 = MassProfile('M33', 0)

# creates an array of radii from 0.1 to 30 kpc with 100 intervals
r = np.linspace(0.1, 30, 100) * u.kpc

# initalizes a plot of Mass profiles
fig1, ax1 = plt.subplots(1,3,figsize=(9,3.5), sharey=True, layout='constrained', dpi=1000)

# sets x axis limits for each subplot
for ax in ax1:
    ax.set_xlim(0,30)

# sets title and axis lables for subplots
fig1.suptitle('Galactic Mass Profiles', fontsize=14, fontweight='bold')
fig1.supylabel(r'$\log(M/M_{\odot})$', fontsize=10)
fig1.supxlabel('Radius from COM [kpc]', fontsize=10)

# Plots total and component mass enclosed against radius for Milky Way
ax1[0].semilogy(r,MW.MassEnclosedTotal(r),c='k') # total
ax1[0].semilogy(r,MW.MassEnclosed(1, r),linestyle='--', c='k') # dark matter
ax1[0].semilogy(r,MW.MassEnclosed(2, r),linestyle=':', c='k') # disk stars
ax1[0].semilogy(r,MW.MassEnclosed(3, r),linestyle='-.', c='k') # bulge stars
ax1[0].set_title('Milky Way') # sets title of subplot

# Plots total and component mass enclosed against radius for M31
ax1[1].semilogy(r,M31.MassEnclosedTotal(r),c='k') # total
ax1[1].semilogy(r,M31.MassEnclosed(1, r),linestyle='--', c='k') # dark matter
ax1[1].semilogy(r,M31.MassEnclosed(2, r),linestyle=':', c='k') # disk stars
ax1[1].semilogy(r,M31.MassEnclosed(3, r),linestyle='-.', c='k') # bulge stars
ax1[1].set_title('Andromeda (M31)') # sets title of subplot

# Plots total and component mass enclosed against radius for M33
ax1[2].semilogy(r,M33.MassEnclosedTotal(r),c='k') # total
ax1[2].semilogy(r,M33.MassEnclosed(1, r),linestyle='--', c='k') # dark matter
ax1[2].semilogy(r,M33.MassEnclosed(2, r),linestyle=':', c='k') # disk stars
ax1[2].set_title('Triangulum (M33)') # sets title of subplot

# defines variables for total Halo mass of each galaxy as found in HW 3
Halo_MW = 1.975e12
Halo_M31 = 1.921e12
Halo_M33 = 0.187e12

# Sets boundaries and initial guesses for scale height and total Halo mass for each galaxy
    # bounds were set up so total halo mass barely changes
p0MW = [0.5,Halo_MW]
boundsMW = ([1e-5,Halo_MW-1],[np.inf,Halo_MW+1])
p0M31 = [0.5,Halo_M31]
boundsM31 = ([1e-5,Halo_M31-1],[np.inf,Halo_M31+1])
p0M33 = [0.5,Halo_M33]
boundsM33 = ([1e-5,Halo_M33-1],[np.inf,Halo_M33+1])

# creates total enclosed mass arrays for each galaxy
MWmasses = MW.MassEnclosedTotal(r).value  
M31masses = M31.MassEnclosedTotal(r).value  
M33masses = M33.MassEnclosedTotal(r).value  

# Fit the Hernquist model to total enclosed mass arrays
poptMW, pcovMW = sco.curve_fit(MW.HernquistMass, r.value, MWmasses, p0=p0MW, bounds=boundsMW)
poptM31, pcovM31 = sco.curve_fit(MW.HernquistMass, r.value, M31masses, p0=p0M31, bounds=boundsM31)
poptM33, pcovM33 = sco.curve_fit(MW.HernquistMass, r.value, M33masses, p0=p0M33, bounds=boundsM33)

# Print the fitted parameters
print(f"Fitted parameter (a_MW): {poptMW[0]:.2e}")
print(f"Uncertainties in fitted parameter: {np.sqrt(np.diag(pcovMW)[0]):.2e}")
print(f"\nFitted parameter (a_M31): {poptM31[0]:.2e}")
print(f"Uncertainties in fitted parameter: {np.sqrt(np.diag(pcovM31)[0]):.2e}")
print(f"\nFitted parameter (a_M33): {poptM33[0]:.2e}")
print(f"Uncertainties in fitted parameter: {np.sqrt(np.diag(pcovM33)[0]):.2e}")

# Plot the Hernquist mass with the fitted parameters
ax1[0].semilogy(r, MW.HernquistMass(r.value, *poptMW), label="Hernquist Fit", linestyle='-', c='r')
ax1[1].semilogy(r, MW.HernquistMass(r.value, *poptM31), label="Hernquist Fit", linestyle='-', c='r')
ax1[2].semilogy(r, MW.HernquistMass(r.value, *poptM33), label="Hernquist Fit", linestyle='-', c='r')

# Add text for best fit Hernquist scale heighs
ax1[0].text(0.03, 0.99, r"$a_{\text{MW}}"+ f" = {poptMW[0]*u.kpc:.3f}$", transform=ax1[0].transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', color='k')
ax1[1].text(0.03, 0.99, r"$a_{\text{M31}}"+ f" = {poptM31[0]*u.kpc:.3f}$", transform=ax1[1].transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', color='k')
ax1[2].text(0.03, 0.99, r"$a_{\text{M33}}"+ f" = {poptM33[0]*u.kpc:.3f}$", transform=ax1[2].transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', color='k')

# Creates handles for linestyle legend
handles = [
    Line2D([0], [0], color='k', linestyle='-', label='Total Mass'),
    Line2D([0], [0], color='k', linestyle='--', label='Dark Matter'),
    Line2D([0], [0], color='k', linestyle=':', label='Disk Stars'),
    Line2D([0], [0], color='k', linestyle='-.', label='Bulge Stars'),
    Line2D([0], [0], color='r', linestyle='-', label='Hernquist Fit')
    ]

# Add the legend to the entire figure
ax1[2].legend(handles=handles, loc='lower right', bbox_to_anchor=(1.0, 0.0), fontsize=10)

fig1.savefig('Mass_Profiles.png')


# initializes a plot for roation curves for each galaxy
fig2, ax2 = plt.subplots(1,3,figsize=(9,3.5), sharey=True, layout='constrained', dpi=1000)

# sets x and y axis limits for each subplot
for ax in ax2:
    ax.set_xlim(0,30)
    ax.set_ylim(0,2000)

# sets title and axis lables
fig2.suptitle('Galactic Rotation Curves', fontsize=14, fontweight='bold')
fig2.supylabel(r'$\text{Circular Velocity} \, [\mathrm{km/s}]$', fontsize=10)
fig2.supxlabel('Radius from COM [kpc]', fontsize=10)

# plots total circular velocity and component velocities for Milky Way
ax2[0].plot(r,MW.CircularVelocityTotal(r),c='k') # total
ax2[0].plot(r,MW.CircularVelocity(1, r),linestyle='--', c='k') # dark matter
ax2[0].plot(r,MW.CircularVelocity(2, r),linestyle=':', c='k') # disk stars
ax2[0].plot(r,MW.CircularVelocity(3, r),linestyle='-.', c='k') # bulge stars
ax2[0].set_title('Milky Way') # titles subplot

# plots total circular velocity and component velocities for M31
ax2[1].plot(r,M31.CircularVelocityTotal(r),c='k') # total
ax2[1].plot(r,M31.CircularVelocity(1, r),linestyle='--', c='k') # dark matter
ax2[1].plot(r,M31.CircularVelocity(2, r),linestyle=':', c='k') # disk stars
ax2[1].plot(r,M31.CircularVelocity(3, r),linestyle='-.', c='k') # bulge stars
ax2[1].set_title('Andromeda (M31)') # titles subplot

# plots total circular velocity and component velocities for M33
ax2[2].plot(r,M33.CircularVelocityTotal(r),c='k') # total
ax2[2].plot(r,M33.CircularVelocity(1, r),linestyle='--', c='k') # dark matter
ax2[2].plot(r,M33.CircularVelocity(2, r),linestyle=':', c='k') # disk stars
ax2[2].set_title('Triangulum (M33)') # titles subplot

# plots modeled hernquist circular velocities using the fit scale heights found previously
ax2[0].plot(r, MW.HernquistVCirc(r.value, *poptMW), label="Hernquist Fit", linestyle='-', c='r')
ax2[1].plot(r, MW.HernquistVCirc(r.value, *poptM31), label="Hernquist Fit", linestyle='-', c='r')
ax2[2].plot(r, MW.HernquistVCirc(r.value, *poptM33), label="Hernquist Fit", linestyle='-', c='r')

# creates the legend
ax2[2].legend(handles=handles, loc='upper right', bbox_to_anchor=(1, 1), fontsize=10)

# adds scale height used to the top left corner of each subplot
ax2[0].text(0.1, 0.99, r"$a_{\text{MW}}"+ f" = {poptMW[0]*u.kpc:.3f}$", transform=ax2[0].transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', color='k')
ax2[1].text(0.03, 0.99, r"$a_{\text{M31}}"+ f" = {poptM31[0]*u.kpc:.3f}$", transform=ax2[1].transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', color='k')
ax2[2].text(0.03, 0.99, r"$a_{\text{M33}}"+ f" = {poptM33[0]*u.kpc:.3f}$", transform=ax2[2].transAxes, fontsize=8, verticalalignment='top', horizontalalignment='left', color='k')

# saves the figure
fig2.savefig('Rotation_Curves.png')
            