#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 12:24:28 2025

@author: petershea
"""
# Imports required functions and libraries
import numpy as np
import astropy.units as u
from astropy.table import Table
from astropy.io import ascii
from ReadFile import Read

def ComponentMass(filename, particle_type):
    '''Returns the total mass of the designated component of a given galaxy
    
    Inputs:
        filename: str, filename of the given galaxy ex. MW_000.txt
        particle_type: integer;1= Dark matter, 2= Disk Star, 3= Halo/Bulge Star
    
    returns:
        M_tot: astropy quantity [10^12 M_sun], the total mass of the given 
            component
    '''
    
    # Reads in the data contained in the given file
    data = Read(filename)[2]
    # Reads in the number of particles contained in the given file
    n_particles = Read(filename)[1]
    
    # Determines the frist index of the desired particle type, if there are no 
    # particles of that type zero is returned
    try:
        key_index = np.where(data['type'] == particle_type)[0][0]
    except:
        return 0 * u.M_sun * 10**12
    
    # Defines a mutable index to navigate data
    index = key_index
    
    # Defines a variable to contain total mass of the given component
    M_tot = 0
    
    # Navigates through data table using the index from the first instance of a 
    # particle type through the last. Might cause issues in the future if data 
    # isn't always organized by particle type. Can fix by removing particle 
    # type check from while and creating an if statement. Alternatively, could 
    # limit data to only particles of given type then navigate through.
    while index < n_particles and data['type'][index] == particle_type:
        M_tot += data['m'][index] # Adds particle mass to mass total
        index += 1 # advances the index
    
    # Applies astropy unit 10^12 solar masses
    M_tot = M_tot * (u.M_sun * 10**10)
    
    # Defines new unit of specifically 10^12 solar masses
    TM_sun = u.M_sun * 1e12
    
    # coverts from Solar masses to 10^12 Solar masses as the unit
    M_tot = M_tot.to(TM_sun)
    
    # Rounds mass total to 3 decimal places
    M_tot = np.round(M_tot,3)
    
    # Returns total component mass
    return M_tot

# Uses ComponentMass function to determine halo, disk, and bulge mass of the 
# Milky Way. Only value and not astropy units was used in the table
MW_M_Halo = ComponentMass('MW_000.txt', 1).value
MW_M_Disk = ComponentMass('MW_000.txt', 2).value
MW_M_Bulge = ComponentMass('MW_000.txt', 3).value

# Uses ComponentMass function to determine halo, disk, and bulge mass of the 
# M31. Only value and not astropy units was used in the table
M31_M_Halo = ComponentMass('M31_000.txt', 1).value
M31_M_Disk = ComponentMass('M31_000.txt', 2).value
M31_M_Bulge = ComponentMass('M31_000.txt', 3).value

# Uses ComponentMass function to determine halo, disk, and bulge mass of the 
# M33. Only value and not astropy units was used in the table
M33_M_Halo = ComponentMass('M33_000.txt', 1).value
M33_M_Disk = ComponentMass('M33_000.txt', 2).value
M33_M_Bulge = ComponentMass('M33_000.txt', 3).value

# Calculates the total hass of the halo, disk, and bulge of the three major 
# objects of the local group
Local_M_Halo = MW_M_Halo + M31_M_Halo + M33_M_Halo
Local_M_Disk = MW_M_Disk + M31_M_Disk + M33_M_Disk
Local_M_Bulge = MW_M_Bulge + M31_M_Bulge + M33_M_Bulge

# Calculates the total mass of the Milky Way, M31, M33, and the local group by 
# adding Halo, Disk, and Bulge masses
MW_M_tot = MW_M_Halo + MW_M_Disk + MW_M_Bulge
M31_M_tot = M31_M_Halo + M31_M_Disk + M31_M_Bulge
M33_M_tot = M33_M_Halo + M33_M_Disk + M33_M_Bulge
Local_M_tot = MW_M_tot + M31_M_tot + M33_M_tot

# Calculates the baryon fraction of the Milky Way, M31, M33, and the local 
# group
MW_fbar = np.round((MW_M_Disk + MW_M_Bulge) / MW_M_tot,4)
M31_fbar = np.round((M31_M_Disk + M31_M_Bulge) / M31_M_tot,4)
M33_fbar = np.round((M33_M_Disk + M33_M_Bulge) / M33_M_tot,4)
Local_fbar = np.round((Local_M_Disk + Local_M_Bulge) / Local_M_tot,4)

# arrays containing the data for each column
Galaxy_names = ['Milky Way','Andromeda (M31)','Triangulum (M33)', # Names
                r'\hline Local Group']
Halo_mass = [MW_M_Halo, M31_M_Halo, M33_M_Halo, Local_M_Halo] # Halo Mass
Disk_mass = [MW_M_Disk, M31_M_Disk, M33_M_Disk, Local_M_Disk] # Disk Mass
Bulge_mass = [MW_M_Bulge, M31_M_Bulge, M33_M_Bulge, Local_M_Bulge] # Bulge Mass
Total_mass = [MW_M_tot, M31_M_tot, M33_M_tot, Local_M_tot] # Total Mass
Baryon_fraction = [MW_fbar, M31_fbar, M33_fbar, Local_fbar] # Baryon Fraction


# Array of strings that will become the titles of the columns, name [units]
Column_names = ['Galaxy Name',
                r'Halo Mass [$10^{12} M_{\odot}$]',
                r'Disk Mass [$10^{12} M_{\odot}$]',
                r'Bulge Mass [$10^{12} M_{\odot}$]',
                r'Total Mass [$10^{12} M_{\odot}$]',
                'Baryon Fraction',]

# Creates the 2D array of data using the arrays containing mass components
Data = [Galaxy_names,
        Halo_mass,
        Disk_mass,
        Bulge_mass,
        Total_mass,
        Baryon_fraction]

# The Dtype associated with the data in each column of the table, not necessary
# if the information in each column has the same dtype
d_types = ['U20',
           'f4',
           'f4',
           'f4',
           'f4',
           'f4']

def Table_Write(Data, Column_names, d_types):
    '''Creates an astropy table and returns the coresponding LaTex code for it
    
    Inputs:
        Data:  2D array, contains the information that will be turned into the 
            table, each subarray should have only one dtype
        Column_names: array, a list of strings that will become the names of 
            the columns of the table, len(Data)
        d_types: array, a list of strings corresponding to the dtype of the 
            information stored in each subarray (U20: String, f4: float, 
                                                 i4: integer)
    Returns:
        t: astropy table created 
        
        Writes LaTex code for the table of mass components, needs to be edited 
        to ensure style and format is good and LaTex document is created
    '''
    
    # Creates an astropy table using the three inputs
    t = Table(Data,names=Column_names,dtype=d_types)
    
    # Determines style of the LaTex table created, AA is Astropysics style, 
    # col_align is how the data is aligned in the table l-left, c-center. This 
    # part should be tailored for each table individually
    latex_dict = ascii.latex.latexdicts['AA']
    latex_dict.update({'col_align': 'lccccc'})
    
    # Writes the LaTex code used to create the table to the terminal
    ascii.write(t,
                Writer=ascii.Latex, # Determines output is in LaTex
                latexdict=latex_dict, # Uses latex_dict to define table style
                overwrite=True) # Allows the script to be edited and run again 
    
    # Returns the Table
    return(t)

# Initializes the Table_Write function, uncomment to do so
Table_Write(Data, Column_names, d_types)
    
# Output was edited to better control style and format
    