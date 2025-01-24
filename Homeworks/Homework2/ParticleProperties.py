#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 22:17:59 2025

@author: petershea
"""

import numpy as np
import astropy.units as u
from ReadFile import Read


def ParticleInfo(filename, particle_type, particle_number):   
    '''
    This function provides the 3D Distance, 3D velocity magnitude, and mass of 
    a given particle given its  number, type, and the file.
    Inputs:
        filename: string,in the form 'MW_***.txt'
        particle_type: integer;1= Dark matter, 2= Disk Star, 3= Halo/Bulge Star
        particle_number: integer, the position of the particle selected
    
    Returns:
        Distance: index 0, distance of the particle from origin in kpc
        Velocity Magnitude: index 1, 3D velocity of the particle in km/s
        Mass: index 2, mass of the selected particle in solar masses
    '''
    
    #Reads in the data from the txt file using the read function
    data = Read(filename)[2]
    
    #Determines the frist index of the desired particle type
    key_index = np.where(data['type'] == particle_type)[0][0]
    #translates both the particle type and particle number into proper index
    index = particle_number - 1 + key_index
    
    #defines x, y, and z values for the particle using the correct index
    x = data['x'][index]
    y = data['y'][index]
    z = data['z'][index]
    
    #calcuulates distance using typical formula and defined position values
    dist = np.sqrt(x**2 + y**2 + z**2) * u.kpc

    #defines vx, vy, and vz values for the particle using the correct index
    vx = data['vx'][index]
    vy = data['vy'][index]
    vz = data['vz'][index]
    
    #calculates velocity magnitude using defined velocity components
    vmag = np.sqrt(vx**2 + vy**2 + vz**2) * u.km / u.s
    
    #defines mass using the correct index
    m = data['m'][key_index] * u.M_sun
    
    #returns distance and velocity rounded to 3 decimals, as well as mass
    return np.round(dist,3), np.round(vmag,3), m


# Test Print Statement for Question 5
print('Distance:', ParticleInfo('MW_000.txt', 2, 100)[0],
      '\n         ', np.round(ParticleInfo('MW_000.txt',2,100)[0].to(u.lyr),3),
      '\nVelocity Magnitude:', ParticleInfo('MW_000.txt', 2, 100)[1],
      '\nMass:', ParticleInfo('MW_000.txt', 2, 100)[2])

