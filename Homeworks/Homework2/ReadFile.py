#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 21:15:23 2025

@author: petershea
"""

import astropy.units as u
import numpy as np

def Read(filename):
    '''
    This function reads in data from the given file
    Inputs:
        filename: string,in the form 'MW_***.txt'
    Returns:
        time: time the simulation step occurs at in Myr
        n_particles: the number of particles in the simulation at this step
        data: the remaining data from teh file organized as a table 
    '''
    
    #opens the called file in read mode
    f = open(filename,'r')
    
    #reads the first line, splits contents, and saves time with correct units
    line1 = f.readline()
    label, value = line1.strip().split()
    time = float(value) * u.Myr
    
    #reads the second line, splits contents, and saves number of particles
    line2 = f.readline()
    label, value = line2.strip().split()
    n_particles = float(value)
    
    #closes file
    f.close()
    
    #defines a table for the rest of the data skipping the header
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    #returns, time step in Myrs, number of particles, and data table
    return time, n_particles, data
