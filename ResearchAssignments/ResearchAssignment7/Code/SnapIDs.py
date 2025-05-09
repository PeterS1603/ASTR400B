#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 11:17:06 2025

@author: petershea
"""

import numpy as np
from OrbitCOM import readOrbit
from OrbitCOM import relative

########## Snap IDs of Interest ##########
def find_close_encounter_IDs():
    ''' Function that determines the snap ids of each close encounter and 
    several after. This will help determine which HiRes files to sftp.
    
    Parameters:
        None.
        
    Returns:
        snap_ids: ndarray; snap ids of each close encounter and 5 and 10 ids 
            after. Also includes the final snap id to see end of simulation 
            state
    '''
    
    MW_data = readOrbit('Orbits/Orbit_MW.txt')
    M31_data = readOrbit('Orbits/Orbit_M31.txt')
    
    #time = MW_data[0]
    
    # of MW and M31
    MW_M31_relPos = relative(MW_data[1],M31_data[1])
    
    delta = 0
    close_encounters = np.array([])
    
    min_index = np.argmin(MW_M31_relPos)
    
    for i in range(len(MW_M31_relPos)-1):
        delta_new = MW_M31_relPos[i+1]-MW_M31_relPos[i]
        
        if delta_new * delta < 0 and i-1 < min_index:
            close_encounters = np.append(close_encounters,i)
        
        delta = delta_new

    
    
    for i in close_encounters:
        close_encounters = np.append(close_encounters,i+1)
        close_encounters = np.append(close_encounters,i+2)
    
    close_encounters = np.sort(close_encounters)
    
    snap_ids = close_encounters*5
    
    snap_ids = np.append(snap_ids,801.)
    
    return(snap_ids)

snap_ids = find_close_encounter_IDs()

def sftp_command(galaxy, snap_ids):
    f = open('sftp_snap_get.txt','w')
    
    f.write(f'cd ../astr400b/HighRes/{galaxy}')
    
    for i in snap_ids:
        snap = int(i)
        f.write(f'\nget {galaxy}_{snap}.txt')
    
    f.write('\nquit')

    f.close()
    
#sftp_command('M33', snap_ids)
    






