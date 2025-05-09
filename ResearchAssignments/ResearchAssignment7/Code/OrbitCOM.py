

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
import scipy.signal
from scipy.optimize import curve_fit


# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass




def OrbitCOM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM 
    pos and vel as a function of time.
    inputs:
        galaxy: string; name of target: MW, M31, or M33
        start: int; the first snap id in the desired range, inclusive
        end: int; the upper bounds snap id in desired range, exclusive
        n: int; step size of snap ids
          
    outputs: 
        Creates txt file containing x, y, z, vx, vy, and vz for the given 
        galaxy at time interval t. File Called Orbit_{galaxy}.txt
    """
    
    # compose the filename for output
    fileout = f'Orbit_{galaxy}.txt'
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    volDec = 2
    volDec_M33 = 4
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end, n)
    
    if len(snap_ids) == 0:
        print('Invalid Array')
        return
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids),7])
    
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        ilbl = '000' + str(snap_id)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        # Creates relative path and filename
        filename=f"{galaxy}/{galaxy}_{ilbl}.txt"

        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename, 2)
        
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        # Uses volDec = 2 for MW and M31
        if galaxy != 'M33':
            COM_pos = COM.COM_P(volDec,delta)
            # Uses volDec_M33 = 4 for M33
        else:
            COM_pos = COM.COM_P(volDec_M33,delta)
        
        # Stores COM velocity
        COM_vel = COM.COM_V(*COM_pos)
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        # stores time converted to Gyr
        t = COM.time.to(u.Gyr)
        # stores time, x, y, z, vx, vy, vz values (not units)
        orbit[i] = t.value, *COM_pos.value, *COM_vel.value
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
"""
# Initializes important variables for OrbitCOM function
start = 0
end = 801
n = 5
galaxies = np.array(['MW','M31','M33'])
"""
# loops through galaxies array to run OrbitCOM for all 3 objects
#for galaxy in galaxies:
    #print(f'\t{galaxy}')
    #OrbitCOM(galaxy, start, end, n)


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt

# defines a new function to read in 
def readOrbit(filename):
    '''function to read in data from Orbit_{galaxy}.txt file created by 
    OrbitCOM function
    
    inputs: 
        filename: string; name of desired file "Orbit_{galaxy}.txt"
        
    Returns:
        t: astropy quantities [Gyr]; an array of time values for each step of 
            the simulation
        COMP: astropy quantities [kpc]; a 2d array of x, y, and z position
            components of the COM at each step of the simulation
        COMV: astropy quantities [km/s]; a 2d array of vx, vy, and vz velocity 
            components of the COM at each step of the simulation
    '''
    # Uses gennfromtxt to create data table from filename
    data = np.genfromtxt(filename, comments='#', names=True)
    
    # Defines the time array from the data table
    time = data['t']
    
    #  Defines the COM pos and COM vel arrays from the data 
    COMP = np.array([data['x'],data['y'],data['z']])
    COMV = np.array([data['vx'],data['vy'],data['vz']])
    
    return time, COMP, COMV
"""
# reads in COM orbit data for each galaxy using the readOrbit function
MW_data = readOrbit('Orbit_MW.txt')
M31_data = readOrbit('Orbit_M31.txt')
M33_data = readOrbit('Orbit_M33.txt')
"""
# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  

def relative(vector1, vector2):
    '''Function determines the difference between two vectors 
    
    inputs:
        vector1: n dimension dnarray
        vector1: n dimension dnarray
        
    '''
    # Determines the difference between the two vectors
    rel_vector = vector2 - vector1
    
    # Calculates the magnitude of the difference
    mag = np.linalg.norm(rel_vector, axis=0)
    
    return mag
"""
# Determine the magnitude of the relative position and velocities 
# defines a time array using an arbitary galaxy 
time = MW_data[0]

# of MW and M31
MW_M31_relPos = relative(MW_data[1],M31_data[1])
MW_M31_relVel = relative(MW_data[2],M31_data[2])

# of M33 and M31
M33_M31_relPos = relative(M31_data[1],M33_data[1])
M33_M31_relVel = relative(M31_data[2],M33_data[2])

# Plot the Orbit of the galaxies 
#################################

# Creates plot of relative position of MW and M31 in one subplot and M31 and 
# M33 in another
fig1, ax1 = plt.subplots(1,2,figsize=(6,3.75),layout='constrained',dpi=1000)

# Adds Titles and Axis lables to the figure
fig1.suptitle('Relative COM Position', fontsize=16)
fig1.supxlabel('Simulation Time [Gyr]')
fig1.supylabel('COM Separation [kpc]')

# plots data and asigns sublables to subplots
ax1[0].plot(time, MW_M31_relPos)
ax1[0].set_title('MW-M31')
ax1[1].plot(time,M33_M31_relPos)
ax1[1].set_title('M31-M33')

# Saves the figure
#fig1.savefig('Relative_Position.png')


# Plot the orbital velocities of the galaxies 
#################################

# Creates plot of relative velocities of MW and M31 in one subplot and M31 and 
# M33 in another
fig2, ax2 = plt.subplots(1,2,figsize=(6,3.75),layout='constrained', dpi=1000)

# Adds Titles and Axis lables to the figure

fig2.suptitle('Relative COM Velocity', fontsize=16)
fig2.supxlabel('Simulation Time [Gyr]')
fig2.supylabel('Velocity Magnitude [km/s]')

# plots data and asigns sublables to subplots
ax2[0].plot(time,MW_M31_relVel)
ax2[0].set_title('MW-M31')
ax2[1].plot(time,M33_M31_relVel)
ax2[1].set_title('M31-M33')

# Saves the figure
#fig2.savefig('Relative_Velocity.png')




# Answers to Question 4
########################

### Part 1 ###
print('Problem 1')
print('''\tThe Milkyway and Andromeda appear to have 2 close encounters before they 
\tmerge. This can be seen in 2 local minima in the relative position of the COM
\tof each object before an increase in the COM separation until finally there 
\tis no increase and the separation remains at or close to 0.
''')
### Part 2 ###
print('\nProblem 2')
print('''\tIt is clear to see that as the separation of objects decreases their relative
\tvelocity increases. Thus there is an inverse relatinship between separation
\tand relative velocity as time evolves which is expected as the objects would 
\taccelerate towards each other.
''')
### Part 3 ###
print('\nProblem 3')
# finds the index of the minimum separation of the two COMs which coincides with the merger
min_index = np.argmin(MW_M31_relPos)
print(f'\tTime of MW-M31 Merge: {time[min_index]}')
print('''\tWhen MW and M31 merge M33's orbit decays significantly shown by the 
\tsignificant decrease in the separation of of M31 and M33 between the maximum 
\tseparation before and after the merger.
''')
### Part 4 ###

print('\nProblem 4')

local_maxima = scipy.signal.argrelextrema(M33_M31_relPos, np.greater)[0] # finds local maxima

t_min = 6 # Gyrs

# the indexs of appocenters after 6 Gyrs 
local_maxima = local_maxima[time[local_maxima] > t_min]

# determines the separation values at apocenter
Apocenters = M33_M31_relPos[local_maxima]
# determines the time when these appocenters occur
Apo_time = time[local_maxima]

# Prints Apocenters after t_min
print(f"\tApocenters after {t_min} Gyrs:")
for xm, ym in zip(Apo_time, Apocenters):
    print(f"\t\tt = {xm:.2f} [Gry], Separation = {ym:.2f} [kpc]")

# Compute average orbitaal period
if len(Apo_time) > 1:
    # determines the period between each apocenter and averages
    orbital_period = np.mean(np.diff(Apo_time)) 
else:
    orbital_period = float('nan')  # Handle case with insufficient data

print(f'\tAvg Orbital Period: {orbital_period:.3f} [Gyr]')

# Compute average decay rate
if len(Apocenters) > 1:
    decay_rate = np.mean(Apocenters[1:] / Apocenters[:-1]) # divides all apocenters
else:
    decay_rate = float('nan')  # Handle case with insufficient data

print(f'\tAverage Decay Rate: {decay_rate:.3f}')

# lambda function for the exponential decay of the orbit
Orbit_decay = lambda t, decay, Orb_per: Apocenters[-1]* decay ** ((t - Apo_time[-1]) / Orb_per)
# applies reasonable bounds to parameters
bounds = ([0,np.min(np.diff(Apo_time))],[1,np.max(np.diff(Apo_time))])
# uses curve fit to determine appropriate orbital decay constant and orbital period
popt, perr = curve_fit(Orbit_decay, Apo_time, Apocenters, bounds=bounds, p0=[decay_rate, orbital_period])

# Extract best-fit decay rate and orbital period
decay_fit, Orb_per_fit = popt

print(f'\tFit Orbital Period: {Orb_per_fit:.3f} [Gyr]')
print(f'\tFit Decay Rate: {decay_fit:.3f}')

# plots found Apocenters and Orbital decay fit to ensure accuracy
#ax1[1].scatter(Apo_time,Apocenters)
#ax1[1].plot(time,Orbit_decay(time,*popt))

merg_thresh = 5 #kpc     Merge Threshold
t_init = Apo_time[-1] # Gyr
delta_t = 0.001 # Gyr

# while loop to determine time of merge starting from the last apocenter
t = t_init
while Orbit_decay(t,*popt) > merg_thresh:
    t += delta_t

print(f"\tM33 merges with MW-M31 remnant at t ≈ {t:.3f} Gyr")
"""




