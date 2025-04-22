
# Homework 5 Solutions
# Mass Profiles
# G. Besla & H. Foote


# import modules
import numpy as np
import astropy.units as u
import astropy.constants as const

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules from previous homeworks
from ReadFile import Read
from CenterOfMass2 import CenterOfMass



class MassProfile:
    '''Class that measures and plots mass profiles and rotation curves of
    simulation snapshots'''
    
    def __init__(self, galaxy, snap):
        ''' This class reads snapshots and plots the mass profiles 
        and rotation curves of galaxies.

        PARAMETERS
        ----------
        galaxy : `str; 'MW', 'M31', or 'M33'`
                Name of the galaxy to read in
        snap : `int`
            Number of the snapshot to read in
        '''
        
        # Determine Filename
        # add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        # create filenames
        self.filename=f'{galaxy}/{galaxy}_{ilbl}.txt'
        
        if galaxy != "M33":
            self.volDec = 2
        else:
            self.volDec = 4
        
        # read the particle data                                                                                             
        self.time, self.total, self.data = Read(self.filename)

        # store the mass, positions, velocities of all particles                                
        self.m = self.data['m']#*u.Msun
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
    
        # store galaxy name
        self.gname = galaxy
    
        # converting G to units of kpc*km^2/s^2/Msun
        self.G = const.G.to(u.kpc*u.km**2/u.s**2/u.Msun) 
    
    
    def massEnclosed(self, ptype, radii):
        '''This method computes and returns the mass profile of the galaxy
        based on the specified particle type.

        PARAMETERS
        ----------
        ptype : `int; 1, 2, or 3`
            particle type
        radii : `np.ndarray`
            array of radius bin edges, in kpc

        RETURNS
        -------
        m_enc : `np.ndarray`
            array containing the mass within the radii specified 
            by r, in Msun
        '''
    
        # Determine the COM position using Disk Particles
        # Disk Particles afford the best centroiding.
        # Create a COM object
        com = CenterOfMass(self.filename,2)
        # Store the COM position of the galaxy
        # Set Delta = whatever you determined to be a 
        #good value in Homework 4.
        com_pos = com.COM_P(self.volDec,0.1)
            
        # create an array to store indexes of particles of desired Ptype                                                
        index = np.where(self.data['type'] == ptype)

        # Store positions of particles of given ptype from the COMP. 
        xG = self.x[index] - com_pos[0]
        yG = self.y[index] - com_pos[1]
        zG = self.z[index] - com_pos[2]
            
        # Compute the mag. of the 3D radius
        rG = np.sqrt(xG**2 + yG**2 + zG**2)
            
        # store mass of particles of a given ptype
        mG = self.m[index]
            
        # Array to store enclosed mass as a function of the 
        #input radius array
        m_enc = np.zeros(np.size(radii))
        # equivalently: 
        # m_enc = np.zeros_like(radii)
    
        # loop through the radii array
        for i in range(np.size(radii)):
            # Only want particles within the given radius
            indexR = np.where(rG <  radii[i]*u.kpc)
            m_enc[i] = np.sum(mG[indexR])         
        
        # return the array of enclosed mass with appropriate units
        return m_enc*u.Msun*1e10
        
    
    def massEnclosedTotal(self, radii):    
        '''This method computes and returns the mass profile of 
        the galaxy based on ALL particles.

        PARAMETERS
        ----------
        radii : `np.ndarray`
            array of radius bin edges, in kpc

        RETURNS
        -------
        m_enc : `np.ndarray`
            array containing the mass within the radii
            specified by r, in Msun
        '''
     
        # Sum up all the mass of each component.
        m_enc = self.massEnclosed(1,radii) + self.massEnclosed(2,radii) + self.massEnclosed(3,radii)
    
        # Recall that M33 only has 2 components!  No bulge
        if (self.gname == 'M33'):
            m_enc = self.massEnclosed(1,radii)+ self.massEnclosed(2,radii)  
          
        return m_enc
    
        
        
    def hernquistMass(self, r, a, mhalo):
        ''' This method returns the mass enclosed within a radius based on
        the analytic Hernquist density profile.

        PARAMETERS
        ----------
        r : `float` 
            radius to compute mass within in kpc
        a : `float`
            Hernquist profile scale radius in kpc
        mhalo : `astropy.Quantity`
            total halo mass in Msun

        RETURNS
        -------
        m_enc : `astropy.Quantity'
            mass enclosed by r in Msun
        '''

        # adding units
        r = r*u.kpc
        a = a*u.kpc
        
        # compute numerator and denominator separately
        A = mhalo * r**2
        B = (a + r)**2
        
        return A/B
       
        
        
        
    def circularVelocity(self, ptype, radii):
        '''This method computes and returns the rotation curve of the galaxy
        based on the specified particle type.
    
        PARAMETERS
        ----------
        ptype : `int; 1, 2, or 3`
            particle type
        radii : `np.ndarray`
            array of radius bin edges, in kpc

        RETURNS
        -------
        v_circ : `np.ndarray`
            array containing the circular orbital velocity at
            the radii specified by r, in km/s
        '''
    
        # compute the mass enclosed for the given ptype
        m_enc = self.massEnclosed(ptype,radii)
        
        # Determine the circular speed as a function of 
        # input radius assuming spherical symmetry 
        # note that radius needs units 
        # This will return units of kpc/Gyr
        v_circ = np.sqrt(self.G*m_enc/(radii*u.kpc))
        
        # round and return the array in km/s
        return np.around(v_circ.to(u.km/u.s), 2)
    
    
    def circularVelocityTotal(self, radii):
        '''This method computes and returns the rotation curve 
        of the galaxy based on ALL particles.

        PARAMETERS
        ----------
        radii : `np.ndarray`
            array of radius bin edges, in kpc

        RETURNS
        -------
        v_circ : `np.ndarray`
            array containing the circular orbital velocity at
            the radii specified by r, in km/s

        '''
    
        # compute the total mass enclosed
        m_enc = self.massEnclosedTotal(radii)
        
        # Determine the circular speed 
        # note that radii need units . 
        #This will return units of kpc/Gyr 
        v_circ = np.sqrt(self.G*m_enc/(radii*u.kpc))
        
        # round and return the circular speed in km/s
        return np.around(v_circ.to(u.km/u.s), 2)
    
    
    
    def hernquistVCirc(self, radii, a, mhalo):
        '''This method returns the mass enclosed within a radius based on
        the analytic Hernquist density profile.

        PARAMETERS
        ----------
        r : `float` 
            radius to compute mass within in kpc
        a : `float`
            Hernquist profile scale radius in kpc
        Mhalo : `astropy.Quantity`
            total halo mass in Msun

        RETURNS
        -------
        v_circ : `np.ndarray'
            circular orbital velocity at r in km/s
        '''
           
        # Store the enclosed mass 
        m_enc = self.hernquistMass(radii,a,mhalo)
    
        # Circular speed using enclosed mass 
        # Note: You could also write this analytically without 
        #first calling HernquistMass
        v_circ = np.sqrt(self.G*m_enc/(radii*u.kpc))
        
        # Return circular speed, rounded and in km/s
        return np.around(v_circ.to(u.km/u.s), 2)
    