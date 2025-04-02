
# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 




# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass2 import CenterOfMass


# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

from OrbitCOM import readOrbit

# # M33AnalyticOrbit


class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): # **** add inputs
        """
        Parameters:
            filename: str; a filename to save the analytic orbit of M33
            
        Returns:
            None
        """

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33COM = CenterOfMass("M33_000.txt",2)
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        M33COM_pos = M33COM.COM_P(4,0.1).value
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33COM_vel = M33COM.COM_V(*M33COM.COM_P(4,0.1)).value
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31COM = CenterOfMass("M31_000.txt",2)
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        M31COM_pos = M31COM.COM_P(2,0.1).value
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31COM_vel = M31COM.COM_V(*M31COM.COM_P(2,0.1)).value
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = M33COM_pos - M31COM_pos
        self.v0 = M33COM_vel - M31COM_vel
        
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5

        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass("M31_000.txt", 2).value * 1e12
        
        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1

        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass("M31_000.txt", 3).value * 1e12
        
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62

        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = ComponentMass("M31_000.txt", 1).value * 1e12
    
    
    def HernquistAccel(self, M, r_a, r): # it is easiest if you take as an input the position VECTOR 
        """ This function calculates acceleration of the halo and bulge using a
        hernquist model of the halo
        
        Parameters:
            M: float [Msun]; mass of the galaxy component which the 
            acceleration is being calculated for 
            r_a: float [kpc]; scale height of the galaxy component
            r: ndarray [kpc]; position vector for the COM of M33 [x,y,z]
            
        Returns:
            Hern: ndarray [kpc/Gyr^2]; acceleration vector for the COM of M33
        """
        
        ### **** Store the magnitude of the position vector
        rmag = np.linalg.norm(r)
        
        ### *** Store the Acceleration
        numerator = - self.G * M * r
        denominator = rmag * ((r_a + rmag)**2)
        Hern =  numerator / denominator #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r):# it is easiest if you take as an input a position VECTOR  r 
        """
        This function calculates the acceleration of the disk using the 
        Miyamoto Nagai Model
        
        Parameters:
            M: float [Msun]; mass of the disk
            r_d: float [kpc]; scale height of the disk
            r: ndarray [kpc]; position vector for the COM of M33 [x,y,z]
        
        """

        z_d = r_d/5 # Defines z_d based on the disk scale height

        R = np.sqrt(r[0]**2 + r[1]**2) # defines R using x and y component of r
        B = r_d + np.sqrt(r[2]**2 + z_d**2) # defines B using given equation
        
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        
        numerator = - self.G * M # calculates numerator using globally defined G
        denominator = (R**2 + B**2)**1.5 # calculates denominator
        vectors = r * np.array([1,1,B/np.sqrt(r[2]**2 + z_d**2)]) # determines product of r and given vector
        
        MiyamotoProfile = numerator/denominator * vectors # simplified equation
       
        return MiyamotoProfile
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r): # input should include the position vector, r
        """
        This function determines the total acceleration due to M31 each 
        component
        
        Paramaters:
            r: ndarray [kpc]; position vector for the COM of M33 [x,y,z]
            
        Returns:
            a_tot: ndarray [kpc/Gyr^2]; total acceleration of the COM of M33
        
        """

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        Halo_a = self.HernquistAccel(self.Mhalo, self.rhalo, r)
        Bulge_a = self.HernquistAccel(self.Mbulge, self.rbulge, r)
        Disk_a = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r)
        
        a_tot = np.sum((Halo_a, Bulge_a, Disk_a), axis=0)
            
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return a_tot
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """
        This function uses the leapfrog method to calculate an analytical model 
        of M33's orbit
        
        Parameters:
            dt: float [Gyr]; the interval of time over which M33's position is 
                calculated
            r: ndarray [kpc]; position vector for the COM of M33 [x,y,z]
            v: ndarray [kpc/Gyr]; velocity vector for the COM of M33 [vx,vy.vz]
        
        Returns:
            rnew: ndarray [kpc]; the new position of M33's COM at the next time 
                interval of the analytical model
            vnew: ndarray [kpc/Gyr]; the new velocity of M33's COM at the next 
                time interval of the analytical model
        """
        
        # predict the position at the next half timestep
        rhalf = r + v*(dt/2)
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + self.M31Accel(rhalf)*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew*(dt/2)
        
        return rnew, vnew# **** return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """
        This function uses the leapfrog method to create the analytical model 
            with the given time interval 
            
        t0: float [Gyr]; The initial time of the analytical model
        dt: float [Gyr]; The size of interval used by the leapfrog method
        tmax: float [Gyr]; The end of the analytical model
        
        """

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros([int(tmax/dt)+2, 7])
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t <= tmax):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t += dt
           
            # **** store the new time in the first column of the ith row
            orbit[i][0] = t
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            orbit[i,1:4], orbit[i,4:7] = self.LeapFrog(dt, orbit[i-1,1:4], orbit[i-1,4:7])
         
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            
            
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            
            
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1
        
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

# defines the time interval of the analytical orbit model
dt = 0.1

# creates an instance of the analytical orbit class and creates a model from 0 
# to 10 gyrs
#M33AnalyticOrbit(f"M33AnalyticOrbit_{dt}.txt").OrbitIntegration(0, dt, 10)

# stores the position, velocity, and time information for the analytical model 
# and simulation
M33orbit_analytic = readOrbit(f"M33AnalyticOrbit_{dt}.txt")
M33orbit_sim_pos = readOrbit("Orbit_M33.txt")[1] - readOrbit("Orbit_M31.txt")[1]
M33orbit_sim_vel = readOrbit("Orbit_M33.txt")[2] - readOrbit("Orbit_M31.txt")[2]
sim_time = readOrbit("Orbit_M33.txt")[0]

# creates a figure to plot the position and velocity of the analytical model vs 
# simulation 
fig, ax = plt.subplots(1,2,figsize=(8,4.5), layout='constrained')

ax[0].plot(M33orbit_analytic[0],np.linalg.norm(M33orbit_analytic[1], axis=0), label='Analytic Model')
ax[0].plot(sim_time,np.linalg.norm(M33orbit_sim_pos, axis=0), label='Simulation Data')

ax[1].plot(M33orbit_analytic[0],np.linalg.norm(M33orbit_analytic[2], axis=0), label='Analytical Model')
ax[1].plot(sim_time,np.linalg.norm(M33orbit_sim_vel, axis=0), label='Simulation Data')

fig.suptitle("Analytic Model Vs Simulation Data of M33's Orbit", fontsize=16)
fig.supxlabel("Time [Gyrs]")

ax[0].set_ylabel("M33 COM Position Mognitude [kpc]")
ax[0].set_title("M33 - M31 Position")
ax[1].set_ylabel("M33 COM velocity Mognitude [kpc/Gyr]")
ax[1].set_title("M33 - M31 Velocity")



########## HW Questions ##########
print("########## Analysis ##########")
print("Question 2")
print('''\tWhile up until the 1st Gyr, the analytic model of M33's orbit 
      \tmatches quite closely with the simulaltion data, after that point the
      \ttwo models greatly differ from each other. In the 10 Gyrs for which
      \tthe analytic model ran, M33 did not complete an entire orbit around
      \tM31 and orbited to a much greater distance. In comparison, the 
      \tsimulation suggests M33's orbit to be far more chaotic, have a much 
      \tshorter period, and to be decaying indicating an eventual merger.''')
      
print('\nQuestion 3')
print('''\tThere are likely two significant factors which lead to this 
      \tsignificant difference in orbits between the simulation and the model.
      \tThese are: 1) the absense of the Milkyway merger from the analytic 
      \torbit model (neglecting the added mass from the milky way after the 
      \tmerger around ~6 Gyrs), and 2) the lack of dynamical friction in the 
      \tsystem (by modeling M33's orbit as a point mass at the COM, we neglect 
      \tmuch of the extended structure and the role it plays in the orbital 
      \tmechanics).''')

print('\nQuestion 4')
print('''\tIn order to include the Milky Way in these calculations we could add
      \tanother pointmass to the analytic model representing the MW's COM. 
      \tWhile this would neglect more dynamical friction, it would solve the 
      \tproblem of mass of the cetral body not increasing enough.''')


