
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
from CenterOfMass import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

# # M33AnalyticOrbit




class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): # **** add inputs
        """
        Class to calculate the analytical orbit of M33 around M31
        
        Parameters
        ----------
        filename : String
            The name of the output file that will be written to.
            
        Output
        ------
        None.
        """

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33COM = CenterOfMass("M33_000.txt", 2)

        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        M33COM_p = M33COM.COM_P(0.1)

        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33COM_v = M33COM.COM_V(M33COM_p[0], M33COM_p[1], M33COM_p[2]).value
        
        M33COM_p = M33COM_p.value
        
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31COM = CenterOfMass("M31_000.txt", 2)

        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        M31COM_p = M31COM.COM_P(0.1)

        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31COM_v = M31COM.COM_V(M31COM_p[0], M31COM_p[1], M31COM_p[2]).value
        
        M31COM_p = M31COM_p.value
        
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = M33COM_p - M31COM_p
        self.v0 = M33COM_v - M31COM_v
        
        
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5 # kpc

        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass("M31_000.txt", 2) * 1e12 # Msun
        
        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1 # kpc

        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass("M31_000.txt", 3) * 1e12 # Msun
        
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 25 # kpc

        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = ComponentMass("M31_000.txt", 1) * 1e12
    
    
    def HernquistAccel(self, M, r_a, r): # it is easiest if you take as an input the position VECTOR 
        """ 
        Function that calculate the acceleration from a Hernquist potential
        
        Parameters
        ----------
        M : float
            The total mass of the particle type in the galaxy
        r_a : float
            The Hernquist scale length
        r : np.ndarray of floats
            The position vector
            
        Output
        ------
        Hern : np.ndarray of floats
            The Hernquist acceleration vector of each particle
        """
        
        ### **** Store the magnitude of the position vector
        rmag = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
        
        ### *** Store the Acceleration
        b = -self.G * M
        c = rmag * (r_a + rmag)**2
        
        Hern =  (b/c) * r #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r):# it is easiest if you take as an input a position VECTOR  r 
        """ 
        Function that calculates the acceleration of disk particles using the Miyamoto-Nagai
        1975 profile potential
        
        Parameters
        ----------
        r : np.ndarray of floats
            The position vector
        
        Output
        ------
        Miya : The Miyamoto-Nagai acceleration vector of each disk particle
        """

        
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        R = np.sqrt(r[0]**2 + r[1]**2)
        B = r_d + np.sqrt(r[2]**2 + (r_d/5)**2)
        
        v = np.array([1,1, B/(np.sqrt(r[2]**2 + (r_d/5)**2))])
        
        b = -self.G * self.Mdisk
        c = (R**2 + B**2)**1.5
        
        Miya = (b/c) * r * v
       
        return Miya
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r): # input should include the position vector, r
        """
        Function that adds together the acceleration components of the bulge, disk, and halo
        
        Parameters
        ----------
        r : np.ndarray of floats
            the position vector
        
        Output
        ------
        total : np.ndarray of floats
            the vector of the sum of each galaxy component's acceleration
        """

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        bulgeAccel = self.HernquistAccel(self.Mbulge, self.rbulge, r)
        haloAccel = self.HernquistAccel(self.Mhalo, self.rhalo, r)
        diskAccel = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r)
        
        total = bulgeAccel + haloAccel + diskAccel
        
        # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return total
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """ 
        Function that time steps the position and velocity vectors in order to get
        a new position and velocity vector after the time step.
        
        Parameters
        ----------
        dt : float
            the time interval for integration
        r : np.array of floats
            the starting position vector
        v : np.array of floats
            the starting velocity vector
            
        Output
        ------
        rnew : np.array of floats
            the calculated position vector after timestepping
        vnew : np.array of floats
            the calculated velocity vector after timestepping
        """
        
        # predict the position at the next half timestep
        rhalf = r + v * dt/2
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        ahalf = self.M31Accel(rhalf)
        vnew = v + ahalf * dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew * dt/2
        
        return rnew, vnew # **** return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """ 
        Function to integrate the orbit of M33
        """

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        rows = int(tmax/dt) + 2
        orbit = np.zeros((rows, 7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t < tmax):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t = t + dt
           
            # **** store the new time in the first column of the ith row
            orbit[i,0] = t
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            rnew,vnew = self.LeapFrog(dt, orbit[i-1, 1:4], orbit[i-1, 4:7])
         
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            orbit[i, 1:4] = tuple(rnew)
            orbit[i, 4:7] = tuple(vnew)
            
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i = i+1
        
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

# Part 5 Analysis

M33 = M33AnalyticOrbit("M33_AnalyticalOrbit.txt")

M33.OrbitIntegration(0, 0.01, 10)

# Read in the data files for the orbits of each galaxy
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt

M31_AnaData = np.genfromtxt("M33_AnalyticalOrbit.txt", names=True)
M31_data = np.genfromtxt("Orbit_M31.txt", names=True)
M33_data = np.genfromtxt("Orbit_M33.txt", names=True)

t_M31A  = M31_AnaData['t']
x_M31A  = M31_AnaData['x']
y_M31A  = M31_AnaData['y']
z_M31A  = M31_AnaData['z']
vx_M31A = M31_AnaData['vx']
vy_M31A = M31_AnaData['vy']
vz_M31A = M31_AnaData['vz']

t_M31  = M31_data['t']
x_M31  = M31_data['x']
y_M31  = M31_data['y']
z_M31  = M31_data['z']
vx_M31 = M31_data['vx']
vy_M31 = M31_data['vy']
vz_M31 = M31_data['vz']

t_M33  = M33_data['t']
x_M33  = M33_data['x']
y_M33  = M33_data['y']
z_M33  = M33_data['z']
vx_M33 = M33_data['vx']
vy_M33 = M33_data['vy']
vz_M33 = M33_data['vz']


# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
def VectorDifference(a1, b1, c1, a2, b2, c2):
    """
    This function takes the difference between two vectors and returns the magnitude
    of that difference
    Inputs
    ----------
    a1, a2: float or array of floats
    b1, b2: float or array of floats
    c1, c2: float or array of floats
        All are components of 2 different vectors

    Outputs
    -------
    mag: float
        The magnitude of the distance between the two given vectors

    """
    aDiff = a2-a1
    bDiff = b2-b1
    cDiff = c2-c1
    mag = np.sqrt(aDiff**2 + bDiff**2 + cDiff**2)
    return mag

# Determine the magnitude of the relative position and velocities 

# of M33 and M31
M33_M31_separation = VectorDifference(x_M33, y_M33, z_M33, x_M31, y_M31, z_M31)
M33_M31_vel = VectorDifference(vx_M33, vy_M33, vz_M33, vx_M31, vy_M31, vz_M31)

M33_M31_AnaSeparation = VectorDifference(0, 0, 0, x_M31A, y_M31A, z_M31A)
M33_M31_AnaVel = VectorDifference(0, 0, 0, vx_M31A, vy_M31A, vz_M31A)
# Zeros since Analytical data should already be the difference

# Plot the Orbit of the galaxies 
#################################
fig, ax = plt.subplots(1,2, figsize=(25,10))

ax[0].plot(t_M31A, M33_M31_AnaSeparation, color='blue',  label='MW-M31 Analytical')
ax[0].plot(t_M31, M33_M31_separation, color='red', label='M33-M31')
ax[0].set_xlabel('Time (Gyr)')
ax[0].set_ylabel('Separation (kpc)')
ax[0].set_title('Separation vs Time')
ax[0].legend()

# Plot the orbital velocities of the galaxies 
#################################
ax[1].plot(t_M31A, M33_M31_AnaVel, color='blue',  label='MW-M31 Analytical')
ax[1].plot(t_M31, M33_M31_vel, color='red', label='M33-M31')
ax[1].set_xlabel('Time (Gyr)')
ax[1].set_ylabel('Relative Speed (km/s)')
ax[1].set_title('Relative Velocity vs Time')
ax[1].legend()

plt.savefig("Hmwk7_Plots.png")

