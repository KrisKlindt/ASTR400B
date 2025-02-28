

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass




def OrbitCOM(galaxy, start, end, n=5):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
        galaxy: str
            The name of the galaxy
        start: int
            The number of the first snap shot to be read in
        end: int
            The number of the last snap shot to be read in
        n: int
            The interval over which the COM will be returned
          
    outputs: 
        Doesn't return anything, but creates and writes to a file that is stored
    """
    
    # compose the filename for output
    fileout = "Orbit_" + galaxy + ".txt"
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    if galaxy == "M33":
        volDec = 4.0
    else:
        volDec = 2.0
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end+1, n)
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros((len(snap_ids), 7))
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        # Reconstruct filename from inputs
        ilbl = "000" + str(snap_id)
        ilbl = ilbl[-3:]
        filename = f"{galaxy}/{galaxy}_{ilbl}.txt"
        
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename, 2)
        
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM_pos = COM.COM_P(delta, volDec)
        COM_v = COM.COM_V(COM_pos[0], COM_pos[1], COM_pos[2])
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        orbit[i, 0] = COM.time.to(u.Gyr).value
        orbit[i, 1] = COM_pos[0].value # x
        orbit[i, 2] = COM_pos[1].value # y
        orbit[i, 3] = COM_pos[2].value # z
        orbit[i, 4] = COM_v[0].value # vx
        orbit[i, 5] = COM_v[1].value # vy
        orbit[i, 6] = COM_v[2].value # vz
        
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

OrbitCOM("MW", 0, 800)
OrbitCOM("M31", 0, 800)
OrbitCOM("M33", 0, 800)


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt

MW_data = np.genfromtxt("Orbit_MW.txt", names=True)
M31_data = np.genfromtxt("Orbit_M31.txt", names=True)
M33_data = np.genfromtxt("Orbit_M33.txt", names=True)

t_MW  = MW_data['t']
x_MW  = MW_data['x']
y_MW  = MW_data['y']
z_MW  = MW_data['z']
vx_MW = MW_data['vx']
vy_MW = MW_data['vy']
vz_MW = MW_data['vz']

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

# of MW and M31
MW_M31_separation = VectorDifference(x_MW, y_MW, z_MW, x_M31, y_M31, z_M31)
MW_M31_vel = VectorDifference(vx_MW, vy_MW, vz_MW, vx_M31, vy_M31, vz_M31)

# of M33 and M31
M33_M31_separation = VectorDifference(x_M33, y_M33, z_M33, x_M31, y_M31, z_M31)
M33_M31_vel = VectorDifference(vx_M33, vy_M33, vz_M33, vx_M31, vy_M31, vz_M31)


# Plot the Orbit of the galaxies 
#################################
fig, ax = plt.subplots(1,2, figsize=(25,10))

ax[0].plot(t_MW, MW_M31_separation, color='blue',  label='MW-M31')
ax[0].plot(t_M31, M33_M31_separation, color='red', label='M33-M31')
ax[0].set_xlabel('Time (Gyr)')
ax[0].set_ylabel('Separation (kpc)')
ax[0].set_title('Separation vs Time')
ax[0].legend()

# Plot the orbital velocities of the galaxies 
#################################
ax[1].plot(t_MW, MW_M31_vel, color='blue',  label='MW-M31')
ax[1].plot(t_M31, M33_M31_vel, color='red', label='M33-M31')
ax[1].set_xlabel('Time (Gyr)')
ax[1].set_ylabel('Relative Speed (km/s)')
ax[1].set_title('Relative Velocity vs Time')
ax[1].legend()

plt.savefig("Hmwk6_Plots.png")

