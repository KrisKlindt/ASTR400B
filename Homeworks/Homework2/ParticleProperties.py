# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 20:20:38 2025

@author: krist
"""

# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from ReadFile import Read # import Read to get data from file

def ParticleInfo(filename, pType, pNum):
    """
    This function takes a particle type and the number of that particle type and
    returns the magnitude of the distance, velocity, and mass of that particle in
    units of kpc, km/s, and solar masses, respectively
    Imputs:
        filename : (string) The name of the file to be opened
        pType : (float) The type of particle (1 for Dark Matter, 2 for Disk Star, 3 for Bulge Stars)
        pNum : (int) The number of the particle

    Returns:
        magDistance: (astropy quantity) The magnitude of the given particle's distance from the galactic center (in kpc)
        magVelocity: (astropy quantity) The magnitude of the particle's total velocity (in km/s)
        m: (astropy quantity) The mass of the particle (in M_sun)
    """
    
    time, numP, data = Read("MW_000.txt") # gets the time, total number of particles, and data from the file
    
    indices = np.where(data['type'] == pType) # gets the indices where the particle type is equal to what is given
    index = indices[0][0] + pNum - 1 # finds the index of the given particle number
    
    m = data['m'][index] # gets the mass of the given particle
    x, y, z = data['x'][index], data['y'][index], data['z'][index] # gets the x, y, and z values of the particle
    vx, vy, vz = data['vx'][index], data['vy'][index], data['vz'][index] # gets the x, y, and z values of the particle
    
    m = m*10**10*u.M_sun # converts mass into correct solar mass units
    x, y, z = x*u.kpc, y*u.kpc, z*u.kpc # gives units of kpc to x, y, z
    vx, vy, vz = vx*u.km/u.s, vy*u.km/u.s, vz*u.km/u.s # gives units of km/s to vx, vy, vz
    
    magDistance = np.around(np.sqrt(x**2 + y**2 + z**2), 3) # calculates the magnitude of the distance from galactic center
    magVelocity = np.around(np.sqrt(vx**2 + vy**2 + vz**2), 3) # calculates the magnitude of the velocity
    
    return magDistance, magVelocity, m
    
distance, velocity, mass = ParticleInfo("MW_000.txt", 2.0, 100) # calls the function

print("Distance: " + str(distance))
print("Velocity: " + str(velocity))
print("Mass: " + str(mass))

lightyears = np.around(distance.to(u.lyr), 3)
print(str(lightyears))

