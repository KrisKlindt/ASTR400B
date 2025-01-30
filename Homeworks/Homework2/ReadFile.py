# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 15:34:39 2025

@author: krist
"""

# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units

def Read(filename):
    """
    This functions takes in a file name and organizes the data within the file
    based on the type of data.
    Inputs: 
        filename: (string) Is the name of the file that is to be opened
    Returns:
        time: (astropy quantity) Is the time in Myr given in the file
        numParticles: (float) Is the number of particles in the data
        data: The organized array of data from the file
    """
    
    file = open(filename, 'r') # open the file
    
    line1 = file.readline() # reads in the first line of the file
    label, value = line1.split() # takes the label and value of the time
    time = float(value)*u.Myr # puts the time in units of Myr
    
    line2 = file.readline() # reads in the second line of the file
    particleLabel, numParticles = line2.split() # takes the label and value of 
                                                # the number of particles
    numParticles = float(numParticles) # converts numParticles to a float for future use
    
    file.close() # closes the file
    
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)
    # Puts the data contained within the file in arrays based on the column in
    # which the data is in, uses the fourth line of the data as headers
    
    return time, numParticles, data 
    
"""
time, numP, data = Read("MW_000.txt")

print(data['x'][0])
print(data['y'][5])
print(data['z'][4])
print(data['vx'][2])
print(data['vy'][1])
print(data['vz'][3])
"""