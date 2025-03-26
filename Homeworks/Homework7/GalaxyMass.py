# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 16:14:22 2025

@author: krist
"""

# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from ReadFile import Read # import Read to read in the file

def ComponentMass(filename, pType):
    """
    This function will calculate and return the total mass of the galaxy component
    that is specified in the inputs
    Inputs:
        filename: (string) The name of the file to be opened
        pType: (float) Determines which particle is specified
            (1.0 for Halo Type, 2.0 for Disk Type, 3.0 for Bulge Type)
    Outputs:
        totalMass: (float) The total mass of the specified particle type
            within the galaxy (units of 10^12 Msun)
    """
    
    time, numP, data = Read(filename) # gets the time, total number of particles, and data from the file

    indices = np.where(data['type'] == pType) # gets the indices where the particle type is equal to what is given
    
    masses = data['m'][indices] # get a list of all the masses of each specified particle
    
    totalMass = 0 # initialize the total mass
    for m in masses:
        totalMass += m # adds up each individual mass to the total mass
        
    totalMass = np.round((totalMass * 10**(-2)), 3) # since mass was originally in units of
    # 10^10 Msun, need to multiply by 10^-2 to get totalMass in units of 10^12 Msun
    
    return totalMass

"""
MWHaloMass = ComponentMass("MW_000.txt", 1.0) # gets the total mass of MW halo particles
MWDiskMass = ComponentMass("MW_000.txt", 2.0) # gets the total mass of MW disk particles
MWBulgeMass = ComponentMass("MW_000.txt", 3.0) # gets the total mass of MW bulge particles

M31HaloMass = ComponentMass("M31_000.txt", 1.0) # gets the total mass of M31 halo particles
M31DiskMass = ComponentMass("M31_000.txt", 2.0) # gets the total mass of M31 disk particles
M31BulgeMass = ComponentMass("M31_000.txt", 3.0) # gets the total mass of M31 bulge particles

M33HaloMass = ComponentMass("M33_000.txt", 1.0) # gets the total mass of M33 halo particles
M33DiskMass = ComponentMass("M33_000.txt", 2.0) # gets the total mass of M33 disk particles
M33BulgeMass = ComponentMass("M33_000.txt", 3.0) # gets the total mass of M33 bulge particles

print(MWHaloMass)
print(MWDiskMass)
print(MWBulgeMass)

print(M31HaloMass)
print(M31DiskMass)
print(M31BulgeMass)

print(M33HaloMass)
print(M33DiskMass)
print(M33BulgeMass)

MWTotalMass = MWHaloMass+MWDiskMass+MWBulgeMass # adds together all the mass in MW
print(MWTotalMass)
M31TotalMass = M31HaloMass+M31DiskMass+M31BulgeMass # adds together all the mass in M31
print(M31TotalMass)
M33TotalMass = M33HaloMass+M33DiskMass+M33BulgeMass # adds together all the mass in M33
print(M33TotalMass)

# f_bar = BaryonicMass/TotalMass

MWBaryonMass = MWDiskMass+MWBulgeMass # adds togther all the baryonic mass in MW
print((MWBaryonMass/MWTotalMass)) # print f_bar for MW
M31BaryonMass = M31DiskMass+M31BulgeMass # adds together all the baryonic mass in M31
print((M31BaryonMass/M31TotalMass)) # print f_bar for M31
M33BaryonMass = M33DiskMass+M33BulgeMass # adds together all the baryonic mass in M33
print((M33BaryonMass/M33TotalMass)) # print f_bar for M33

LGMass = MWTotalMass+M31TotalMass+M33TotalMass # calculates total mass of local group
print(LGMass)

LGBaryonMass = MWBaryonMass+M31BaryonMass+M33BaryonMass # calculates baryonic mass of local group
print(LGBaryonMass/LGMass) # print f_bar of local group
"""