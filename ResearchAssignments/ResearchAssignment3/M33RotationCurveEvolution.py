# -*- coding: utf-8 -*-

# M33RotationCurveEvolution

# The topic of my research project is how the internal stellar kinematics of galaxies
# evolve due to tides from a massive host.
# 
# In particular, I will be looking at how M33's stellar rotation curve evolves as the
# Milky Way and M31 merge.
#
# One method I will use to do this is by overplotting the analytical stellar rotation curve 
# derived from the mass profile at several points in M33's orbit, in particular its initial 
# rotation curve and the rotation curve near each relative pericenter and apocenter of M33's orbit,
# as well as one more point that looks like something interesting could be occuring.

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt

# my modules
from ReadFile import Read
from CenterOfMass2 import CenterOfMass
from MassProfile import MassProfile
from GalaxyMass import ComponentMass

# Snapshots: 0, 65, 190, 300, 380, 420, 475, 525, 580, 625

# Function to compute a rotation curve given snapshot and Jacobi Radius
def rotationCurve(snap, rj):
    """
    

    Parameters
    ----------
    snap : TYPE
        DESCRIPTION.
    rj : TYPE
        DESCRIPTION.

    Returns
    -------
    Vcirc : np.ndarray
        An array containing the circular velocity corresponding to a radius, in km/s
    """
    # create a mass profile object for M33
    M33 = MassProfile('M33', snap) # MassProfile gets the right string for given snapshot
    
    # create an array of positions
    rr = np.arange(0.01, rj, 0.1) # From 0.01 to rj, spaced by 0.1
    
    # Circular Velocity Profile
    Vcirc = M33.circularVelocityTotal(rr)
    
    return rr, Vcirc

# Function to compute Jacobi Radius
# rj = r * (Msat/2Mhost) ^ 1/3
def jacobiRadius(Msat, Mhost, r):
    """
    

    Parameters
    ----------
    Msat : TYPE
        DESCRIPTION.
    Mhost : TYPE
        DESCRIPTION.
    r : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # group masses
    a = (Msat/(2*Mhost))**(1/3)
    
    # do calculation
    return r*a

# Function to compute the total mass of a galaxy
def galaxyTotalMass(filename):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # Use ComponentMass to get masses for each particle type (units 10^12 Msun)
    
    # Check if galaxy is M33, since M33 doesn't have bulge particles (3)
    if "M33" in filename:
        Mhalo = ComponentMass(filename, 1)
        Mdisk = ComponentMass(filename, 2)
        
        # Add together, return
        return (Mhalo + Mdisk) * 1e12 # to put in units of Msun
    else:
        Mhalo = ComponentMass(filename, 1)
        Mdisk = ComponentMass(filename, 2)
        Mbulge = ComponentMass(filename, 3)
        
        # Add together, return
        return (Mhalo + Mdisk + Mbulge) * 1e12 # to put in units of Msun
    
# Function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position for two 
# galaxies
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

# Function to calculate the distance between the centers of mass of two galaxies
def separation(filename1, filename2):
    """
    

    Parameters
    ----------
    filename1 : TYPE
        DESCRIPTION.
    filename2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # Create CenterOfMass objects for each galaxy, using disk particles
    g1 = CenterOfMass(filename1, 2)
    g2 = CenterOfMass(filename2, 2)
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    if "M33" in filename1:
        volDec1 = 4.0
        volDec2 = 2.0
    elif "M33" in filename2:
        volDec1 = 2.0
        volDec2 = 4.0
    else:
        # not sure if this will be needed, have just in case
        volDec1 = 2.0
        volDec2 = 2.0
    
    # Initialize an instance of CenterOfMass class, using disk particles
    g1COM = CenterOfMass(filename1, 2)
    g2COM = CenterOfMass(filename2, 2)
     
    # Store the COM pos
    g1COM_pos = g1COM.COM_P(delta, volDec1)
    g2COM_pos = g2COM.COM_P(delta, volDec2)
    
    # Calculate magnitude separation between centers of mass
    sep = VectorDifference(g1COM_pos[0].value, g1COM_pos[1].value, g1COM_pos[2].value,
                           g2COM_pos[0].value, g2COM_pos[1].value, g2COM_pos[2].value)
 
    return sep
        
    
    
    
# calculate total masses for M33, M31, Milky Way
massM33 = galaxyTotalMass("M33/M33_000.txt") # change to highRes later
massM31 = galaxyTotalMass("M31/M31_000.txt")
massMW = galaxyTotalMass("MW/MW_000.txt")

massMerged = massM31 + massMW

# calculate Jacobi Radii for each snapshot
# M31 alone is the host 
rj1 = jacobiRadius(massM33, massM31, separation("M33/M33_000.txt", "M31/M31_000.txt")) # change to highRes later
rj2 = jacobiRadius(massM33, massM31, separation("M33/M33_065.txt", "M31/M31_065.txt"))
rj3 = jacobiRadius(massM33, massM31, separation("M33/M33_190.txt", "M31/M31_190.txt"))
rj4 = jacobiRadius(massM33, massM31, separation("M33/M33_300.txt", "M31/M31_300.txt"))
rj5 = jacobiRadius(massM33, massM31, separation("M33/M33_380.txt", "M31/M31_380.txt"))
rj6 = jacobiRadius(massM33, massM31, separation("M33/M33_420.txt", "M31/M31_420.txt"))

# M31, MW have merged, treat as one big host
# Since MW, M31 are merged, COM of one = COM of other
rj7 = jacobiRadius(massM33, massMerged, separation("M33/M33_475.txt", "M31/M31_475.txt"))
rj8 = jacobiRadius(massM33, massMerged, separation("M33/M33_525.txt", "M31/M31_525.txt"))
rj9 = jacobiRadius(massM33, massMerged, separation("M33/M33_580.txt", "M31/M31_580.txt"))
rj10 = jacobiRadius(massM33, massMerged, separation("M33/M33_625.txt", "M31/M31_625.txt"))

# Create an array to hold positions + rotation curves to loop through later
curves = np.zeros(2,10)

# call rotationCurve for each snapshot, add to curves
rr1, rC1 = rotationCurve(0, rj1)
curves[0][0] = rr1
curves[1][0] = rC1

rr2, rC2 = rotationCurve(65, rj2)
curves[0][1] = rr2
curves[1][1] = rC2

rr3, rC3 = rotationCurve(190, rj3)
curves[0][2] = rr3
curves[1][2] = rC3

rr4, rC4 = rotationCurve(300, rj4)
curves[0][3] = rr4
curves[1][3] = rC4

rr5, rC5 = rotationCurve(380, rj5)
curves[0][4] = rr5
curves[1][4] = rC5

rr6, rC6 = rotationCurve(420, rj6)
curves[0][5] = rr6
curves[1][5] = rC6

rr7, rC7 = rotationCurve(475, rj7)
curves[0][6] = rr7
curves[1][6] = rC7

rr8, rC8 = rotationCurve(525, rj8)
curves[0][7] = rr8
curves[1][7] = rC8

rr9, rC9 = rotationCurve(580, rj9)
curves[0][8] = rr9
curves[1][8] = rC9

rr10, rC10 = rotationCurve(625, rj10)
curves[0][0] = rr10
curves[1][0] = rC10

# create array for time of each snapshot
sToT = 10/0.7/1000 # so T is in Gyrs
times = np.array([0, 65 * sToT, 190 * sToT, 300 * sToT, 380 * sToT, 420 * sToT, 475 * sToT,
                  525 * sToT, 580 * sToT, 625 * sToT])

# 


# Create plot
fig = plt.figure(figsize=(12,10))
ax = plt.subplot(111)

# loop through curves, plot each curve
for i in range(10):
    plt.plot(curves[0][i], curves[1][i], label = f"Time = {times[i]:.2f} (Gyrs)")


# Add axis labels and title
plt.title("M33's Rotation Curve Evolution")
plt.xlabel('Radius (kpc)', fontsize=22)
plt.ylabel('Circular Velocity (km/s) ', fontsize=22)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# Save file
plt.savefig('M33_RotationCurveEvolution.png')


