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

# I will then calculate the actual average circular velocities of the stars within radius bins
# and overplot those rotation curves for the same snapshots mentioned above. This should
# show that, as time passes, the actual circular velocity of the stars in the galaxy will
# look less and less like the analytical circular velocity, meaning that M33 is becoming
# less rotationally supported and more dispersion supported, implying that M31 and the
# Milky Way caused M33 to evolve from a dwarf irregular galaxy to a dwarf spheroidal/elliptical
# galaxy

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
from Lab9_Soln import RotateFrame

# Snapshots: 0, 65, 190, 300, 380, 420, 475, 525, 580, 625

# Function to compute a rotation curve given snapshot and Jacobi Radius
def rotationCurve(snap, rj):
    """
    

    Parameters
    ----------
    snap : float
        The desired snapshop number.
    rj : float
        The Jacobi radius (kpc).

    Returns
    -------
    Vcirc : np.ndarray
        An array containing the circular velocity corresponding to a radius, in km/s
    """
    # create a mass profile object for M33
    M33 = MassProfile('M33/M33', snap) # MassProfile gets the right string for given snapshot
    
    # create an array of positions
    rr = np.arange(0.01, rj, 0.1) # From 0.01 to rj, spaced by 0.1
    
    # Circular Velocity Profile
    Vcirc = M33.circularVelocityTotal(rr).value
    
    return rr, Vcirc

# Function to compute Jacobi Radius
# rj = r * (Msat/2Mhost) ^ 1/3
def jacobiRadius(Msat, Mhost, r):
    """
    Function to calculate the Jacobi radius of a satellite galaxy

    Parameters
    ----------
    Msat : float
        The mass of the satellite galaxy (Msun).
    Mhost : float
        The mass of the host galaxy (Msun).
    r : float
        The distance between the centers of mass of the two galaxies (kpc).

    Returns
    -------
    rj : float
        The Jacobi radius (kpc)

    """
    # group masses
    a = (Msat/(2*Mhost))**(1/3)
    
    # do calculation
    rj = r*a
    
    return rj

# Function to compute the total mass of a galaxy
def galaxyMassEnclosed(galaxy):
    """
    Function to compure the total mass of the galaxy

    Parameters
    ----------
    galaxy : string
        The name of the galaxy whose mass will be calculated.

    Returns
    -------
    Tmass : float
        The total mass of the galaxy

    """
    # Use ComponentMass to get masses for each particle type (units 10^12 Msun)
    
    snaps = np.array([0, 65, 190, 300, 380, 420, 475, 525, 580, 625])
    
    for i in range(snaps.size):
        # add a string of the filenumber to the value "000"
        ilbl = '000' + str(snaps[i])
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        
        filename = galaxy + f"_{ilbl}.txt"
    
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
            Tmass = (Mhalo + Mdisk + Mbulge) * 1e12 # to put in units of Msun)
            return Tmass
    
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
def sepAndPosAndVel(filename1, filename2):
    """
    Function to calculate the distance between the centers of mass of the two
    galaxies, as well return the position and velocity vectors of each disk particle in M33

    Parameters
    ----------
    filename1 : string
        The snapshot file to be used to create a CenterOfMass object. Will always be M33
    filename2 : string
        Second snapshot gile to be used to create a CenterOfMass object.

    Returns
    -------
    sep : float
        The magnitude of the separation of the centers of mass of the two galaxies
    r : np.ndarray of floats
        The position vectors of each disk particle in M33
    v : np.ndarray of floats
        The velocity vectors of each disk particle in M33

    """
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
    
    # Store COM velocity for M33
    COMV = g1COM.COM_V(g1COM_pos[0],g1COM_pos[1],g1COM_pos[2])
    
    # Determine positions of disk particles relative to COM 
    # g1 will always be M33
    xD = g1COM.x - g1COM_pos[0].value 
    yD = g1COM.y - g1COM_pos[1].value 
    zD = g1COM.z - g1COM_pos[2].value 

    # Determine velocities of disk particles relatiev to COM motion
    vxD = g1COM.vx - COMV[0].value 
    vyD = g1COM.vy - COMV[1].value 
    vzD = g1COM.vz - COMV[2].value 

    # Arrays for r and v 
    r = np.array([xD,yD,zD]).T # transposed for the Rotate Function later
    v = np.array([vxD,vyD,vzD]).T
    
    # Calculate magnitude separation between centers of mass
    sep = VectorDifference(g1COM_pos[0].value, g1COM_pos[1].value, g1COM_pos[2].value,
                           g2COM_pos[0].value, g2COM_pos[1].value, g2COM_pos[2].value)
    
    
 
    return sep, r, v
    

# Calculate Mass Enclosed of M31 or M31 and MW
# Assume total mass remains the same for all
    
# calculate total masses for M33, M31, Milky Way
massM33 = galaxyMassEnclosed("M33/M33") # change to highRes later
massM31 = galaxyMassEnclosed("M31/M31")
massMW = galaxyMassEnclosed("MW/MW")

massMerged = massM31 + massMW

sep1, r1, v1 = sepAndPosAndVel("M33/M33_000.txt", "M31/M31_000.txt")
sep2, r2, v2 = sepAndPosAndVel("M33/M33_065.txt", "M31/M31_065.txt")
sep3, r3, v3 = sepAndPosAndVel("M33/M33_190.txt", "M31/M31_190.txt")
sep4, r4, v4 = sepAndPosAndVel("M33/M33_300.txt", "M31/M31_300.txt")
sep5, r5, v5 = sepAndPosAndVel("M33/M33_380.txt", "M31/M31_380.txt")
sep6, r6, v6 = sepAndPosAndVel("M33/M33_420.txt", "M31/M31_420.txt")
sep7, r7, v7 = sepAndPosAndVel("M33/M33_475.txt", "M31/M31_475.txt")
sep8, r8, v8 = sepAndPosAndVel("M33/M33_525.txt", "M31/M31_525.txt")
sep9, r9, v9 = sepAndPosAndVel("M33/M33_580.txt", "M31/M31_580.txt")
sep10, r10, v10 = sepAndPosAndVel("M33/M33_625.txt", "M31/M31_625.txt")

# calculate Jacobi Radii for each snapshot
# M31 alone is the host 
rj1 = jacobiRadius(massM33, massM31, sep1) # change to highRes later
rj2 = jacobiRadius(massM33, massM31, sep2)
rj3 = jacobiRadius(massM33, massM31, sep3)
rj4 = jacobiRadius(massM33, massM31, sep4)
rj5 = jacobiRadius(massM33, massM31, sep5)
rj6 = jacobiRadius(massM33, massM31, sep6)

# M31, MW have merged, treat as one big host
# Since MW, M31 are merged, COM of one = COM of other
rj7 = jacobiRadius(massM33, massMerged, sep7)
rj8 = jacobiRadius(massM33, massMerged, sep8)
rj9 = jacobiRadius(massM33, massMerged, sep9)
rj10 = jacobiRadius(massM33, massMerged, sep10)

# Store Jacobi Radii
jacobiRadii = [rj1,rj2,rj3,rj4,rj5,rj6,rj7,rj8,rj9,rj10]

# Create an array to hold positions + rotation curves to loop through later
curves = []
radii = []

# call rotationCurve for each snapshot, add to curves
rr1, rC1 = rotationCurve(0, rj1)
radii.append(rr1)
curves.append(rC1)

rr2, rC2 = rotationCurve(65, rj2)
radii.append(rr2)
curves.append(rC2)

rr3, rC3 = rotationCurve(190, rj3)
radii.append(rr3)
curves.append(rC3)

rr4, rC4 = rotationCurve(300, rj4)
radii.append(rr4)
curves.append(rC4)

rr5, rC5 = rotationCurve(380, rj5)
radii.append(rr5)
curves.append(rC5)

rr6, rC6 = rotationCurve(420, rj6)
radii.append(rr6)
curves.append(rC6)

rr7, rC7 = rotationCurve(475, rj7)
radii.append(rr7)
curves.append(rC7)

rr8, rC8 = rotationCurve(525, rj8)
radii.append(rr8)
curves.append(rC8)

rr9, rC9 = rotationCurve(580, rj9)
radii.append(rr9)
curves.append(rC9)

rr10, rC10 = rotationCurve(625, rj10)
radii.append(rr10)
curves.append(rC10)

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
    plt.plot(radii[i], curves[i], label = f"Time = {times[i]:.2f} (Gyrs)")


# Add axis labels and title
plt.title("M33's Analytical Rotation Curve Evolution")
plt.xlabel('Radius (kpc)', fontsize=22)
plt.ylabel('Circular Velocity (km/s) ', fontsize=22)
plt.legend()

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# Save file
plt.savefig('M33_RotationCurveEvolution.png')


# compute the rotated position and velocity vectors
rn1, vn1 = RotateFrame(r1,v1)
rn2, vn2 = RotateFrame(r2,v2)
rn3, vn3 = RotateFrame(r3,v3)
rn4, vn4 = RotateFrame(r4,v4)
rn5, vn5 = RotateFrame(r5,v5)
rn6, vn6 = RotateFrame(r6,v6)
rn7, vn7 = RotateFrame(r7,v7)
rn8, vn8 = RotateFrame(r8,v8)
rn9, vn9 = RotateFrame(r9,v9)
rn10, vn10 = RotateFrame(r10,v10)

# Put these values into arrays
rotatedRs = [rn1,rn2,rn3,rn4,rn5,rn6,rn7,rn8,rn9,rn10]
rotatedVs = [vn1,vn2,vn3,vn4,vn5,vn6,vn7,vn8,vn9,vn10]

# total magnitude
rtot1 = np.sqrt(rn1[:,0]**2 + rn1[:,1]**2 + rn1[:,2]**2)
rtot2 = np.sqrt(rn2[:,0]**2 + rn2[:,1]**2 + rn2[:,2]**2)
rtot3 = np.sqrt(rn3[:,0]**2 + rn3[:,1]**2 + rn3[:,2]**2)
rtot4 = np.sqrt(rn4[:,0]**2 + rn4[:,1]**2 + rn4[:,2]**2)
rtot5 = np.sqrt(rn5[:,0]**2 + rn5[:,1]**2 + rn5[:,2]**2)
rtot6 = np.sqrt(rn6[:,0]**2 + rn6[:,1]**2 + rn6[:,2]**2)
rtot7 = np.sqrt(rn7[:,0]**2 + rn7[:,1]**2 + rn7[:,2]**2)
rtot8 = np.sqrt(rn8[:,0]**2 + rn8[:,1]**2 + rn8[:,2]**2)
rtot9 = np.sqrt(rn9[:,0]**2 + rn9[:,1]**2 + rn9[:,2]**2)
rtot10 = np.sqrt(rn10[:,0]**2 + rn10[:,1]**2 + rn10[:,2]**2)

# velocity in cylindrical coordinates. 

# radius 
rho1 = np.sqrt(rn1[:,0]**2 + rn1[:,1]**2) 
rho2 = np.sqrt(rn2[:,0]**2 + rn2[:,1]**2) 
rho3 = np.sqrt(rn3[:,0]**2 + rn3[:,1]**2) 
rho4 = np.sqrt(rn4[:,0]**2 + rn4[:,1]**2) 
rho5 = np.sqrt(rn5[:,0]**2 + rn5[:,1]**2) 
rho6 = np.sqrt(rn6[:,0]**2 + rn6[:,1]**2) 
rho7 = np.sqrt(rn7[:,0]**2 + rn7[:,1]**2) 
rho8 = np.sqrt(rn8[:,0]**2 + rn8[:,1]**2) 
rho9 = np.sqrt(rn9[:,0]**2 + rn9[:,1]**2) 
rho10 = np.sqrt(rn10[:,0]**2 + rn10[:,1]**2) 

rhos = [rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8,rho9,rho10]

"""
# radial velocity 
vr1 = (rn1[:,0] * vn1[:,0] + rn1[:,1] * vn1[:,1]) / rho1
vr2 = (rn2[:,0] * vn2[:,0] + rn2[:,1] * vn2[:,1]) / rho2
vr3 = (rn3[:,0] * vn3[:,0] + rn3[:,1] * vn3[:,1]) / rho3
vr4 = (rn4[:,0] * vn4[:,0] + rn4[:,1] * vn4[:,1]) / rho4
vr5 = (rn5[:,0] * vn5[:,0] + rn5[:,1] * vn5[:,1]) / rho5
vr6 = (rn6[:,0] * vn6[:,0] + rn6[:,1] * vn6[:,1]) / rho6
vr7 = (rn7[:,0] * vn7[:,0] + rn7[:,1] * vn7[:,1]) / rho7
vr8 = (rn8[:,0] * vn8[:,0] + rn8[:,1] * vn8[:,1]) / rho8
vr9 = (rn9[:,0] * vn9[:,0] + rn9[:,1] * vn9[:,1]) / rho9
vr10 = (rn10[:,0] * vn10[:,0] + rn10[:,1] * vn10[:,1]) / rho10
"""

# azimuthal velocity
vphi1 = (rn1[:,0] *  vn1[:,1] - rn1[:,1] * vn1[:,0]) / rho1
vphi2 = (rn2[:,0] *  vn2[:,1] - rn2[:,1] * vn2[:,0]) / rho2
vphi3 = (rn3[:,0] *  vn3[:,1] - rn3[:,1] * vn3[:,0]) / rho3
vphi4 = (rn4[:,0] *  vn4[:,1] - rn4[:,1] * vn4[:,0]) / rho4
vphi5 = (rn5[:,0] *  vn5[:,1] - rn5[:,1] * vn5[:,0]) / rho5
vphi6 = (rn6[:,0] *  vn6[:,1] - rn6[:,1] * vn6[:,0]) / rho6
vphi7 = (rn7[:,0] *  vn7[:,1] - rn7[:,1] * vn7[:,0]) / rho7
vphi8 = (rn8[:,0] *  vn8[:,1] - rn8[:,1] * vn8[:,0]) / rho8
vphi9 = (rn9[:,0] *  vn9[:,1] - rn9[:,1] * vn9[:,0]) / rho9
vphi10 = (rn10[:,0] *  vn10[:,1] - rn10[:,1] * vn10[:,0]) / rho10

vphis = [vphi1,vphi2,vphi3,vphi4,vphi5,vphi6,vphi7,vphi8,vphi9,vphi10]
    
# Determine the mean vphi per radius

# Initialize Empty Array for Velocity 
# (same size as radial array)
Vel1= np.zeros(np.size(rr1))
Vel2= np.zeros(np.size(rr2))
Vel3= np.zeros(np.size(rr3))
Vel4= np.zeros(np.size(rr4))
Vel5= np.zeros(np.size(rr5))
Vel6= np.zeros(np.size(rr6))
Vel7= np.zeros(np.size(rr7))
Vel8= np.zeros(np.size(rr8))
Vel9= np.zeros(np.size(rr9))
Vel10= np.zeros(np.size(rr10))

Vels = [Vel1,Vel2,Vel3,Vel4,Vel5,Vel6,Vel7,Vel8,Vel9,Vel10]

# compute the mean vphi in radial bins

for i in range(len(radii)):
    for a in range(len(radii[i])):
        index = np.where((rhos[i] > 0.5*a) & (rhos[i] < 0.5*i+0.5)) # walking out in radial bins
        Vels[i][a] = np.mean(np.abs(vphis[i][index])) # mean velocity


fig = plt.figure(figsize=(12,10))
ax = plt.subplot(111)

# looking at MW edge on along x axis, vy is line of sight velocity

# Add the circular velocity

# loop through Vels, plot each curve
for i in range(10):
    plt.plot(radii[i], Vels[i], label = f"Time = {times[i]:.2f} (Gyrs)")


# Add axis labels and title
plt.title("M33's Mean v$_\phi$ Rotation Curve Evolution")
plt.xlabel('R (kpc)', fontsize=22)
plt.ylabel(r'v$_\phi$ (km/s)', fontsize=22)
plt.legend()

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# Save file
plt.savefig('M33_RotationCurve_Vphi.png')


# The mean cylindrical vphi deviates in the inner regions from the mass-derived 
# rotation curve because of the bar. 
