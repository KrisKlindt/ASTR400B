# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:11:33 2025

@author: krist
"""

# import relevant modules
import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import G
from ReadFile import Read
import astropy.units as u
from CenterOfMass import CenterOfMass

class MassProfile:
# Class to define the mass profile and rotation curve of a given galaxy 
# and simulation snapshot

    def __init__(self, galaxy, snap):
        ''' 
        Initialize the class by reading data from a snapshot. 
            
        PARAMETERS
        ----------
        galaxy : str
            The name of the galaxy (MW, M31, or M33)
        snap : int
            The number of the desired snapshot file
        '''
        # Save the name of the galaxy
        self.gname = galaxy
        
        # Reconstruct filename from inputs
        ilbl = "000" + str(snap)
        ilbl = ilbl[-3:]
        self.filename = f"{galaxy}_{ilbl}.txt"
        
        # read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename)
        
        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc
        self.m = self.data['m']
        
        # convert G to appropriate units to get km/s velocity
        # Store as an instance variable to avoid duplicated code
        self.G = G.to(u.kpc * u.km**2 / u.s**2 / u.Msun)
        
    def MassEnclosed(self, ptype, radii):
        """
        Finds the mass enclosed within each radius in an array of radii
        
        Parameters
        ----------
        ptype : int
            The type of the desired particle (1 = Halo, 2 = Disk, 3 = Bulge)
        radii : array of floats
            The radii for which an enclosed mass is to be calculated
        Returns
        -------
        enclosed_mass: array of astropy quantities
            The masses enclosed with each radius
        """
        
        #create an array to store indexes of particles of desired Ptype                                
        index = np.where(self.data['type'] == ptype)
        # Store values for desired particle type
        x = self.x[index]
        y = self.y[index]
        z = self.z[index]
        m = self.m[index]
        
        # Create a CenterOfMass object
        com = CenterOfMass(self.filename, 2) # 2 for disk particles
        x_COM, y_COM, z_COM = com.COM_P(0.1) # call COM_P to get 3D position for COM
        
        # Center relative to COM
        x_part = x - x_COM
        y_part = y - y_COM
        z_part = z - z_COM
        r_part = np.sqrt(x_part**2 + y_part**2 + z_part**2)  # distance of each particle from COM
        
        # Ensure 'radii' is an astropy quantity in kpc
        radii = radii * u.kpc
        
        # Initialize array for enclosed mass
        enclosed_mass = np.zeros(len(radii))
        
        # Loop over radii and calculate each enclosed mass
        for i in range(len(radii)):
            indices = np.where(r_part < radii[i])
            enclosed_mass[i] = np.sum(m[indices])
        
        # Give correct units
        new_enclosed_mass = enclosed_mass * 1e10 * u.Msun
        return new_enclosed_mass
    
    def MassEnclosedTotal(self, radii):
        """
        Finds the total mass (halo + disk + bulge) enclosed within each radius in
        an array of radii

        Parameters
        ----------
        radii : array of floats
            The radii for which the total enclosed mass will be calculated
        Returns
        -------
        tot_enclosed_mass: array of astropy quantities
            The total enclosed mass at each radius
        """
        # Calculate mass of halo and disk
        m_halo = self.MassEnclosed(1, radii)
        m_disk = self.MassEnclosed(2, radii)
        
        # Account for fact that M33 does not have bulge particles
        if (self.gname != "M33"):
            m_bulge = self.MassEnclosed(3, radii)
        else:
            m_bulge = np.zeros(len(radii))
            
        # Add all masses together
        tot_enclosed_mass = m_halo+m_disk+m_bulge
        return tot_enclosed_mass
    
    def HernquistMass(self, radius, a, Mhalo):
        """ Function that defines the Hernquist dark matter mass profile 
        
        Inputs:
            radius: array of floats
                Galactocentric distance in kpc
            a: float
                scale radius of the Hernquist profile in kpc
            m_halo: float
                total halo mass 
        Ouputs:
            mass:  array of astropy quantities
                total mass within the input radius r in Msun
        """
        a = a * u.kpc
        radius = radius * u.kpc
        
        const = Mhalo #constants
        b = radius**2 /(a + radius)**2
        
        mass = const*b #Hernquist profile
        
        return mass
    
    def CircularVelocity(self, ptype, radii):
        """
        Function that calculates the circular velocity using a given ptype at given
        radii
        
        Parameters
        ----------
        ptype : int
            The type of the desired particle (1 = Halo, 2 = Disk, 3 = Bulge)
        radii : array of floats
            The radii at which the circular velocity is to be calculated
        Returns
        -------
        v_circ : array of astropy quantities
            The circular velocity at each radius in radii
        """
        # Get enclosed mass
        M_encl = self.MassEnclosed(ptype, radii)
        radii = radii * u.kpc
        
        # Calculate circular velocity
        # v = (G*M/r)**(1/2)
        v_circ = np.sqrt(self.G*M_encl/radii)
        #print(v_circ)
        return v_circ
    
    def CircularVelocityTotal(self, radii):
        """
        Function that calculates the circular velocity at given radii using the
        total mass enclosed within

        Parameters
        ----------
        radii : array of floats
            The radii at which the circular velocity will be calculated
        Returns
        -------
        v_circ_total : array of astropy quantities
            The circular velocity at each radius found using the total mass enclosed
            within each radius

        """
        # Calculate total mass enclosed
        total_mass = self.MassEnclosedTotal(radii)
            
        radii = radii * u.kpc
        
        # Calculate circular velocity
        # v = (G*M_total/r)**(1/2)
        v_circ_total = np.sqrt(self.G*total_mass/radii)
        return v_circ_total
    
    def HernquistVCirc(self, radius, a, Mhalo):
        """
        Function that caluclates the circular velocity at given radii using the
        Hernquist mass profile

        Parameters
        ----------
        radius: array of floats
            Galactocentric distance in kpc
        a: float
            scale radius of the Hernquist profile in kpc
        m_halo: float
            total halo mass in units of 1e12 Msun 
        Returns
        -------
        Hern_Vcirc : array of astropy quantities
            The Hernquist circular velocity at each radius
        """
        # Get Hernquist mass
        mass = self.HernquistMass(radius, a, Mhalo)
        # Give radius correct units
        radius = radius * u.kpc
        
        # Calculate circular velocity
        # v = G*M_Hernquist/r
        Hern_Vcirc = np.around(np.sqrt(self.G*mass/radius),2)
        return Hern_Vcirc
        
# Plots    
    
# Create MassProfile objects for each galaxy at Snapshot 0
MW = MassProfile("MW", 0)
M31 = MassProfile("M31", 0)
M33 = MassProfile("M33", 0)

# Define an array of radii from 0.1 to 30 kpc
r = np.linspace(0.1, 30, 100)


# Milky Way Mass Profile
# Mass of each component
MW_Halo_mass = MW.MassEnclosed(1, r)
MW_Disk_mass = MW.MassEnclosed(2, r)
MW_Bulge_mass = MW.MassEnclosed(3, r)

# Total mass
MW_Total_mass = MW.MassEnclosedTotal(r)

# Hernquist
MWTotal_Halo_mass = MW_Halo_mass[-1] # Use the last value (at highest r) as total halo mass 
a_guess_MW = 20  # kpc - Adjust as needed for best fit.

MW_Hern_mass = MW.HernquistMass(r, a_guess_MW, MWTotal_Halo_mass)

plt.figure(figsize=(8,6))
plt.semilogy(r, MW_Halo_mass, label='Halo')
plt.semilogy(r, MW_Disk_mass, label='Disk')
plt.semilogy(r, MW_Bulge_mass, label='Bulge')
plt.semilogy(r, MW_Total_mass, label='Total')
plt.semilogy(r, MW_Hern_mass, '--', label=f'Hernquist Halo (a={a_guess_MW} kpc)')

plt.title('Milky Way Mass Profile')
plt.xlabel('Radius (kpc)')
plt.ylabel('Log10 Enclosed Mass (Msun)')
plt.legend()
plt.savefig("Milky Way Mass Profile.png")



# M31 Mass Profile
# Mass of each component
M31_Halo_mass = M31.MassEnclosed(1, r)
M31_Disk_mass = M31.MassEnclosed(2, r)
M31_Bulge_mass = M31.MassEnclosed(3, r)

# Total mass
M31_Total_mass = M31.MassEnclosedTotal(r)

# Hernquist
M31Total_Halo_mass = M31_Halo_mass[-1] # Use the last value (at highest r) as total halo mass 
a_guess_M31 = 15  # kpc - Adjust as needed for best fit.

M31_Hern_mass = M31.HernquistMass(r, a_guess_M31, M31Total_Halo_mass)

plt.figure(figsize=(8,6))
plt.semilogy(r, M31_Halo_mass, label='Halo')
plt.semilogy(r, M31_Disk_mass, label='Disk')
plt.semilogy(r, M31_Bulge_mass, label='Bulge')
plt.semilogy(r, M31_Total_mass, label='Total')
plt.semilogy(r, M31_Hern_mass, '--', label=f'Hernquist Halo (a={a_guess_M31} kpc)')

plt.title('M31 Mass Profile')
plt.xlabel('Radius (kpc)')
plt.ylabel('Log10 Enclosed Mass (Msun)')
plt.legend()
plt.savefig("M31 Mass Profile.png")



# M33 Mass Profile
# Mass of each component
M33_Halo_mass = M33.MassEnclosed(1, r)
M33_Disk_mass = M33.MassEnclosed(2, r)

# Total mass
M33_Total_mass = M33.MassEnclosedTotal(r)

# Hernquist
M33Total_Halo_mass = M33_Halo_mass[-1] # Use the last value (at highest r) as total halo mass 
a_guess_M33 = 10  # kpc - Adjust as needed for best fit.

M33_Hern_mass = M33.HernquistMass(r, a_guess_M33, M33Total_Halo_mass)

plt.figure(figsize=(8,6))
plt.semilogy(r, M33_Halo_mass, label='Halo')
plt.semilogy(r, M33_Disk_mass, label='Disk')
plt.semilogy(r, M33_Total_mass, label='Total')
plt.semilogy(r, M33_Hern_mass, '--', label=f'Hernquist Halo (a={a_guess_M33} kpc)')

plt.title('M33 Mass Profile')
plt.xlabel('Radius (kpc)')
plt.ylabel('Log10 Enclosed Mass (Msun)')
plt.legend()
plt.savefig("M33 Mass Profile.png")



# Milky Way Rotation Curve
# Circular velocity of each component
MW_Halo_vcirc = MW.CircularVelocity(1, r)
MW_Disk_vcirc = MW.CircularVelocity(2, r)
MW_Bulge_vcirc = MW.CircularVelocity(3, r)

# Total circular velocity
MW_Total_vcirc = MW.CircularVelocityTotal(r)

MW_Hern_vcirc = MW.HernquistVCirc(r, a_guess_MW, MWTotal_Halo_mass)

plt.figure(figsize=(8,6))
plt.plot(r, MW_Halo_vcirc, label='Halo')
plt.plot(r, MW_Disk_vcirc, label='Disk')
plt.plot(r, MW_Bulge_vcirc, label='Bulge')
plt.plot(r, MW_Total_vcirc, label='Total')
plt.plot(r, MW_Hern_vcirc, '--', label=f'Hernquist Halo (a={a_guess_MW} kpc)')

plt.title('Milky Way Rotation Curve')
plt.xlabel('Radius (kpc)')
plt.ylabel('Circular Velocity (km/s)')
plt.legend()
plt.savefig("Milky Way Rotation Curve.png")



# M31 Rotation Curve
# Circular velocity of each component
M31_Halo_vcirc = M31.CircularVelocity(1, r)
M31_Disk_vcirc = M31.CircularVelocity(2, r)
M31_Bulge_vcirc = M31.CircularVelocity(3, r)

# Total circular velocity
M31_Total_vcirc = M31.CircularVelocityTotal(r)

M31_Hern_vcirc = M31.HernquistVCirc(r, a_guess_M31, M31Total_Halo_mass)

plt.figure(figsize=(8,6))
plt.plot(r, M31_Halo_vcirc, label='Halo')
plt.plot(r, M31_Disk_vcirc, label='Disk')
plt.plot(r, M31_Bulge_vcirc, label='Bulge')
plt.plot(r, M31_Total_vcirc, label='Total')
plt.plot(r, M31_Hern_vcirc, '--', label=f'Hernquist Halo (a={a_guess_M31} kpc)')

plt.title('M31 Rotation Curve')
plt.xlabel('Radius (kpc)')
plt.ylabel('Circular Velocity (km/s)')
plt.legend()
plt.savefig("M31 Rotation Curve.png")


# M33 Rotation Curve
# Circular velocity of each component
M33_Halo_vcirc = M33.CircularVelocity(1, r)
M33_Disk_vcirc = M33.CircularVelocity(2, r)
M33_Bulge_vcirc = M33.CircularVelocity(3, r)

# Total circular velocity
M33_Total_vcirc = M33.CircularVelocityTotal(r)

M33_Hern_vcirc = M33.HernquistVCirc(r, a_guess_M33, M33Total_Halo_mass)

plt.figure(figsize=(8,6))
plt.plot(r, M33_Halo_vcirc, label='Halo')
plt.plot(r, M33_Disk_vcirc, label='Disk')
plt.plot(r, M33_Bulge_vcirc, label='Bulge')
plt.plot(r, M33_Total_vcirc, label='Total')
plt.plot(r, M33_Hern_vcirc, '--', label=f'Hernquist Halo (a={a_guess_M33} kpc)')

plt.title('M33 Rotation Curve')
plt.xlabel('Radius (kpc)')
plt.ylabel('Circular Velocity (km/s)')
plt.legend()
plt.savefig("M33 Rotation Curve.png")