
# # In Class Lab 1
# Must be uploaded to your Github repository under a "Labs/Lab1" by midnight thursday

# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of rest (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 



# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants

#4.74*mu*Ro = VLSR + vsun

def VLSR(Ro, mu=6.379, vsun=12.24*u.km/u.s):
    """
    This funtion will compute the velocity at the local standard of rest
    
    VLSR = 4.74*mu*Ro - vsun
        Reid and Brunthaler 2004
    
    Inputs: Ro (astropy units kpc) Distance from the sun to the galactic center
            mu Is the proper motion of SagA* (mas/yr) 
                Default is from Reid and Brunthaler 2004
            vsun (astropy units km/s) Peculiar motion of the sun in the v direction
                (Schonrich+2010)
                
    Outputs: VLSR (astropy units km/s) The local standard of rest
    """
    
    VLSR = 4.74*mu*(Ro/u.kpc)*u.km/u.s - vsun
    return VLSR

# Different values of the distance to the galactic center
RoReid = 8.34 * u.kpc # Reid + 2014
RoAbuter = 8.178 *u.kpc # Abuter + 2019
RoSparke = 7.9 *u.kpc # Sparke and Gallagher Text

#Compute VLSR using Reid 2014
VLSR_Reid = VLSR(RoReid)
print(VLSR_Reid)

#Compute VLSR using Gravtiy Collab
VLSR_Abuter = VLSR(RoAbuter)
print(np.round(VLSR_Abuter))

#Compute VLSR using Sparke and Gallagher
VLSR_Sparke = VLSR(RoSparke)
print(VLSR_Sparke)

# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr

# Orbital period = 2piR/v
def TorbSun(Ro, Vc):
    """
    Function that computes the orbital period of the sun
    T = 2 pi R / V
    Inputs:
        Ro: (astropy quantity) distance to the galactic center from the sun (kpc)
        Vc: (astropy quantity) velocity of the sun in the "v" direction (km/s)
    Outputs:
        T: (astropy quantity) Orbital period (Gyr)
    """
    VkpcGyr = Vc.to(u.kpc/u.Gyr) # converting V to kpc/Gyr
    T = 2*np.pi*Ro / VkpcGyr # Orbital period
    return T

VsunPec = 12.24*u.km/u.s # Peculiar motion
Vsun = VLSR_Abuter + VsunPec # Total motion of the sun in the "v" direction
    
#Orbital Period of the Sun
T_Abuter = TorbSun(RoAbuter,Vsun)
print(np.round(T_Abuter, 3))

# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)

AgeUniverse = 13.8 * u.Gyr
print(AgeUniverse/T_Abuter)

# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# 
# 
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 


Grav = const.G.to(u.kpc**3 / u.Gyr**2 / u.Msun)

# Density Profile rho = VLSR^2 / (4*pi*G*R^2)
# Mass (r) = Integrate rho dV
#          = Integrate rho 4*pi*r^2*dr
#          = Integrate VLSR^2 / 4*pi*r*2 *4*pi*r^2*dr
#          = Integrate VLSR^2/G dr
#          = VLSR^2/G * r

def massIso(r,VLSR):
    """
    This function will compute the dark matter mass enclosed within a given 
    distance, r, assuming an Isothermal Sphere Model
    M(r) = VLSR^2/G * r
    Inputs:
        r : (astropy quantity) distancce from the galactic center (kpc)
        VLSR : (astropy quantity) the velocity at the Local Standard of Rest (km/s)
    Outputs:
        M: (astropy quantity) mass enclosed within r (Msun)
    """
    
    VLSRkpcGyr = VLSR.to(u.kpc/u.Gyr) # translating to kpc/Gyr
    
    M = VLSRkpcGyr**2/Grav*r # Isothermal Sphere Mass Profile
    
    return M

# Compute mass enclosed within Ro using Gravity Collab
mIsoSolar = massIso(RoAbuter, VLSR_Abuter)
print(mIsoSolar)
print(f"{mIsoSolar:.2e}")

# Compute mass enclosed within 260 kpc
mIso260 = massIso(260*u.kpc, VLSR_Abuter)
print(f"{mIso260:.2e}")

# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)

# Potential for a Hernquist Sphere
# Phi = -G*M/(r+a)

# Escape Speed becomes:
# vesc^2 = 2*G*M/(r+a)

# rearrange for M
# M = vesc^2/2/G*(r+a)

def massHernVesc(vesc, r, a=30*u.kpc):
    """
    This function determines the total dark matter mass needed given an escape speed,
    assuming a Hernquist profile
    M = vesc^2/2/G*(r+a)
    Inputs:
        vesc : (astropy quantity) escape speed (or speed of satellite) (km/s)
        a : (astropy quantity) the Hernquist scale length (kpc)
            default value of 30 kpc
        r : (astropy quantity) distance from the galactic center (kpc)

    Outputs:
        M (astropy quantity) mass within r (Msun)
    """

    vescKpcGyr = vesc.to(u.kpc/u.Gyr) # translate to kpc/Gyr
    
    M = vescKpcGyr**2/2/Grav*(r+a)

    return M

Vleo = 196*u.km/u.s # Speed of Leo I Sohn et al.
r = 260*u.kpc

MLeoI = massHernVesc(Vleo, r)
print(f"{MLeoI:.2e}")
