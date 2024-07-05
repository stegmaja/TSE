import numpy as np
from scipy.integrate import solve_ivp
import sys
import orbittools as ot
from stars import SingleStar, InteractingBinaryStar
from scipy.stats import maxwell
import astropy.units as u
from astropy.constants import G,c
import matplotlib.pyplot as plt
from initialconditions import ic
from spinevolution import yp_spins
from galaxyevolution import Gal,yp_galaxy
from common import dot,cross

np.random.seed(ic.seed)

# Default units are Rsun, Msun, Myr
c = c.to(u.Rsun/u.Myr).value
G = G.to(u.Rsun**3/u.Msun/u.Myr**2).value

####################################################################################################
# Useful functions
####################################################################################################

def plot(t,y,m1,m2,m3,logr1,logr2,logr3,filename=str(ic.seed),title='Evolution of triple system'):

    e_in = np.sqrt(y[0]**2+y[1]**2+y[2]**2)
    a_in = y[6]
    e_out = np.sqrt(y[7]**2+y[8]**2+y[9]**2)
    a_out = y[13]

    RL1 = ot.Roche_lobe_radius(m1,m2)
    RL2 = ot.Roche_lobe_radius(m2,m1)
    RL3 = ot.Roche_lobe_radius(m3,m1+m2)

    # Plot in subplots
    fig, axs = plt.subplots(2, 2, figsize=(12,12))

    # Semi-major axes
    axs[0,0].plot(t,a_in,label='Inner')
    axs[0,0].plot(t,a_out,label='Outer')
    axs[0,0].set_title('Semi-major axes')
    axs[0,0].set_xlabel('Time [Myr]')
    axs[0,0].set_ylabel('a [Rsun]')
    axs[0,0].set_yscale('log')
    axs[0,0].legend(loc='center left')
    axs[0,0].set_xlim(0,None)

    # Eccentricities
    axs[0,1].plot(t,1-e_in,label='Inner')
    axs[0,1].plot(t,1-e_out,label='Outer')
    axs[0,1].set_title('Eccentricities')
    axs[0,1].set_xlabel('Time [Myr]')
    axs[0,1].set_ylabel('1-e')
    axs[0,1].set_yscale('log')
    axs[0,1].legend(loc='center left')
    axs[0,1].set_ylim(None,1)
    axs[0,1].set_xlim(0,None)

    # Radii
    axs[1,1].plot(t,10**logr1, label='Primary')
    axs[1,1].plot(t,10**logr2, label='Secondary')
    axs[1,1].plot(t,10**logr3, label='Tertiary')
    axs[1,1,].plot(t,a_in*(1-e_in)*RL1,linestyle='--',label='Primary RL',color='C0')
    axs[1,1].plot(t,a_in*(1+e_in)*RL2,linestyle='--',label='Secondary RL',color='C1')
    axs[1,1].plot(t,a_out*(1-e_out)*RL3,linestyle='--',label='Tertiary RL',color='C2')
    axs[1,1].set_title('Radii')
    axs[1,1].set_xlabel('Time [Myr]')
    axs[1,1].set_ylabel('R [Rsun]')
    axs[1,1].set_yscale('log')
    axs[1,1].legend(loc='center left')
    axs[1,1].set_xlim(0,None)

    # Masses
    axs[1,0].plot(t,m1, label='Primary')
    axs[1,0].plot(t,m2, label='Secondary')
    axs[1,0].plot(t,m3, label='Tertiary')
    axs[1,0].set_title('Masses')
    axs[1,0].set_xlabel('Time [Myr]')
    axs[1,0].set_ylabel('M [Msun]')
    axs[1,0].legend(loc='center left')
    axs[1,0].set_ylim(0,None)
    axs[1,0].set_xlim(0,None)

    # Add title for entire plot
    fig.suptitle(title)

    # Save figure
    plt.tight_layout()
    plt.savefig('./../plots/'+filename+'.png')
    plt.close()

####################################################################################################
# Functions that detect termination events
####################################################################################################
def primary_RL(t,y,star1,star2,star3):
    '''
    Function to test whether the primary star fills its Roche lobe

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    result: RL-overfilling (positive if RL-filling)
    '''

    # Mass of each star at time t
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)

    # Radius of primary star at time t
    r1 = 10**star1.interpolators['logr'](t)

    # Inner semi-major axis
    a_in = y[6]

    # Inner eccentricity
    e_in = np.linalg.norm(y[0:3])

    # Roche lobe factor
    R_L = ot.Roche_lobe_radius(m1,m2)

    return r1 - R_L*a_in*(1-e_in)

def secondary_RL(t,y,star1,star2,star3):
    '''
    Function to test whether the secondary star fills its Roche lobe

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    result: RL-overfilling (positive if RL-filling)
    '''

    # Mass of each star at time t
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)

    # Radius of secondary star at time t
    r2 = 10**star2.interpolators['logr'](t)

    # Inner semi-major axis
    a_in = y[6]

    # Inner eccentricity
    e_in = np.linalg.norm(y[0:3])

    # Roche lobe factor
    R_L = ot.Roche_lobe_radius(m2,m1)

    return r2 - R_L*a_in*(1-e_in)

def tertiary_RL(t,y,star1,star2,star3):
    '''
    Function to test whether the tertiary companion fills its Roche lobe

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    result: RL-overfilling (positive if RL-filling)
    '''

    # Mass of each star at time t
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)
    m3 = star3.interpolators['m'](t)

    # Radius of tertiary star at time t
    r3 = 10**star3.interpolators['logr'](t)

    # Outer semi-major axis
    a_out = y[13]

    # Outer eccentricity
    e_out = np.linalg.norm(y[7:10])

    # Roche lobe factor
    R_L = ot.Roche_lobe_radius(m3,m1+m2)

    return r3 - R_L*a_out*(1-e_out)

def unstable(t,y,star1,star2,star3):
    '''
    Function to test whether the triple system is unstable

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    result: Stability (positive if unstable)
    '''

    # Mass of each star at time t
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)
    m3 = star3.interpolators['m'](t)

    # Inner semi-major axis
    a_in = y[6]

    # Outer semi-major axis
    a_out = y[13]

    # Outer eccentricity
    e_out = np.linalg.norm(y[7:10])

    if ot.check_triple_stability(a_in,a_out,e_out,m1+m2,m3) == True:
        return -1
    else:
        return 1

def primary_SN(t,y,star1,star2,star3):
    '''
    Function to test whether the primary star undergoes a supernova

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    result: Supernova (positive if supernova)
    '''

    if star1.initial_CO:
        return -np.inf
    else:
        return star1.tco - t

def secondary_SN(t,y,star1,star2,star3):
    '''
    Function to test whether the secondary star undergoes a supernova

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    result: Supernova (positive if supernova)
    '''

    if star2.initial_CO:
        return -np.inf
    else:
        return star2.tco - t

def tertiary_SN(t,y,star1,star2,star3):
    '''
    Function to test whether the tertiary star undergoes a supernova

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    result: Supernova (positive if supernova)
    '''

    if star3.initial_CO:
        return -np.inf
    else:
        return star3.tco - t
    
def DCO_merger(t,y,star1,star2,star3):
    '''
    Function to test whether two compact objects in the inner binary merge

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    result: Merger (positive if merger)
    '''

    # Stellar types
    k1 = star1.interpolators['k'](t)
    k2 = star2.interpolators['k'](t)

    if k1 < 13 and k2 < 13:
        return -np.inf
    else:
        return (1e-3*u.AU).to(u.Rsun).value - y[6]
    
def Unphysical(t,y,star1,star2,star3):
    '''
    Function that test unphysical solutions

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    result: Unphysical solution (positive if unphysical)
    '''

    # Inner semi-major axis
    a_in = y[6]

    # Inner eccentricity
    e_in = np.linalg.norm(y[0:3])

    # Outer semi-major axis
    a_out = y[13]

    # Outer eccentricity
    e_out = np.linalg.norm(y[7:10])

    if e_in >= 1 or e_in < 0 or a_in <= 0 or ~np.isfinite(a_in) or ~np.isfinite(e_in) or ~np.isfinite(a_in).all():
        return 1
    else:
        return -1
    
def CustomEvent(t,y,star1,star2,star3):
    '''
    Function that detects a custom event

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    result: Custom event (positive if event)
    '''

    '''
    # Example function that detects the formation of three BHs

    # Stellar types
    k1 = star1.interpolators['k'](t)
    k2 = star2.interpolators['k'](t)
    k3 = star3.interpolators['k'](t)

    if k1 == 14 and k2 == 14 and k3 == 14:
        return 1
    else:
        return -1
    '''

    return -1

####################################################################################################
# Functions that model termination events
####################################################################################################
def model_RLO(t,y,star1,star2):
    ''' 
    Function to model Roche lobe overflow

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object

    Output:
    t_new: Time at which RLO ends
    y_new: Array of variables at time t_new
    event_status: Status describing the outcome of the RLO
    '''

    # Get masses
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)
    m0_1 = star1.interpolators['m0'](t)
    m0_2 = star2.interpolators['m0'](t)

    # Stellar types
    k1 = star1.interpolators['k'](t)
    k2 = star2.interpolators['k'](t)

    # Epochs 
    epoch1 = star1.interpolators['epoch'](t)
    epoch2 = star2.interpolators['epoch'](t)

    # Spins
    ospin1 = np.linalg.norm(y[14:17])/1e6
    ospin2 = np.linalg.norm(y[17:20])/1e6

    # Inner semi-major axis
    a_in = y[6]

    # Inner eccentricity
    e_in = np.linalg.norm(y[0:3])

    # Inner orbital period
    P_in = ot.orbital_period(a_in,m=m1+m2,units=(u.Rsun,u.day,u.Msun))

    binary = InteractingBinaryStar(mass_1=m1,mass_2=m2,mass0_1=m0_1,mass0_2=m0_2,
                                   period=P_in,ecc=e_in,type1=-int(k1),type2=-int(k2),
                                   epoch_1=epoch1,epoch_2=epoch2,
                                   ospin_1=ospin1,ospin_2=ospin2,
                                   tphys=t,max_time=t+50.)
    
    # Get the time at which RLO ends
    t_new = binary.t
    a_new = binary.sep
    e_new = binary.ecc

    # Prepare output
    y_new = y
    y_new[0:3] *= e_new/np.linalg.norm(y[0:3])
    y_new[3:6] *= np.sqrt(1-e_new**2)/np.linalg.norm(y[3:6])
    y_new[6] = a_new
    y_new[13] *= (m1+m2)/binary.m12
    y_new[14:17] *= binary.spin1*1e6/np.linalg.norm(y[14:17])
    y_new[17:20] *= binary.spin2*1e6/np.linalg.norm(y[17:20])
    if binary.k1 == 14:
        y_new[20:23] = y[14:17]/np.linalg.norm(y[14:17])
    if binary.k2 == 14:
        y_new[23:26] = y[17:20]/np.linalg.norm(y[17:20])
    event_status = binary.event_status
    star1 = binary.star1
    star2 = binary.star2

    return t_new,y_new,event_status,star1,star2

def apply_inner_SN(t,y,star1,star2,star3):
    ''' 
    Function to model supernova of inner star

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    t_new: Time at which SN occurs
    y_new: Array of variables at time t_new
    event_status: Status describing the outcome of the SN
    '''

    # Get masses
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)
    m3 = star3.interpolators['m'](t)
    m_SN = star1.m_SN # Mass of the exploding inner star before SN
    dm_SN = star1.dm_SN # Mass-loss during SN
    m_comp = m2 # Inner binary companion acts as companion

    # Inner semi-major axis
    a_in = y[6]

    # Inner eccentricity
    e_in = np.linalg.norm(y[0:3])

    # Inner angles
    lvec_in = y[3:6]/np.linalg.norm(y[3:6])
    evec_in = y[0:3]/np.linalg.norm(y[0:3])
    cos_i_in,Omega_in,omega_in = ot.vectorial_elements_to_orbital_elements(lvec_in, evec_in)

    # Kick magnitude (km/s)
    vkick = maxwell.rvs(scale=star1.sigma)

    # If BH is formed, lower kick by fallback fraction
    if star1.kfinal == 14:
        vkick *= (1-star1.ffb)

    a_in_new,e_in_new,cos_i_in_new,Omega_in_new,omega_in_new,f_in_new,v_com_in_new = ot.apply_kick_to_orbit(a=a_in,
                                                                                                            vkick=vkick,
                                                                                                            m_SN=m_SN,
                                                                                                            dm_SN=dm_SN,
                                                                                                            m_comp=m_comp,
                                                                                                            e=e_in,
                                                                                                            cos_i=cos_i_in,
                                                                                                            Omega=Omega_in,
                                                                                                            omega=omega_in,
                                                                                                            units=(u.Rsun,u.km/u.s,u.Msun),
                                                                                                            verbose=False)
    
    # Now we need to repeat the same for the outer orbit, using the new com velocity as the kick velocity

    # Update masses
    m_SN += m2
    m_comp = m3

    # Outer semi-major axis
    a_out = y[13]

    # Outer eccentricity
    e_out = np.linalg.norm(y[7:10])

    # Outer angles
    lvec_out = y[10:13]/np.linalg.norm(y[10:13])
    evec_out = y[7:10]/np.linalg.norm(y[7:10])
    cos_i_out,Omega_out,omega_out = ot.vectorial_elements_to_orbital_elements(lvec_out, evec_out)

    # Kick magnitude (km/s)
    vkick = np.linalg.norm(v_com_in_new)
    vkick_phi = np.arctan2(v_com_in_new[1],v_com_in_new[0])
    vkick_theta = np.arccos(v_com_in_new[2]/np.linalg.norm(v_com_in_new))

    a_out_new,e_out_new,cos_i_out_new,Omega_out_new,omega_out_new,f_out_new,v_com_out_new = ot.apply_kick_to_orbit(a=a_out,
                                                                                                                   vkick=vkick,
                                                                                                                   m_SN=m_SN,
                                                                                                                   dm_SN=dm_SN,
                                                                                                                   m_comp=m_comp,
                                                                                                                   vkick_phi=vkick_phi, # Note the change in input
                                                                                                                   vkick_theta=vkick_theta, # Note the change in input
                                                                                                                   e=e_out,
                                                                                                                   cos_i=cos_i_out,
                                                                                                                   Omega=Omega_out,
                                                                                                                   omega=omega_out,
                                                                                                                   units=(u.Rsun,u.km/u.s,u.Msun),
                                                                                                                   verbose=False)
    
    # Prepare output
    t_new = t+1e-3
    y_new = y

    lvec_in_new,evec_in_new,_ = ot.orbital_elements_to_vectorial_elements(cos_i_in_new, Omega_in_new, omega_in_new)
    y_new[0:3] = evec_in_new*e_in_new
    with np.errstate(invalid='ignore'): # Ignore warnings about sqrt of negative numbers
        y_new[3:6] = lvec_in_new*np.sqrt(1-e_in_new**2)
    y_new[6] = a_in_new

    lvec_out_new,evec_out_new,_ = ot.orbital_elements_to_vectorial_elements(cos_i_out_new, Omega_out_new, omega_out_new)
    y_new[7:10] = evec_out_new*e_out_new
    with np.errstate(invalid='ignore'): # Ignore warnings about sqrt of negative numbers
        y_new[10:13] = lvec_out_new*np.sqrt(1-e_out_new**2)
    y_new[13] = a_out_new

    # Print if inner orbit is unbound
    if e_in_new >= 1 or e_in_new < 0 or a_in_new <= 0 or ~np.isfinite(a_in_new) or ~np.isfinite(e_in_new) or ~np.isfinite(cos_i_in_new) or ~np.isfinite(Omega_in_new) or ~np.isfinite(omega_in_new) or ~np.isfinite(f_in_new) or ~np.isfinite(a_in_new).all():
        print('Inner orbit gets unbound.',end='\n\n')
        event_status = -1
    # Print if outer orbit is unbound
    elif e_out_new >= 1 or e_out_new < 0 or a_out_new <= 0 or ~np.isfinite(a_out_new) or ~np.isfinite(e_out_new) or ~np.isfinite(cos_i_out_new) or ~np.isfinite(Omega_out_new) or ~np.isfinite(omega_out_new) or ~np.isfinite(f_out_new) or ~np.isfinite(a_out_new).all():
        print('Outer orbit gets unbound.',end='\n\n')
        event_status = -1
    else:
        print('Both orbits remain bound.',end='\n\n')
        event_status = 1

    return t_new,y_new,event_status

def apply_outer_SN(t,y,star1,star2,star3):
    ''' 
    Function to model supernova of outer star

    Input:
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    t_new: Time at which SN occurs
    y_new: Array of variables at time t_new
    event_status: Status describing the outcome of the SN
    '''

    # Get masses
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)
    m_comp = m1+m2 # Inner binary acts as companion
    m_SN = star3.m_SN # Mass of the exploding tertiary star before SN
    dm_SN = star3.dm_SN # Mass-loss during SN

    # Outer semi-major axis
    a_out = y[13]

    # Outer eccentricity
    e_out = np.linalg.norm(y[7:10])

    # Outer angles
    lvec_out = y[10:13]/np.linalg.norm(y[10:13])
    evec_out = y[7:10]/np.linalg.norm(y[7:10])
    cos_i_out,Omega_out,omega_out = ot.vectorial_elements_to_orbital_elements(lvec_out, evec_out)

    # Kick magnitude (km/s)
    vkick = maxwell.rvs(scale=star3.sigma)

    # If BH is formed, lower kick by fallback fractino
    if star3.kfinal == 14:
        print('Fallback fraction:',star3.ffb)
        vkick *= (1-star3.ffb)

    a_out_new,e_out_new,cos_i_out_new,Omega_out_new,omega_out_new,f_out_new,v_com_out_new = ot.apply_kick_to_orbit(a=a_out, 
                                                                                                                   vkick=vkick, 
                                                                                                                   m_SN=m_SN, 
                                                                                                                   dm_SN=dm_SN, 
                                                                                                                   m_comp=m_comp, 
                                                                                                                   e=e_out, 
                                                                                                                   cos_i=cos_i_out, 
                                                                                                                   Omega=Omega_out, 
                                                                                                                   omega=omega_out, 
                                                                                                                   units=(u.Rsun,u.km/u.s,u.Msun), 
                                                                                                                   verbose=False)
    
    # Prepare output
    t_new = t+1e-3
    y_new = y
    
    lvec_out_new,evec_out_new,_ = ot.orbital_elements_to_vectorial_elements(cos_i_out_new, Omega_out_new, omega_out_new)
    y_new[7:10] = evec_out_new*e_out_new
    with np.errstate(invalid='ignore'): # Ignore warnings about sqrt of negative numbers
        y_new[10:13] = lvec_out_new*np.sqrt(1-e_out_new**2)
    y_new[13] = a_out_new

    # Print if orbit is unbound
    if e_out_new >= 1 or e_out_new < 0 or a_out_new <= 0 or ~np.isfinite(a_out_new) or ~np.isfinite(e_out_new) or ~np.isfinite(cos_i_out_new) or ~np.isfinite(Omega_out_new) or ~np.isfinite(omega_out_new) or ~np.isfinite(f_out_new) or ~np.isfinite(a_out_new).all():
        print('Outer orbit gets unbound.',end='\n\n')
        event_status = -1
    else:
        print('Outer orbit remains bound.',end='\n\n')
        event_status = 1

    return t_new,y_new,event_status

####################################################################################################
# Function that models the triple evolution
####################################################################################################
def evolve(t,y,star1,star2,star3):
    '''
    Function to evolve the system of ODEs

    Input: 
    t: Time
    y: Array of variables
    star1: Star 1 object
    star2: Star 2 object
    star3: Star 3 object

    Output:
    yp: Array of derivatives
    '''

    yp = np.zeros_like(y)

    # Get masses
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)
    m3 = star3.interpolators['m'](t)
    m12 = m1+m2
    m123 = m1+m2+m3
    mu1 = m1*m2/m12
    mu2 = m12*m3/m123

    # Mass loss
    dm1 = star1.interpolators['dm'](t)
    dm2 = star2.interpolators['dm'](t)
    dm3 = star3.interpolators['dm'](t)

    ### Orbital frames ###
    ev = y[0:3]
    e = np.linalg.norm(ev)
    jv = y[3:6]
    j = np.linalg.norm(jv)
    Ev = y[7:10]
    E = np.linalg.norm(Ev)
    Jv = y[10:13]
    J = np.linalg.norm(Jv)

    ### Semi-major axes ###
    a = y[6]
    A = y[13]

    ### Angular momenta ###
    L1 = mu1*np.sqrt(G*m12*a)
    L2 = mu2*np.sqrt(G*m123*A)

    n = np.sqrt(G*m12/a**3)

    tLK = m12/m3/n*(A/a)**3*(1-E**2)**(3/2)

    eps_oct = (m2-m1)/m12*a/A*E/(1-E**2)
    eps_GR = 3.*G*m12**2*A**3*(1.-E**2)**(3/2)/a**4/c**2/m3

    ### LK Quad ###

    #0,1,2
    yp[0:3] += 3./4./tLK*(dot(jv,Jv/J)*cross(ev,Jv/J)+2.*cross(jv,ev)-5*dot(ev,Jv/J)*cross(jv,Jv/J))
    #3,4,5
    yp[3:6] += 3./4./tLK*(dot(jv,Jv/J)*cross(jv,Jv/J)-5*dot(ev,Jv/J)*cross(ev,Jv/J))

    #7,8,9
    yp[7:10] += 3./4./tLK/np.sqrt(1-e**2)*L1/L2*(dot(jv,Jv/J)*cross(Ev,jv)-5.*dot(ev,Jv/J)*cross(Ev,ev)-(1/2-3*e**2+25/2*dot(ev,Jv/J)**2-5/2*dot(jv,Jv/J)**2)*cross(Jv/J,Ev))
    #10,11,12
    yp[10:13] += 3./4./tLK*L1/L2*(dot(jv,Jv/J)*cross(Jv/J,jv)-5.*dot(ev,Jv/J)*cross(Jv/J,ev))

    ### LK Oct ###

    #0,1,2
    yp[0:3] += -75*eps_oct/64/tLK*(cross((2*dot(ev,Jv/J)*dot(jv,Jv/J)*ev+(8/5*e**2-1/5-7*dot(ev,Jv/J)**2)*jv),Ev/E)
    +(2*(dot(ev,Ev/E)*dot(jv,Jv/J)+dot(ev,Jv/J)*dot(jv,Ev/E))*ev)
    +2*cross(2*(dot(jv,Jv/J)*dot(jv,Ev/E)-7*dot(ev,Jv/J)*dot(ev,Ev/E))*jv,Jv/J)
    +16/5*dot(ev,Ev/E)*cross(jv,ev))
    #3,4,5
    yp[3:6] += -75*eps_oct/64/tLK*(2*(dot(ev,Ev/E)*dot(jv,Jv/J)+dot(ev,Jv/J)*dot(jv,Ev/E))*jv
    +2*cross((dot(jv,Ev/E)*dot(jv,Jv/J)-7*dot(ev,Ev/E)*dot(ev,Jv/J))*ev,Jv/J)
    +cross(2*dot(ev,Jv/J)*dot(jv,Jv/J)*jv+(8/5*e**2-1/5-7*dot(ev,Jv/J)**2+dot(jv,Jv/J)**2)*ev,Ev/E))

    #7,8,9
    yp[7:10] += -75*eps_oct/64/tLK/np.sqrt(1-E**2)*L1/L2*(2*cross(dot(ev,Jv/J)*dot(jv,Ev)*Ev/E+dot(jv,Jv/J)*dot(ev,Ev)*Ev/E+(1-E**2)/E*dot(ev,Jv/J)*dot(jv,Jv/J)*Jv/J,jv)
    +cross((2*dot(jv,Ev)*dot(jv,Jv/J)*Ev/E-14*dot(ev,Ev)*dot(ev,Jv/J)*Ev/E+(1-E**2)/E*(8/5*e**2-1/5-7*dot(ev,Jv/J)**2+dot(jv,Jv/J)**2))*Jv/J,ev)
    -cross(2*(1/5-8/5*e**2)*dot(ev,Ev/E)*Ev+14*dot(ev,Jv/J)*dot(jv,Ev/E)*dot(jv,Jv/J)*Ev+7*dot(ev,Ev/E)*(8/5*e**2-1/5-7*dot(ev,Jv/J)**2+dot(jv,Jv/J)**2)*Ev,Jv/J))
    #10,11,12
    yp[10:13] += -75*eps_oct/64/tLK*L1/L2*(2*cross(dot(ev,Jv/J)*dot(jv,Ev/E)*Jv/J+dot(ev,Ev/E)*dot(jv,Jv/J)*Jv/J+dot(ev,Jv/J)*dot(jv,Jv/J)*Ev/E,jv)
    +cross(2*dot(jv,Ev/E)*dot(jv,Jv/J)*Jv/J-14*dot(ev,Ev/E)*dot(ev,Jv/J)*Jv/J
    +(8/5*e**2-1/5-7*dot(ev,Jv/J)**2+dot(jv,Jv/J)**2)*Ev/E,ev))

    ### GW ###

    yp[0:3] += -304.*ev*G**3*m1*m2*m12*(1.+121.*e**2/304.)/15./c**5/a**4/(1.-e**2)**2.5
    yp[3:6] +=  304.*e *G**3*m1*m2*m12*(1.+121.*e**2/304.)/15./c**5/a**4/(1.-e**2)**2.5*e/np.sqrt(1-e**2)*jv/j
    yp[6]   += -64.*G**3*m1*m2*m12*(1.+73.*e**2/24.+37.*e**4/96.)/5./c**5/a**3/(1.-e**2)**3.5

    ### 1PN ###

    yp[0:3] += eps_GR/tLK*cross(jv/j,ev)/(1-e**2)

    ### Mass loss ###

    yp[6] += a*np.abs(dm1)/m12
    yp[6] += a*np.abs(dm2)/m12

    yp[13] += A*np.abs(dm1)/m123
    yp[13] += A*np.abs(dm2)/m123
    yp[13] += A*np.abs(dm3)/m123

    ### Spin evolution ###

    yp += yp_spins(t,y,star1,star2,star3)

    ### Galactic tides ###

    if ic.galactic_tides:
        yp += yp_galaxy(t,y,star1,star2,star3)

    return yp

####################################################################################################
# Main routine
####################################################################################################
if __name__ == '__main__':

    # Initial conditions
    m1 = ic.m1
    m2 = ic.m2
    m3 = ic.m3
    a_in = ic.a1
    a_out = ic.a2
    e_in = ic.e1
    e_out = ic.e2
    cos_i_in = ic.cosi1
    omega_in = ic.omega1
    Omega_in = ic.Omega1
    cos_i_out = ic.cosi2
    omega_out = ic.omega2
    Omega_out = ic.Omega2

    lvec_in,evec_in,nvec_in = ot.orbital_elements_to_vectorial_elements(cos_i_in, Omega_in, omega_in)
    lvec_out,evec_out,nvec_out = ot.orbital_elements_to_vectorial_elements(cos_i_out, Omega_out, omega_out)

    print(ic)

    y0 = np.zeros(26)

    y0[0:3] = evec_in*e_in
    y0[3:6] = lvec_in*np.sqrt(1-e_in**2)
    y0[6] = a_in

    y0[7:10] = evec_out*e_out
    y0[10:13] = lvec_out*np.sqrt(1-e_out**2)
    y0[13] = a_out

    # Star objects
    print('Evolve single stars (ignore entries for dummy secondary)',end='\n\n')
    star1 = SingleStar(mass0_1=m1)
    star2 = SingleStar(mass0_1=m2)
    star3 = SingleStar(mass0_1=m3)

    # Stellar spins
    y0[14:17] = lvec_in*star1.interpolators['ospin'](0)*1e6
    y0[17:20] = lvec_in*star2.interpolators['ospin'](0)*1e6

    # Terminating events
    primary_RL.terminal = True
    secondary_RL.terminal = True
    tertiary_RL.terminal = True
    unstable.terminal = True
    primary_SN.terminal = True
    secondary_SN.terminal = True
    tertiary_SN.terminal = True
    DCO_merger.terminal = True
    Unphysical.terminal = True
    CustomEvent.terminal = True
    events = [primary_RL,secondary_RL,tertiary_RL,unstable,primary_SN,secondary_SN,tertiary_SN,DCO_merger,Unphysical,CustomEvent]

    # Do a short burn-in integration
    sol = solve_ivp(evolve, [0,1e-6], y0, args=(star1,star2,star3), method=ic.method, events=events, atol=ic.atol, rtol=ic.rtol, max_step=ic.max_step)

    if sol.status == -1:
        print('Integration failed at burn-in. Likely unrealistic initial conditions')
        sys.exit()
    elif sol.status == 1:
        print('Termination event at burn-in. Likely unrealistic initial conditions')
        sys.exit()
    else:
        sol.status = 2 # Set dummy status to initiate the while loop below

    t_sol = sol.t
    y_sol = sol.y

    m1_sol = star1.interpolators['m'](t_sol)
    m2_sol = star2.interpolators['m'](t_sol)
    m3_sol = star3.interpolators['m'](t_sol)
    logr1_sol = star1.interpolators['logr'](t_sol)
    logr2_sol = star2.interpolators['logr'](t_sol)
    logr3_sol = star3.interpolators['logr'](t_sol)

    while sol.status != 0:
        # Get the last solution
        t0 = t_sol[-1]
        y0 = y_sol[:,-1]

        # Do a short integration
        sol = solve_ivp(evolve, [t0,ic.max_time], y0, args=(star1,star2,star3), method=ic.method, events=events, atol=ic.atol, rtol=ic.rtol, max_step=ic.max_step)

        # Append the solution
        t_sol = np.concatenate((t_sol,sol.t))
        y_sol = np.concatenate((y_sol,sol.y),axis=1)
        m1_sol = np.append(m1_sol,star1.interpolators['m'](sol.t))
        m2_sol = np.append(m2_sol,star2.interpolators['m'](sol.t))
        m3_sol = np.append(m3_sol,star3.interpolators['m'](sol.t))
        logr1_sol = np.append(logr1_sol,star1.interpolators['logr'](sol.t))
        logr2_sol = np.append(logr2_sol,star2.interpolators['logr'](sol.t))
        logr3_sol = np.append(logr3_sol,star3.interpolators['logr'](sol.t))

        # If integration terminated, check which termination event
        if sol.status == 1:
            # Check which entry of sol.t_events is not empty
            for i in range(len(sol.t_events)):
                if sol.t_events[i].size != 0:
                    break

            print('A termination event occured.',end='\n\n')

            if i == 0:
                print('Primary Roche lobe overflow at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Primary Roche lobe overflow')
                t_new,y_new,event_status,star1,star2 = model_RLO(sol.t_events[i][0],sol.y[:,-1],star1,star2)
            elif i == 1:
                print('Secondary Roche lobe overflow at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Secondary Roche lobe overflow')
                t_new,y_new,event_status,star1,star2 = model_RLO(sol.t_events[i][0],sol.y[:,-1],star1,star2)
            elif i == 2:
                print('Tertiary Roche lobe overflow at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Tertiary Roche lobe overflow')
                sys.exit()
            elif i == 3:
                print('Unstable at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Unstable')
                sys.exit()
            elif i == 4:
                print('Primary supernova at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Primary supernova')
                t_new,y_new,event_status = apply_inner_SN(sol.t_events[i][0],sol.y[:,-1],star1,star2,star3)
                if star1.interpolators['k'](t_new) == 14:
                    y_new[20:23] = y_new[14:17]/np.linalg.norm(y_new[14:17]) # BH spins
            elif i == 5:
                print('Secondary supernova at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Secondary supernova')
                t_new,y_new,event_status = apply_inner_SN(sol.t_events[i][0],sol.y[:,-1],star2,star1,star3)
                if star2.interpolators['k'](t_new) == 14:
                    y_new[23:26] = y_new[17:20]/np.linalg.norm(y_new[17:20]) # BH spins
            elif i == 6:
                print('Tertiary supernova at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Tertiary supernova')
                t_new,y_new,event_status = apply_outer_SN(sol.t_events[i][0],sol.y[:,-1],star1,star2,star3)
            elif i == 7:
                print('DCO merger at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='DCO merger')
                sys.exit()
            elif i == 8:
                print('Unphysical solution at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Unphysical solution')
                sys.exit()
            elif i == 9:
                print('Custom event at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Custom event')
                sys.exit()

            if event_status == -1:
                print('Cannot continue integration after termination event.')
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Cannot continue integration after termination event.')
                sys.exit()
            else:
                t_sol = np.append(t_sol,t_new)
                y_sol = np.append(y_sol,y_new.reshape(-1,1),axis=1)
                m1_sol = np.append(m1_sol,star1.interpolators['m'](t_new))
                m2_sol = np.append(m2_sol,star2.interpolators['m'](t_new))
                m3_sol = np.append(m3_sol,star3.interpolators['m'](t_new))
                logr1_sol = np.append(logr1_sol,star1.interpolators['logr'](t_new))
                logr2_sol = np.append(logr2_sol,star2.interpolators['logr'](t_new))
                logr3_sol = np.append(logr3_sol,star3.interpolators['logr'](t_new))

        elif sol.status == -1:
            print('Integration step failed.')
            plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Integration step failed.')
            sys.exit()

    print(sol.message)
    plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Finished integration')