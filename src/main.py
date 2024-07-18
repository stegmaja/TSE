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
from common import dot,cross

if ic.galactic_tides: 
    from galaxyevolution import yp_galaxy
np.random.seed(ic.seed)

# Default units are Rsun, Msun, Myr
c = c.to(u.Rsun/u.Myr).value
G = G.to(u.Rsun**3/u.Msun/u.Myr**2).value

####################################################################################################
# Useful functions
####################################################################################################

def plot(t,y,m1,m2,m3,logr1,logr2,logr3,filename=str(ic.seed)+'_'+str(ic.Z).replace('.',''),title='Evolution of triple system'):
    print('Plotting...')

    
    # Show only 5000 values
    every = max(1,len(t)//5000)
    t = t[::every]
    y = y[:,::every]
    

    e_in = np.sqrt(y[0]**2+y[1]**2+y[2]**2)
    a_in = y[6]
    e_out = np.sqrt(y[7]**2+y[8]**2+y[9]**2)
    a_out = y[13]

    RL1 = ot.Roche_lobe_radius(m1,m2)
    RL2 = ot.Roche_lobe_radius(m2,m1)
    RL3 = ot.Roche_lobe_radius(m3,m1+m2)

    max_e_in = np.max(e_in)
    max_e_out = np.max(e_out)

    
    m1 = m1[::every]
    m2 = m2[::every]
    m3 = m3[::every]
    logr1 = logr1[::every]
    logr2 = logr2[::every]
    logr3 = logr3[::every]
    RL1 = RL1[::every]
    RL2 = RL2[::every]
    RL3 = RL3[::every]
    

    # Plot in subplots
    fig, axs = plt.subplots(4, 2, figsize=(16,12), sharex=True)

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
    axs[0,1].axhline(1-max_e_in,linestyle='--',color='C0')
    axs[0,1].axhline(1-max_e_out,linestyle='--',color='C1')
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

    # Cos(Angle) between spin and orbital angular momentum
    axs[2,0].plot(t,(y[3]*y[14]+y[4]*y[15]+y[5]*y[16])/np.sqrt(y[3]**2+y[4]**2+y[5]**2)/np.sqrt(y[14]**2+y[15]**2+y[16]**2),label='Primary')
    axs[2,0].plot(t,(y[3]*y[17]+y[4]*y[18]+y[5]*y[19])/np.sqrt(y[3]**2+y[4]**2+y[5]**2)/np.sqrt(y[17]**2+y[18]**2+y[19]**2),label='Secondary')
    axs[2,0].plot(t,(y[14]*y[17]+y[15]*y[18]+y[16]*y[19])/np.sqrt(y[14]**2+y[15]**2+y[16]**2)/np.sqrt(y[17]**2+y[18]**2+y[19]**2),label='Primary-Secondary')
    axs[2,0].set_title('Stellar spin-orbit alignment')
    axs[2,0].set_xlabel('Time [Myr]')
    axs[2,0].set_ylabel('cos(i)')
    axs[2,0].legend(loc='center left')
    axs[2,0].set_ylim(-1,1)
    axs[2,0].set_xlim(0,None)

    # Cos(Angle) between BH spin and orbital angular momentum
    axs[2,1].plot(t,(y[3]*y[20]+y[4]*y[21]+y[5]*y[22])/np.sqrt(y[3]**2+y[4]**2+y[5]**2)/np.sqrt(y[20]**2+y[21]**2+y[22]**2),label='Primary')
    axs[2,1].plot(t,(y[3]*y[23]+y[4]*y[24]+y[5]*y[25])/np.sqrt(y[3]**2+y[4]**2+y[5]**2)/np.sqrt(y[23]**2+y[24]**2+y[25]**2),label='Secondary')
    axs[2,1].plot(t,(y[20]*y[23]+y[21]*y[24]+y[22]*y[25])/np.sqrt(y[20]**2+y[21]**2+y[22]**2)/np.sqrt(y[23]**2+y[24]**2+y[25]**2),label='Primary-Secondary')
    axs[2,1].set_title('BH spin-orbit alignment')
    axs[2,1].set_xlabel('Time [Myr]')
    axs[2,1].set_ylabel('cos(i)')
    axs[2,1].legend(loc='center left')
    axs[2,1].set_ylim(-1,1)
    axs[2,1].set_xlim(0,None)

    # Cos(angle) of inner and outer orbit
    axs[3,0].plot(t,y[5]/np.sqrt(y[3]**2+y[4]**2+y[5]**2),label='Inner')
    axs[3,0].plot(t,y[12]/np.sqrt(y[10]**2+y[11]**2+y[12]**2),label='Outer')
    axs[3,0].set_title('Orbital alignment')
    axs[3,0].set_xlabel('Time [Myr]')
    axs[3,0].set_ylabel('cos(i)')
    axs[3,0].legend(loc='center left')
    axs[3,0].set_ylim(-1,1)
    axs[3,0].set_xlim(0,None)

    # Add title for entire plot
    fig.suptitle(title)

    # Save figure
    plt.tight_layout()
    plt.savefig('./../plots/'+filename+'.png')
    plt.close()

    print('Plot saved as ./../plots/'+filename+'.png')

def store(t,y,star1,star2,star3,filename=str(ic.seed)+'_'+str(ic.Z).replace('.',''),status='Initial'):


    # Stellar types
    k1 = star1.interpolators['k'](t)
    k2 = star2.interpolators['k'](t)
    k3 = star3.interpolators['k'](t)

    # Masses0
    m10 = star1.interpolators['m0'](t)
    m20 = star2.interpolators['m0'](t)
    m30 = star3.interpolators['m0'](t)
    
    # Masses
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)
    m3 = star3.interpolators['m'](t)

    # Radii
    r1 = 10**star1.interpolators['logr'](t)
    r2 = 10**star2.interpolators['logr'](t)
    r3 = 10**star3.interpolators['logr'](t)

    # Epochs
    epoch1 = star1.interpolators['epoch'](t)
    epoch2 = star2.interpolators['epoch'](t)
    epoch3 = star3.interpolators['epoch'](t)

    # opsins
    ospin1 = np.linalg.norm(y[14:17])
    ospin2 = np.linalg.norm(y[17:20])
    ospin3 = np.linalg.norm(y[20:23])

    # Luminosity
    lum1 = 10**star1.interpolators['loglum'](t)
    lum2 = 10**star2.interpolators['loglum'](t)
    lum3 = 10**star3.interpolators['loglum'](t)

    # Envelope mass
    menv1 = star1.interpolators['menv'](t)
    menv2 = star2.interpolators['menv'](t)
    menv3 = star3.interpolators['menv'](t)

    # Envelope radius
    renv1 = star1.interpolators['renv'](t)
    renv2 = star2.interpolators['renv'](t)
    renv3 = star3.interpolators['renv'](t)

    # Core mass
    massc1 = star1.interpolators['massc'](t)
    massc2 = star2.interpolators['massc'](t)
    massc3 = star3.interpolators['massc'](t)

    # Core radius
    radc1 = star1.interpolators['radc'](t)
    radc2 = star2.interpolators['radc'](t)
    radc3 = star3.interpolators['radc'](t)
    
    # Save data
    output = [t,status,y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10],y[11],y[12],y[13],y[14],y[15],y[16],y[17],y[18],y[19],y[20],y[21],y[22],y[23],y[24],y[25],
              k1,k2,k3,m10,m20,m30,m1,m2,m3,r1,r2,r3,epoch1,epoch2,epoch3,ospin1,ospin2,ospin3,lum1,lum2,lum3,menv1,menv2,menv3,renv1,renv2,renv3,massc1,massc2,massc3,radc1,radc2,radc3,
              ic.m1,ic.m2,ic.m3,ic.a1,ic.a2,ic.e1,ic.e2,ic.cosi1,ic.cosi2,ic.omega1,ic.omega2,ic.Omega1,ic.Omega2,ic.method,ic.rtol,ic.atol,ic.max_time,ic.max_step,ic.Z,ic.neta,ic.bwind,ic.hewind,ic.alpha1,ic.lamb,ic.ceflag,ic.tflag,ic.ifflag,ic.wdflag,ic.bhflag,ic.nsflag,ic.piflag,ic.mxns,ic.idum,ic.pts1,ic.pts2,ic.pts3,ic.sigma1,ic.sigma2,ic.beta,ic.xi,ic.acc2,ic.epsnov,ic.eddfac,ic.gamma,ic.seed,ic.galactic_tides,ic.x_MW,ic.y_MW,ic.z_MW,ic.vx_MW,ic.vy_MW,ic.vz_MW,ic.stellar_tides,ic.CE]


    # Append to the end of the file
    with open('./../data/'+filename+'.csv', 'ab') as f:
        f.writelines([b','.join([str(i).encode() for i in output])+b'\n'])

def evolve_remaining_inner_binary(t,y,star1,star2,star3):
    # Print how the remaining binary evolves without the tertiary
    print('Evolve remaining inner binary',end='\n\n')

    # Stellar masses
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)

    # Stellar types
    k1 = star1.interpolators['k'](t)
    k2 = star2.interpolators['k'](t)

    # Inner eccentricity
    e1 = np.linalg.norm(y[0:3])

    # Inner semi-major axis
    a1 = y[6]

    # Spins
    ospin1 = np.linalg.norm(y[14:17])/1e6
    ospin2 = np.linalg.norm(y[17:20])/1e6

    # Epochs
    epoch1 = star1.interpolators['epoch'](t)
    epoch2 = star2.interpolators['epoch'](t)

    # Masses0
    m0_1 = star1.interpolators['m0'](t)
    m0_2 = star2.interpolators['m0'](t)

    P_in = ot.orbital_period(a1,m=m1+m2,units=(u.Rsun,u.day,u.Msun))

    _ = InteractingBinaryStar(mass_1=m1,mass_2=m2,mass0_1=m0_1,mass0_2=m0_2,
                    period=P_in,ecc=e_in,type1=-int(k1),type2=-int(k2),
                    epoch_1=epoch1,epoch_2=epoch2,
                    ospin_1=ospin1,ospin_2=ospin2,
                    tphys=t,max_time=ic.max_time,
                    just_print=True)

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
    result: Stability (zero if unstable)
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
        return 0

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
    result: Supernova (negative if supernova)
    '''

    if ic.primary_SN:
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
    result: Supernova (negative if supernova)
    '''

    if ic.secondary_SN:
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
    result: Supernova (negative if supernova)
    '''

    if ic.tertiary_SN:
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
    result: Unphysical solution (zero if unphysical)
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
        return 0
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

    '''
    Check for an LMXB triple (zero if LMXB triple)
    '''
    '''
    # Stellar types
    k1 = star1.interpolators['k'](t)
    k2 = star2.interpolators['k'](t)
    k3 = star3.interpolators['k'](t)

    # Radii
    r2 = 10**star2.interpolators['logr'](t)

    # Masses
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)
    m3 = star3.interpolators['m'](t)

    # Semi-major axes
    a_in = y[6]
    a_out = y[13]

    # Inner orbital period
    P_in = ot.orbital_period(a_in,m=m1+m2,units=(u.Rsun,u.day,u.Msun))

    # Eccentricities
    e_in = np.linalg.norm(y[0:3])

    # Roche lobe radii
    RL2 = ot.Roche_lobe_radius(m2,m1)

    # Check for a LMXB triple

    if k1 == 14 and k2 < 10 and k3 < 10 and P_in >= 1 and P_in <= 10 and m1 >= 9 and m1 < 11 and a_out > 3e3 and r2 - RL2*a_in*(1-e_in) > 1 and m3>.5 and m3<1.5:
        return 0
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
                                   tphys=t)
    
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
    y_new[14:17] *= binary.ospin1*1e6/np.linalg.norm(y[14:17])
    y_new[17:20] *= binary.ospin2*1e6/np.linalg.norm(y[17:20])
    if binary.k1 == 14:
        y_new[20:23] = y[14:17]/np.linalg.norm(y[14:17])*G*star1.interpolators['m'](t_new)**2/c
    if binary.k2 == 14:
        y_new[23:26] = y[17:20]/np.linalg.norm(y[17:20])*G*star2.interpolators['m'](t_new)**2/c
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

    # BH flags
    if star1.kfinal == 14 and ic.bhflag == 0:
        vkick = 0
    elif star1.kfinal == 14 and ic.bhflag == 1:
        vkick *= (1-star1.ffb)
    elif (star1.kfinal == 14 or star1.kfinal == 13) and ic.bhflag == 2:
        vkick *= (1-star1.ffb)
    elif ic.bhflag == 3:
        raise ValueError('bhflag = 3 not yet implemented.')

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
        event_status = -2
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

    # BH flags
    if star3.kfinal == 14 and ic.bhflag == 0:
        vkick = 0
    elif star3.kfinal == 14 and ic.bhflag == 1:
        vkick *= (1-star3.ffb)
    elif (star3.kfinal == 14 or star3.kfinal == 13) and ic.bhflag == 2:
        vkick *= (1-star3.ffb)
    elif ic.bhflag == 3:
        raise ValueError('bhflag = 3 not yet implemented.')

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

    yp[6] += a*np.abs(dm1)/m12*1e6
    yp[6] += a*np.abs(dm2)/m12*1e6

    yp[13] += A*np.abs(dm1)/m123*1e6
    yp[13] += A*np.abs(dm2)/m123*1e6
    yp[13] += A*np.abs(dm3)/m123*1e6

    ### Spin evolution ###

    if ic.stellar_tides:
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

    # Print how the inner binary would have evolve if it were isolated
    print('Evolve inner binary as if it was isolated',end='\n\n')
    P_in = ot.orbital_period(ic.a1,m=ic.m1+ic.m2,units=(u.Rsun,u.day,u.Msun))
    _ = InteractingBinaryStar(mass_1=ic.m1,mass_2=ic.m2,mass0_1=ic.m1,mass0_2=ic.m2,
                                period=P_in,ecc=ic.e1,type1=1,type2=1,just_print=True)

    # Initial state vector
    y0 = np.zeros(26)

    y0[0:3] = evec_in*e_in
    y0[3:6] = lvec_in*np.sqrt(1-e_in**2)
    y0[6] = a_in

    y0[7:10] = evec_out*e_out
    y0[10:13] = lvec_out*np.sqrt(1-e_out**2)
    y0[13] = a_out

    # Star objects
    print(' ')
    print('Evolve triple system',end='\n\n')
    print('Evolve single stars (ignore entries for dummy secondary)',end='\n\n')
    star1 = SingleStar(mass0_1=m1)
    star2 = SingleStar(mass0_1=m2)
    star3 = SingleStar(mass0_1=m3)

    # Stellar spins
    y0[14:17] = lvec_in*star1.interpolators['ospin'](0)*1e6
    y0[17:20] = lvec_in*star2.interpolators['ospin'](0)*1e6

    # Store initial conditions
    store(t=0,y=y0,star1=star1,star2=star2,star3=star3,status='Initial')

    # Terminating events
    CustomEvent.terminal = True
    primary_RL.terminal = True
    secondary_RL.terminal = True
    tertiary_RL.terminal = True
    unstable.terminal = True
    primary_SN.terminal = True
    secondary_SN.terminal = True
    tertiary_SN.terminal = True
    DCO_merger.terminal = True
    Unphysical.terminal = True
    events = [CustomEvent,primary_RL,secondary_RL,tertiary_RL,unstable,primary_SN,secondary_SN,tertiary_SN,DCO_merger,Unphysical]
    event_label = ['Custom event','Primary Roche lobe overflow','Secondary Roche lobe overflow','Tertiary Roche lobe overflow','Unstable','Primary supernova','Secondary supernova','Tertiary supernova','DCO merger','Unphysical']

    # Do a short burn-in integration
    sol = solve_ivp(evolve, [0,1e-3], y0, args=(star1,star2,star3), method=ic.method, events=events, atol=ic.atol, rtol=ic.rtol, max_step=ic.max_step)

    if sol.status == -1:
        print('Integration failed at burn-in. Likely unrealistic initial conditions')
        store(t=0,y=y0,star1=star1,star2=star2,star3=star3,status='Failed at burn-in')
        sys.exit()
    elif sol.status == 1:
        print('Termination event at burn-in. Likely unrealistic initial conditions')
        print(sol.t_events)
        store(t=0,y=y0,star1=star1,star2=star2,star3=star3,status='Failed at burn-in')
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
                print('Custom event at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Custom event')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Custom event')
                sys.exit()
            elif i == 1:
                print('Primary Roche lobe overflow at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Primary Roche lobe overflow')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Begin primary RLO')
                t_new,y_new,event_status,star1,star2 = model_RLO(sol.t_events[i][0],sol.y[:,-1],star1,star2)
                if event_status == -1:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End primary RLO (no bpp)')
                    sys.exit()
                elif event_status == -2:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End primary RLO (merger)')
                    sys.exit()
                elif event_status == -3:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End primary RLO (not intact)')
                    sys.exit()
                else:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End primary RLO')
            elif i == 2:
                print('Secondary Roche lobe overflow at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Secondary Roche lobe overflow')
                t_new,y_new,event_status,star1,star2 = model_RLO(sol.t_events[i][0],sol.y[:,-1],star1,star2)
                if event_status == -1:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End secondary RLO (no bpp)')
                    sys.exit()
                elif event_status == -2:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End secondary RLO (merger)')
                    sys.exit()
                elif event_status == -3:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End secondary RLO (not intact)')
                    sys.exit()
                else:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='End secondary RLO')
            elif i == 3:
                print('Tertiary Roche lobe overflow at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Tertiary Roche lobe overflow')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Tertiary RLO')
                sys.exit()
            elif i == 4:
                print('Unstable at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Unstable')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Unstable')
                sys.exit()
            elif i == 5:
                ic.primary_SN = True
                print('Primary supernova at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Primary supernova')
                t_new,y_new,event_status = apply_inner_SN(sol.t_events[i][0],sol.y[:,-1],star1,star2,star3)
                if star1.interpolators['k'](t_new) == 14:
                    y_new[20:23] = y_new[14:17]/np.linalg.norm(y_new[14:17])*G*star1.interpolators['m'](t_new)**2/c # BH spins
                if event_status == -1:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Primary SN (inner binary unbound)')
                    sys.exit()
                elif event_status == -2:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Primary SN (outer binary unbound)')
                    evolve_remaining_inner_binary(t_new,y_new,star1,star2,star3)
                    sys.exit()
                else:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Primary SN')
            elif i == 6:
                ic.secondary_SN = True
                print('Secondary supernova at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Secondary supernova')
                t_new,y_new,event_status = apply_inner_SN(sol.t_events[i][0],sol.y[:,-1],star2,star1,star3)
                if star2.interpolators['k'](t_new) == 14:
                    y_new[23:26] = y_new[17:20]/np.linalg.norm(y_new[17:20])*G*star2.interpolators['m'](t_new)**2/c # BH spins
                if event_status == -1:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Secondary SN (inner binary unbound)')
                    sys.exit()
                elif event_status == -2:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Secondary SN (outer binary unbound)')
                    evolve_remaining_inner_binary(t_new,y_new,star1,star2,star3)
                    sys.exit()
                else:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Secondary SN')
            elif i == 7:
                ic.tertiary_SN = True
                print('Tertiary supernova at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Tertiary supernova')
                t_new,y_new,event_status = apply_outer_SN(sol.t_events[i][0],sol.y[:,-1],star1,star2,star3)
                if event_status == -1:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Tertiary SN (outer binary unbound)')
                    evolve_remaining_inner_binary(t_new,y_new,star1,star2,star3)
                    sys.exit()
                else:
                    store(t=t_new,y=y_new,star1=star1,star2=star2,star3=star3,status='Tertiary SN')
            elif i == 8:
                print('DCO merger at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='DCO merger')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='DCO merger')
                sys.exit()
            elif i == 9:
                print('Unphysical solution at',sol.t_events[i][0])
                plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Unphysical solution')
                store(t=sol.t_events[i][0],y=sol.y_events[i][-1],star1=star1,star2=star2,star3=star3,status='Unphysical solution')
                sys.exit()

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
            store(t=sol.t[-1],y=sol.y[:,-1],star1=star1,star2=star2,star3=star3,status='Integration step failed')
            sys.exit()

    print(sol.message)
    plot(t_sol,y_sol,m1_sol,m2_sol,m3_sol,logr1_sol,logr2_sol,logr3_sol,title='Finished integration')
    store(t=t_sol[-1],y=y_sol[:,-1],star1=star1,star2=star2,star3=star3,status='Finished integration')