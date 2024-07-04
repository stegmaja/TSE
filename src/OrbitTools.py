import numpy as np
from astropy.constants import G,c
import astropy.units as u
from scipy.optimize import root
from scipy.integrate import quad

def orbital_elements_to_vectors(a, e, cos_i, Omega, omega, f, m=1, units=(u.AU,u.km/u.s,u.Msun)):
    ''' 
    Input:
        a: semi-major axis (AU)
        e: eccentricity
        cos_i: Cos of inclination
        Omega: longitude of the ascending node
        omega: argument of periapsis
        M: true anomaly
        m: total mass (Msun)

    Output:
        rvec: relative position vector (AU)
        vvec: relative velocity vector (km/s)
    '''

    # Raise error if input is not a number
    if(not all(isinstance(x,(int,float)) for x in [a,e,cos_i,Omega,omega,f,m])):
        raise ValueError('All input parameters must be numbers')

    m *= G.to(units[0]*units[1]**2/units[2]).value

    p = (1-e**2)*a
    r = p/(1+e*np.cos(f))
    u1 = np.array([np.cos(Omega),
                   np.sin(Omega),
                   0])
    u2 = np.array([-cos_i*np.sin(Omega),
                   cos_i*np.cos(Omega),
                   np.sqrt(1-cos_i**2)])
    rvec = r*(u1*np.cos(f+omega)+u2*np.sin(f+omega))
    vvec = np.sqrt(m/p)*(-u1*(e*np.sin(omega)+np.sin(f+omega))+u2*(e*np.cos(omega)+np.cos(f+omega)))

    return rvec,vvec

def orbital_elements_to_vectorial_elements(cos_i, Omega, omega):
    ''' 
    Input:
        a: semi-major axis (AU)
        e: eccentricity
        cos_i: Cos of inclination
        Omega: longitude of the ascending node
        omega: argument of periapsis
        m: total mass (Msun)

    Output:
        lvec: unit vector in the direction of the angular momentum
        evec: unit vector in the direction of the eccentricity vector
        nvec: unit vector in the direction of the line of nodes
    '''

    # Raise error if input is not a number
    if(not all(isinstance(x,(int,float)) for x in [cos_i,Omega,omega])):
        raise ValueError('All input parameters must be numbers')
    
    # Raise error if cos_i is not between -1 and 1
    if(cos_i<-1 or cos_i>1):
        raise ValueError('Cosine of inclination must be between -1 and 1')
    
    sin_i = np.sqrt(1-cos_i**2)

    # From Merritt's book, E.q. (4.59)

    lvec = np.array([sin_i*np.sin(Omega),
                     -sin_i*np.cos(Omega),
                     cos_i])
    
    evec = np.array([np.cos(omega)*np.cos(Omega)-np.sin(omega)*np.sin(Omega)*cos_i,
                     np.cos(omega)*np.sin(Omega)+np.sin(omega)*np.cos(Omega)*cos_i,
                     np.sin(omega)*sin_i])
    
    nvec = np.array([-np.sin(omega)*np.cos(Omega)-np.cos(omega)*np.sin(Omega)*cos_i,
                     -np.sin(omega)*np.sin(Omega)+np.cos(omega)*np.cos(Omega)*cos_i,
                     np.cos(omega)*sin_i])

    return lvec,evec,nvec

def vectorial_elements_to_orbital_elements(lvec, evec):
    ''' 
    Input:
        lvec: unit vector in the direction of the angular momentum
        evec: unit vector in the direction of the eccentricity vector

    Output:
        cos_i: Cos of inclination
        Omega: longitude of the ascending node
        omega: argument of periapsis
    '''

    # Raise error if lvec, evec, nvec are not 3D
    if(len(lvec)!=3 or len(evec)!=3):
        raise ValueError('Input vectors must be 3D')
    
    # Raise error if input is not a numpy array containing numbers
    if(not isinstance(lvec,np.ndarray) or not isinstance(evec,np.ndarray) or not all(isinstance(x,(int,float)) for x in lvec) or not all(isinstance(x,(int,float)) for x in evec)):
        raise ValueError('Input vectors must be numpy arrays containing numbers')

    cos_i = lvec[2]

    n = np.cross([0,0,1],lvec)
    N = np.linalg.norm(n)

    Omega = np.arccos(n[0]/N)
    if(n[1]<0): Omega = 2*np.pi-Omega

    omega = np.arccos(np.dot(n,evec)/(N))
    if(evec[2]<0): omega = 2*np.pi-omega

    return cos_i,Omega,omega

def orbital_angular_momentum(a, e, m1=1, m2=1, units=(u.AU,u.km/u.s,u.Msun)):
    '''
    Input:
        a: semi-major axis (AU)
        e: eccentricity
        m1: mass of the primary (Msun)
        m2: mass of the secondary (Msun)
        units: units of the input/output parameters

    Output:
        L: orbital angular momentum
    '''

    # Raise error if input is not a number
    if(not all(isinstance(x,(int,float)) for x in [a,e,m1,m2])):
        raise ValueError('a,e,m1,m2 must be numbers')
    
    # Raise error if a is negative
    if(a<0):
        raise ValueError('a must be positive')
    
    # Raise error if m1, m2 are negative or if e is not between 0 and 1
    if(m1<0 or m2<0 or e<0 or e>=1):
        raise ValueError('m1,m2 must be positive and e must be between 0 and 1')
    
    m1 *= G.to(units[0]*units[1]**2/units[2]).value
    m2 *= G.to(units[0]*units[1]**2/units[2]).value

    m = m1+m2
    mu = m1*m2/m

    L = mu*np.sqrt(m*a*(1-e**2))

    return L
'''
def merger_time(a, e, m1=1, m2=1, F='Numerical_Integration', units=(u.AU,u.yr,u.Msun)):
    
    a *= units[0]
    m1 *= units[2]
    m2 *= units[2]

    m = m1+m2
    mu = m1*m2/m

    g = lambda e : e**(12/19)/(1-e**2)*(1+(121/304)*e**2)**(870/2299)

    if e==0:
        t = 5/256*(c**5/G**3)*a**4/mu/m**2
        return t.to(units[1]).value
    
    elif F=='Numerical_Integration':
        F0 = 48/19/g(e)**4*quad(lambda e : g(e)**4*(1-e**2)**(5/2)/e/(1+121/304*e**2),0,e)[0]

    elif F=='Low_eccentricity':
        F0 = e**(48/19)/g(e)**4

    elif F=='High_eccentricity':
        F0 = 768/429*(1-e**2)**(7/2)

    else:
        raise ValueError('F must be Numerical_Integration, Low_eccentricity or High_eccentricity')
    
    T0 = orbital_period(a,m=m,units=units)

    return t.to(units[1]).value
'''
def vectors_to_orbital_elements(rvec, vvec, m=1, units=(u.AU,u.km/u.s,u.Msun)):
    ''' 
    Input:
        rvec: relative position vector (AU)
        vvec: relative velocity vector (km/s)
        m: total mass (Msun)

    Output:
        a: semi-major axis (AU)
        e: eccentricity
        cos_i: Cos of inclination
        Omega: longitude of the ascending node
        omega: argument of periapsis
        f: true anomaly
    '''

    # Raise error if vectors are not 3D
    if(len(rvec)!=3 or len(vvec)!=3):
        raise ValueError('Input vectors must be 3D')
    
    # Raise error if input is not a numpy array containing numbers or mass is not a number
    if(not isinstance(rvec,np.ndarray) or not isinstance(vvec,np.ndarray) or not isinstance(m,(int,float))):
        raise ValueError('Input vectors must be numpy arrays and mass must be a number')

    m *= G.to(units[0]*units[1]**2/units[2]).value

    r = np.linalg.norm(rvec)
    v = np.linalg.norm(vvec)
    h = np.cross(rvec,vvec)
    H = np.linalg.norm(h)
    n = np.cross([0,0,1],h)
    N = np.linalg.norm(n)

    evec = np.cross(vvec,h)/m - rvec/r

    e = np.linalg.norm(evec)

    a = 1/(2/r-v**2/m)

    cos_i = h[2]/H

    Omega = np.arccos(n[0]/N)
    if(n[1]<0): Omega = 2*np.pi-Omega

    f = np.arccos(np.dot(rvec,evec)/(r*e))
    if(np.dot(rvec,vvec)<0): f = 2*np.pi-f

    omega = np.arccos(np.dot(n,evec)/(N*e))
    if(evec[2]<0): omega = 2*np.pi-omega

    return a,e,cos_i,Omega,omega,f

def get_true_anomaly(e,M=None):
    '''
    Input:
        e: eccentricity
        M: mean anomaly (if not proveded, a random value is chosen between 0 and 2*pi)

    Output:
        f: true anomaly
    '''

    # Raise error if input is not a number
    if(not isinstance(e,(int,float))):
        raise ValueError('Eccentricity must be a number')
    
    # Raise error if e is not between 0 and 1
    if(e<0 or e>=1):
        raise ValueError('Eccentricity must be between 0 and 1')
    
    # Raise error if M is not None or a number
    if(M is not None and not isinstance(M,(int,float))):
        raise ValueError('Mean anomaly must be None or a number')

    # Calculate true anomaly from random mean anomaly
    if M is None:
        M = np.random.uniform(0,2*np.pi)
    sol = root(lambda E : E-e*np.sin(E)-M, np.random.uniform(0,2*np.pi))
    E = sol.x[0]
    beta = e/(1+np.sqrt(1-e**2))
    f = E+2*np.arctan(beta*np.sin(E)/(1-beta*np.cos(E)))

    return f

def apply_kick_to_orbit(a, vkick, m_SN, dm_SN, m_comp, dm_comp=0, vkick_phi=None, vkick_theta=None, e=0, cos_i=1, Omega=0, omega=0, f=None, v_com=np.zeros(3), units=(u.AU,u.km/u.s,u.Msun), verbose=False):
    '''
    Input:
        a: semi-major axis (AU)
        vkick: kick velocity (km/s)
        m_SN: mass of the exploding star before supernova (Msun)
        dm_SN: mass loss of the supernova (Msun)
        m_comp: mass of the companion (Msun)
        dm_comp: mass loss of the companion, e.g., through winds (Msun)
        vkick_phi: azimuthal angle of the kick velocity (rad)
        vkick_theta: polar angle of the kick velocity (rad)
        e: eccentricity
        cos_i: Cos of inclination
        Omega: longitude of the ascending node
        omega: argument of periapsis
        f: true anomaly
        v_com: centre of mass velocity (km/s)
        units: units of the input parameters
        verbose: print output

    Output:
        a_new: semi-major axis (AU)
        e_new: eccentricity
        cos_i_new: Cos of inclination
        Omega_new: longitude of the ascending node
        omega_new: argument of periapsis
        f_new: true anomaly
        v_com_new: centre of mass velocity (km/s)   
    '''

    '''
    # Raise error if input is not a number or None
    if(not all(isinstance(x,(int,float)) or x is None for x in [f,vkick_phi,vkick_theta])):
        raise ValueError('f,vkick_phi,vkick_theta must be numbers or None')
    
    # Raise error if input is not a number
    if(not all(isinstance(x,(int,float)) for x in [a,vkick,m_SN,dm_SN,m_comp,dm_comp])):
        raise ValueError('a,vkick,m_SN,dm_SN,m_comp,dm_comp must be numbers')
    
    # Raise error if verbose is not a boolean   
    if(not isinstance(verbose,bool)):
        raise ValueError('verbose must be a boolean')
    
    # Raise error if v_com is not a numpy array containing numbers
    if(not isinstance(v_com,np.ndarray) or not all(isinstance(x,(int,float)) for x in v_com)):
        raise ValueError('v_com must be a numpy array containing numbers')
    
    # Raise error if a, vkick, m_SN, dm_SN, m_comp, dm_comp are negative or if e is not between 0 and 1
    if(a<0 or vkick<0 or m_SN<0 or dm_SN<0 or m_comp<0 or dm_comp<0 or e<0 or e>=1):
        raise ValueError('a,vkick,m_SN,dm_SN,m_comp,dm_comp must be positive and e must be between 0 and 1')
    
    # Raise error if dm_SN, dm_comp are greater than m_SN, m_comp
    if(dm_SN>m_SN or dm_comp>m_comp):
        raise ValueError('dm_SN,dm_comp must be less than m_SN,m_comp')
    '''

    if f is None:
        f = get_true_anomaly(e)

        if verbose:
            print('True anomaly not provided. Calculated from random mean anomaly.',end='\n\n')

    if vkick_phi is None:
        vkick_phi = np.random.uniform(0,2*np.pi)

        if verbose:
            print('Azimuthal angle of the kick velocity not provided. Randomly chosen:',vkick_phi,end='\n\n')

    if vkick_theta is None:
        vkick_theta = np.arccos(np.random.uniform(-1,1))

        if verbose:
            print('Polar angle of the kick velocity not provided. Randomly chosen:',vkick_theta,end='\n\n')

    # Orbital vectors before supernova
    rvec_old,vvec_old = orbital_elements_to_vectors(a, e, cos_i, Omega, omega, f, m=m_SN+m_comp, units=units)

    # Velocity vectors of each star before supernova
    vvec1_old = vvec_old*m_comp/(m_SN+m_comp) + v_com
    vvec2_old = -vvec_old*m_SN/(m_SN+m_comp) + v_com

    # Apply kick to component 1
    vkick_vec = np.array([vkick*np.sin(vkick_theta)*np.cos(vkick_phi),
                          vkick*np.sin(vkick_theta)*np.sin(vkick_phi),
                          vkick*np.cos(vkick_theta)])
    vvec1_new = vvec1_old + vkick_vec 

    # Calculate new relative velocity
    vvec_new = vvec1_new - vvec2_old

    # Calculate new centre of mass velocity
    v_com_new = ((m_SN-dm_SN)*vvec1_new + (m_comp-dm_comp)*vvec2_old)/(m_SN+m_comp-dm_SN-dm_comp)

    # Calculate new orbital elements
    a_new,e_new,cos_i_new,Omega_new,omega_new,f_new = vectors_to_orbital_elements(rvec_old, vvec_new, m=m_SN+m_comp-dm_SN-dm_comp, units=units)

    # Print output
    if verbose:
        print('Old masses:',m_SN,m_comp)
        print('Old semi-major axis:',a)
        print('Old eccentricity:',e)
        print('Old inclination:',np.arccos(cos_i))
        print('Old relative velocity:',np.linalg.norm(vvec_old))
        print('Old com velocity:',np.linalg.norm(v_com),end='\n\n')

        print('Kick velocity:',vkick,end='\n\n')

        # Print if orbit is unbound
        if e_new >= 1 or e_new < 0 or a_new <= 0 or ~np.isfinite(a_new) or ~np.isfinite(e_new) or ~np.isfinite(cos_i_new) or ~np.isfinite(Omega_new) or ~np.isfinite(omega_new) or ~np.isfinite(f_new) or ~np.isfinite(v_com_new).all():
            print('Orbit gets unbound.',end='\n\n')
        else:
            print('Orbit remains bound.',end='\n\n')

        print('New masses:',m_SN-dm_SN,m_comp-dm_comp)
        print('New semi-major axis:',a_new)
        print('New eccentricity:',e_new)
        print('New inclination:',np.arccos(cos_i_new))
        print('New relative velocity:',np.linalg.norm(vvec_new))
        print('New com velocity:',np.linalg.norm(v_com_new),end='\n\n')

        print('Units:',units[0],',',units[1],',',units[2],end='\n\n')

    return a_new,e_new,cos_i_new,Omega_new,omega_new,f_new,v_com_new

def check_triple_stability(a_in,a_out,e_out,m_in,m_out):
    '''
    Input:
        a_in: semi-major axis of inner binary (AU)
        a_out: semi-major axis of outer binary (AU)
        e_out: eccentricity of outer binary
        m_in: mass of inner binary (Msun)
        m_out: mass of tertiary companion (Msun)
    
    Output:
        stable: True if triple is stable, False otherwise
    '''

    # Raise error if input is not a number
    #if(not all(isinstance(x,(int,float)) for x in [a_in,a_out,e_out,m_in,m_out])):
    #    raise ValueError('All input parameters must be numbers')
    
    # Raise error if input is not a number
    #if(a_in<0 or a_out<0 or m_in<0 or m_out<0 or e_out<0 or e_out>=1):
    #    raise ValueError('a_in,a_out,m_in,m_out must be positive and e_out must be between 0 and 1')

    # Check if inner binary is stable
    stable = a_out/a_in > 2.8/(1-e_out)*((m_in+m_out)/m_in*(1+e_out)/np.sqrt(1-e_out))**(2/5)

    return stable

def orbital_period(a,m=1,units=(u.AU,u.yr,u.Msun)):
    '''
    Input:
        a: semi-major axis (AU)
        m: total mass (Msun)
        units: units of the input/output parameters

    Output:
        T: orbital period (yr)
    '''

    # Raise error if input is not a number
    if(not all(isinstance(x,(int,float)) for x in [a,m])):
        raise ValueError('a,m must be numbers')
    
    # Raise error if a is negative
    if(a<0):
        raise ValueError('a must be positive')
    
    # Raise error if m is negative
    if(m<0):
        raise ValueError('m must be positive')

    m *= G.to(units[0]**3/units[1]**2/units[2]).value

    T = 2*np.pi*np.sqrt(a**3/m)

    return T

def semi_major_axis(T,m=1,units=(u.AU,u.yr,u.Msun)):
    '''
    Input:
        T: orbital period (yr)
        m: total mass (Msun)
        units: units of the input/output parameters

    Output:
        a: semi-major axis (AU)
    '''

    # Raise error if input is not a number
    if(not all(isinstance(x,(int,float)) for x in [T,m])):
        raise ValueError('T,m must be numbers')
    
    # Raise error if T is negative
    if(T<0):
        raise ValueError('T must be positive')
    
    # Raise error if m is negative
    if(m<0):
        raise ValueError('m must be positive')

    m *= G.to(units[0]**3/units[1]**2/units[2]).value

    a = (T**2/(4*np.pi**2)*m)**(1/3)

    return a

def Kozai_timescale(a_in,a_out,e_out,m_in,m_out,units=(u.AU,u.yr,u.Msun)):
    '''
    Input:
        a_in: semi-major axis of inner binary (AU)
        a_out: semi-major axis of outer binary (AU)
        e_out: eccentricity of outer binary
        m_in: mass of inner binary (Msun)
        m_out: mass of tertiary companion (Msun)
        units: units of the input/output parameters

    Output:
        T_Kozai: Kozai timescale (yr)
    '''

    # Raise error if input is not a number
    if(not all(isinstance(x,(int,float)) for x in [a_in,a_out,e_out,m_in,m_out])):
        raise ValueError('All input parameters must be numbers')
    
    # Raise error if a_in, a_out, m_in, m_out are negative or if e_out is not between 0 and 1
    if(a_in<0 or a_out<0 or m_in<0 or m_out<0 or e_out<0 or e_out>=1):
        raise ValueError('a_in,a_out,m_in,m_out must be positive and e_out must be between 0 and 1')
    
    m_in *= G.to(units[0]**3/units[1]**2/units[2]).value
    m_out *= G.to(units[0]**3/units[1]**2/units[2]).value

    # Calculate mean motion
    n = np.sqrt(m_in/a_in**3)

    # Calculate Kozai timescale
    T_Kozai = m_in/n/m_out*(a_out/a_in)**3*(1-e_out**2)**(3/2)

    return T_Kozai

def Roche_lobe_radius(m1,m2):
    '''
    Input:
        m1: mass of the primary
        m2: mass of the secondary
    Output:
        R_L: Roche lobe radius (units of separation)
    '''

    # Raise error if input is not a number
    #if(not all(isinstance(x,(int,float)) for x in [m1,m2])):
    #    print(m1.dtype,m2.dtype)
    #    raise ValueError('m1,m2 must be numbers')
    
    # Raise error if m1, m2 are negative
    #if(m1<0 or m2<0):
    #    raise ValueError('m1,m2 must be positive')
    
    q = m1/m2

    R_L = 0.49*q**(2/3)/(0.6*q**(2/3)+np.log(1+q**(1/3)))

    return R_L

def get_triple_vectors(a_in=None, e_in=None, cos_i_in=None, Omega_in=None, omega_in=None, f_in=None,
                       a_out=None, e_out=None, cos_i_out=None, Omega_out=None, omega_out=None, f_out=None, 
                       log_a_in_min=0, log_a_in_max=2, log_a_out_min=0, log_a_out_max=4,
                       e_in_alpha=1, e_out_alpha=1,
                       rcom=np.zeros(3), vcom=np.zeros(3),
                       check_stable=True, m1=1, m2=1, m3=1, 
                       units=(u.AU,u.km/u.s,u.Msun),
                       verbose=False):
    '''
    Input:
        a_in: semi-major axis of inner binary (AU)
        e_in: eccentricity of inner binary
        cos_i_in: Cos of inclination of inner binary
        Omega_in: longitude of the ascending node of inner binary
        omega_in: argument of periapsis of inner binary
        f_in: true anomaly of inner binary
        a_out: semi-major axis of outer binary (AU)
        e_out: eccentricity of outer binary
        cos_i_out: Cos of inclination of outer binary
        Omega_out: longitude of the ascending node of outer binary
        omega_out: argument of periapsis of outer binary
        f_out: true anomaly of outer binary
        log_a_in_min: minimum log of semi-major axis of inner binary
        log_a_in_max: maximum log of semi-major axis of inner binary
        log_a_out_min: minimum log of semi-major axis of outer binary
        log_a_out_max: maximum log of semi-major axis of outer binary
        e_in_alpha: power-law index for inner binary eccentricity
        e_out_alpha: power-law index for outer binary eccentricity
        rcom: centre of mass position (AU)
        vcom: centre of mass velocity (km/s)
        check_stable: check if triple is stable
        m1: mass of the primary (Msun)
        m2: mass of the secondary (Msun)
        m3: mass of the tertiary (Msun)
        units: units of the input/output parameters

    Output:
        xvec1: position vector of the primary (AU)
        xvec2: position vector of the secondary (AU)
        xvec3: position vector of the tertiary (AU)
        vvec1: velocity vector of the primary (km/s)
        vvec2: velocity vector of the secondary (km/s)
        vvec3: velocity vector of the tertiary (km/s)
    '''

    
    # Raise error if input is not a number or None
    if(not all(isinstance(x,(int,float)) or x is None for x in [a_in,e_in,cos_i_in,Omega_in,omega_in,f_in,a_out,e_out,cos_i_out,Omega_out,omega_out,f_out])):
        raise ValueError('All input parameters must be numbers or None')
    
    # Raise error if units is not a tuple of 3 elements
    if(not isinstance(units,tuple) or len(units)!=3):
        raise ValueError('units must be a tuple of 3 elements')
    
    # Raise error if rcom, vcom are not numpy arrays containing numbers
    if(not isinstance(rcom,np.ndarray) or not isinstance(vcom,np.ndarray) or not all(isinstance(x,(int,float)) for x in rcom) or not all(isinstance(x,(int,float)) for x in vcom)):
        raise ValueError('rcom,vcom must be numpy arrays containing numbers')
    
    # Raise error if check_stable is not a boolean
    if(not isinstance(check_stable,bool)):
        raise ValueError('check_stable must be a boolean')
    
    # Raise error if log_a_in_min, log_a_in_max, log_a_out_min, log_a_out_max, e_in_alpha, e_out_alpha are not numbers
    if(not all(isinstance(x,(int,float)) for x in [log_a_in_min,log_a_in_max,log_a_out_min,log_a_out_max,e_in_alpha,e_out_alpha])):
        raise ValueError('log_a_in_min,log_a_in_max,log_a_out_min,log_a_out_max,e_in_alpha,e_out_alpha must be numbers')
    
    # Raise error if log_a_out_max is smaller than log_a_in_min
    if(log_a_out_max<log_a_in_min):
        raise ValueError('log_a_out_max must be greater than log_a_in_min')
    
    # Total masses
    m_in = m1+m2
    m_out = m_in+m3
    
    # Generate random inner binary parameters if not provided
    if a_in is None:
        a_in = 10**np.random.uniform(log_a_in_min,log_a_in_max)

    if e_in is None:
        e_in = np.random.uniform(0,1)**(1/(e_in_alpha+1)) # Power-law distribution

    if cos_i_in is None:
        cos_i_in = np.random.uniform(-1,1)

    if Omega_in is None:
        Omega_in = np.random.uniform(0,2*np.pi)

    if omega_in is None:
        omega_in = np.random.uniform(0,2*np.pi)

    if f_in is None:
        f_in = get_true_anomaly(e_in)

    # Record which of the outer binary parameters were None and must not be changed
    a_out_was_None = a_out is None
    e_out_was_None = e_out is None
    cos_i_out_was_None = cos_i_out is None
    Omega_out_was_None = Omega_out is None
    omega_out_was_None = omega_out is None
    f_out_was_None = f_out is None  

    if check_stable and not a_out_was_None and e_out_was_None:
        stable = check_triple_stability(a_in,a_out,e_out,m_in,m3)
        raise ValueError('Provided outer binary yield no stable triple system') # Note that one could have resampled a_in
    
    else:
        stable = False

        counter = 0

        while not stable and counter<1e5:

            # Generate random outer binary parameters if not provided
            if a_out_was_None:
                a_out = 10**np.random.uniform(log_a_out_min,log_a_out_max)
                
            if e_out_was_None:
                e_out = np.random.uniform(0,1)**(1/(e_out_alpha+1))

            if cos_i_out_was_None:
                cos_i_out = np.random.uniform(-1,1)

            if Omega_out_was_None:
                Omega_out = np.random.uniform(0,2*np.pi)

            if omega_out_was_None:
                omega_out = np.random.uniform(0,2*np.pi)

            if f_out_was_None:
                f_out = get_true_anomaly(e_out)

            if check_stable:
                stable = check_triple_stability(a_in,a_out,e_out,m_in,m3)
            else:
                stable = True

            counter += 1

        if counter==1e5:
            raise ValueError('No stable triple system found after 1e5 iterations')

    # If verbose, print the generated parameters
    if verbose:
        print('Masses:')
        print('Primary:',m1,units[2])
        print('Secondary:',m2,units[2])
        print('Tertiary:',m3,units[2],end='\n\n')

        print('Generated inner binary parameters:')
        print('Semi-major axis:',a_in,units[0])
        print('Eccentricity:',e_in)
        print('Cos Inclination:',cos_i_in)
        print('Longitude of the ascending node:',Omega_in)
        print('Argument of periapsis:',omega_in)
        print('True anomaly:',f_in,end='\n\n')

        print('Generated outer binary parameters:')
        print('Semi-major axis:',a_out,units[0])
        print('Eccentricity:',e_out)
        print('Cos Inclination:',cos_i_out)
        print('Longitude of the ascending node:',Omega_out)
        print('Argument of periapsis:',omega_out)
        print('True anomaly:',f_out,end='\n\n')

    # Calculate relative vectors
    rvec_in,vvec_in = orbital_elements_to_vectors(a_in, e_in, cos_i_in, Omega_in, omega_in, f_in, m=m_in, units=units)
    rvec_out,vvec_out = orbital_elements_to_vectors(a_out, e_out, cos_i_out, Omega_out, omega_out, f_out, m=m_out, units=units)

    # Calculate positions and velocities
    xvec_in_com = m3/m_out*rvec_out # Centre of mass position of inner binary
    xvec3 = -m_in/m_out*rvec_out # Position of tertiary
    
    vvec_in_com = m3/m_out*vvec_out # Centre of mass velocity of inner binary
    vvec3 = -m_in/m_out*vvec_out # Position of tertiary

    xvec1 = xvec_in_com + m2/m_in*rvec_in # Position of primary
    xvec2 = xvec_in_com - m1/m_in*rvec_in # Position of secondary

    vvec1 = vvec_in_com + m2/m_in*vvec_in # Velocity of primary
    vvec2 = vvec_in_com - m1/m_in*vvec_in # Velocity of secondary

    # Check if centre of mass is at origin and at rest
    if not np.allclose(np.sum([m1*xvec1,m2*xvec2,m3*xvec3],axis=0),np.zeros(3)):
        #raise ValueError('Centre of mass is not at origin')
        print('Centre of mass is not at origin')
    
    if not np.allclose(np.sum([m1*vvec1,m2*vvec2,m3*vvec3],axis=0),np.zeros(3)):
        #raise ValueError('Centre of mass is not at rest')
        print('Centre of mass is not at rest')
    
    # Add centre of mass position and velocity
    xvec1 += rcom
    xvec2 += rcom
    xvec3 += rcom

    vvec1 += vcom
    vvec2 += vcom
    vvec3 += vcom

    return xvec1,xvec2,xvec3,vvec1,vvec2,vvec3

def power_law_sample(alpha, xmin, xmax):
    '''
    Input:
        alpha: power-law index
        xmin: minimum value
        xmax: maximum value
    Output:
        x: random value
    '''

    # Raise error if input is not a number
    if(not all(isinstance(x,(int,float)) for x in [alpha,xmin,xmax])):
        raise ValueError('All input parameters must be numbers')
    
    # Raise error if xmin is greater than xmax
    if(xmin>xmax):
        raise ValueError('xmin must be less than xmax')
    
    # Raise error if alpha is -1
    if(alpha==-1):
        raise ValueError('alpha must not be -1')
    
    stamm_func = lambda x: x**(alpha+1)/(alpha+1)

    Norm = 1/(stamm_func(xmax) - stamm_func(xmin))
    
    return ((alpha+1)*np.random.uniform()/Norm+xmin**(alpha+1))**(1/(alpha+1))

print('OrbitTools.py loaded.',end='\n\n')

