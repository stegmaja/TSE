import numpy as np
import astropy.units as u
from astropy.constants import G,c
import orbittools as ot
from common import dot,cross,f1,f2,f3,f4,f5

# Default units are Rsun, Msun, Myr
G = G.to(u.Rsun**3/u.Msun/u.Myr**2).value
c = c.to(u.Rsun/u.Myr).value

def yp_spins(t,y,star1,star2,star3):
    yp = np.zeros_like(y)

    ev = y[0:3]
    e = np.linalg.norm(ev)
    jv = y[3:6]
    j = np.linalg.norm(jv)
    qv = cross(jv/j,ev/e) # Is normalized

    # Get masses
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)
    m12 = m1+m2
    mu1 = m1*m2/m12

    # Get radii
    r1 = 10**star1.interpolators['logr'](t)
    r2 = 10**star2.interpolators['logr'](t)

    # Get stellar types
    k1 = star1.interpolators['k'](t)
    k2 = star2.interpolators['k'](t)

    # Get luminosities
    lum1 = 10**star1.interpolators['loglum'](t)
    lum2 = 10**star2.interpolators['loglum'](t)

    # Envelope masses
    menv1 = star1.interpolators['menv'](t)
    menv2 = star2.interpolators['menv'](t)

    # Envelope radii
    renv1 = star1.interpolators['renv'](t)
    renv2 = star2.interpolators['renv'](t)

    # Core masses
    massc1 = star1.interpolators['massc'](t)
    massc2 = star2.interpolators['massc'](t)

    # Core radii
    radc1 = star1.interpolators['radc'](t)
    radc2 = star2.interpolators['radc'](t)

    # Spins
    ospin1 = np.linalg.norm(y[14:17])
    ospin2 = np.linalg.norm(y[17:20])

    # Spin vectors
    Ospin1 = y[14:17]
    Ospin2 = y[17:20]

    # Make sure that envelope radii are not larger than stellar radii
    renv1 = min(renv1,r1)
    renv2 = min(renv2,r2)

    # Make sure that envelope masses are not larger than stellar masses
    menv1 = min(menv1,m1)
    menv2 = min(menv2,m2)

    # Arrays
    m_prim = np.array([m1,m2])
    r_prim = np.array([r1,r2])
    k_prim = np.array([k1,k2])
    m_comp = np.flip(m_prim)
    lum_prim = np.array([lum1,lum2])
    menv_prim = np.array([menv1,menv2])
    renv_prim = np.array([renv1,renv2])
    ospin_prim = np.array([ospin1,ospin2])
    Ospin_prim = np.array([Ospin1,Ospin2])

    # Semi-major axis
    a = y[6]

    # Orbital period
    P_in = ot.orbital_period(a,m=m1+m2,units=(u.Rsun,u.Myr,u.Msun))
    n = 2*np.pi/P_in

    # Angular momentum
    h1 = mu1*np.sqrt(G*m12*a)*np.sqrt(1-e**2)

    for i in range(2):

        # Apsidal motion constant
        kAM = 0.014 # Could be changed based on stellar type

        # Determine which tides
        if (k_prim[i] == 1 and m_prim[i] > 1.2) or k_prim[i] == 4 or k_prim[i] == 7:
            # Radiative tides
            E2 = 1.592e-9*m_prim[i]**2.84
            kAM_over_T = E2*(1+m_comp[i]/m_prim[i])**(5/6)*r_prim[i]*np.sqrt(G*m_prim[i]/a**5)
        elif k_prim[i] < 10:
            # Convective tides
            mr23yr = 0.4311 # Prefactor from mobse to convert tau_conv to yr
            mr23Myr = mr23yr*1e6 # Convert to Myr
            tau_conv = mr23Myr*(menv_prim[i]*renv_prim[i]*(r_prim[i]-1/2*renv_prim[i])/3/lum_prim[i])**(1/3)
            Ptid = 1/abs(1/P_in-ospin_prim[i]/2/np.pi)
            f_conv = min(1,(Ptid/2/tau_conv)**2)
            kAM_over_T = 2/21*f_conv/tau_conv*menv_prim[i]/m_prim[i]
        else:
            # No tides
            continue

        # Viscous timescale
        tV = 3*(1+kAM)**2/kAM_over_T

        # Tidal timescale
        tf = tV/9*(a/r_prim[i])**8*m_prim[i]**2/m_comp[i]/(m_prim[i]+m_comp[i])/(1+2*kAM)**2

        # Tidal evolution

        # Dissipation function
        yp[0:3] += -1/tf*(dot(Ospin_prim[i],ev)/2/n*f2(e)*jv/j + 9*f1(e)*ev - 11/2*dot(Ospin_prim[i],jv/j)/n*f2(e)*ev)
        yp[3:6] += -j/tf*(dot(Ospin_prim[i],ev)/2/n*f5(e)*ev - Ospin_prim[i]/2/n*f3(e) + f4(e)*jv/j - dot(Ospin_prim[i],jv/j)/2/n*f2(e)*jv/j)

        # Non-dissipative tides
        C = kAM/n*(m_prim[i]+m_comp[i])/m_prim[i]*(r_prim[i]/a)**5
        X = -C/(1-e**2)**2*dot(Ospin_prim[i],jv/j)*dot(Ospin_prim[i],ev/e)
        Y = -C/(1-e**2)**2*dot(Ospin_prim[i],jv/j)*dot(Ospin_prim[i],qv)
        Z = 1/2*C/(1-e**2)**2*(2*dot(Ospin_prim[i],jv/j)**2-dot(Ospin_prim[i],qv)**2-dot(Ospin_prim[i],ev/e)**2)

        yp[0:3] += e*(Z*qv-Y*jv/j)
        yp[3:6] += j*(-X*qv+Y*ev/e)

        if i==0:
            I = 0.1*(m1-massc1)*r1**2 + 0.21*massc1*radc1**2
            yp[14:17] += -1/I*h1*(-X*qv+Y*ev/e)
            yp[14:17] += -1/I*h1*(-1/tf)*(dot(Ospin_prim[i],ev)/2/n*f5(e)*ev - Ospin_prim[i]/2/n*f3(e) + f4(e)*jv/j - dot(Ospin_prim[i],jv/j)/2/n*f2(e)*jv/j)
        else:
            I = 0.1*(m2-massc2)*r2**2 + 0.21*massc2*radc2**2
            yp[17:20] += -1/I*h1*(-X*qv+Y*ev/e)
            yp[17:20] += -1/I*h1*(-1/tf)*(dot(Ospin_prim[i],ev)/2/n*f5(e)*ev - Ospin_prim[i]/2/n*f3(e) + f4(e)*jv/j - dot(Ospin_prim[i],jv/j)/2/n*f2(e)*jv/j)

    ### BH spins ###

    S1v = y[20:23]
    S2v = y[23:26]

    ### SO ###

    Seff = (1+3*m2/4/m1)*S1v+(1+3*m1/4/m2)*S2v
    
    yp[0:3] += 2*G/c**2/a**3/j**3*cross(Seff-3*dot(Seff,jv/j)*jv/j,ev)
    yp[3:6] += 2*G/c**2/a**3/j**3*cross(Seff,jv)
    yp[14:17] += 2*G**(3/2)*mu1*m12**(1/2)/c**2/a**(5/2)/j**3*(1+3*m2/4/m1)*cross(jv,S1v)
    yp[17:20] += 2*G**(3/2)*mu1*m12**(1/2)/c**2/a**(5/2)/j**3*(1+3*m1/4/m2)*cross(jv,S2v)
    
    ### SS ###

    S0 = (1+m2/m1)*S1v+(1+m1/m2)*S2v

    yp[0:3] += 3*np.sqrt(G)/4/c**2/m12**(3/2)/a**(7/2)/j**4*cross(5*dot(S0,jv/j)**2*jv/j-2*dot(S0,jv/j)*S0-dot(S0,S0)*jv/j,ev)
    yp[3:6] += -3*np.sqrt(G)/2/c**2/m12**(3/2)/a**(7/2)/j**4*dot(S0,jv/j)*cross(S0,jv)
    yp[14:17] += G*mu1/2/c**2/m12/a**3/j**3*(1+m2/m1)*cross(S0-3*dot(S0,jv/j)*jv/j,S1v)
    yp[17:20] += G*mu1/2/c**2/m12/a**3/j**3*(1+m1/m2)*cross(S0-3*dot(S0,jv/j)*jv/j,S2v)

    return yp