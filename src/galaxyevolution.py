import numpy as np
import astropy.units as u
from astropy.constants import G
import gala.potential as gp

# Milky Way potential
pot = gp.MilkyWayPotential2022()

# Default units are Rsun, Msun, Myr
G = G.to(u.Rsun**3/u.Msun/u.Myr**2).value

def dot(a,b): 
    return np.sum(a*b)

def cross(a,b):
    return np.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])

def yp_spins(t,y,star1,star2,star3):
    yp = np.zeros_like(y)

    # Get masses
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)
    m3 = star3.interpolators['m'](t)
    m12 = m1+m2
    m123 = m1+m2+m3

    ### Orbital frames ###
    ev = y[0:3]
    jv = y[3:6]
    Ev = y[7:10]
    Jv = y[10:13]

    ### Semi-major axes ###
    a = y[6]
    A = y[13]

    ### Cartesian unit vectors ###
    n = np.array([[1,0,0],[0,1,0],[0,0,1]])

    for i in range(3):
        for j in range(3):
            yp[0:3] += Phi[i,j]*(dot(n[j],jv)*cross(ev,n[i])-5*dot(n[j],ev)*cross(jv,n[i])+np.kronecker(i,j)*cross(jv,ev))
            yp[3:6] += Phi[i,j]*(dot(n[j],jv)*cross(jv,n[i])-5*dot(n[j],ev)*cross(ev,n[i]))
            yp[7:10] += Phi[i,j]*(dot(n[j],Jv)*cross(Ev,n[i])-5*dot(n[j],Ev)*cross(Jv,n[i])+np.kronecker(i,j)*cross(Jv,Ev))
            yp[10:13] += Phi[i,j]*(dot(n[j],Jv)*cross(Jv,n[i])-5*dot(n[j],Ev)*cross(Ev,n[i]))

    yp[0:6] *= a**(3/2)/2/np.sqrt(m12)
    yp[7:13] *= A**(3/2)/2/np.sqrt(m123)

    return yp