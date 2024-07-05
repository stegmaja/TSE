import numpy as np
import astropy.units as u
from astropy.constants import G
import gala.potential as gp
import gala.dynamics as gd
import gala.integrate as gi
from gala.units import galactic
from initialconditions import ic
from scipy import interpolate
from common import dot,cross

# Default units are Rsun, Msun, Myr
G = G.to(u.Rsun**3/u.Msun/u.Myr**2).value

class Galaxy:
    def __init__(self, pot = gp.MilkyWayPotential2022(units=galactic)):
        ics = gd.PhaseSpacePosition(pos=[ic.x_MW,ic.y_MW,ic.z_MW] * u.kpc,
                                    vel=[ic.vx_MW,ic.vy_MW,ic.vz_MW] * u.km/u.s)
        orbit = gp.Hamiltonian(pot).integrate_orbit(ics, t1=0, t2=ic.max_time, dt=1., Integrator=gi.DOPRI853Integrator)
        self.x = interpolate.interp1d(orbit.t, orbit.x)
        self.y = interpolate.interp1d(orbit.t, orbit.y)
        self.z = interpolate.interp1d(orbit.t, orbit.z)
        self.pot = pot

    def hessian(self, t):
        x = self.x(t)
        y = self.y(t)
        z = self.z(t)
        w = np.array([x,y,z])
        return self.pot.hessian(w).value
    
Gal = Galaxy()

def yp_galaxy(t,y,star1,star2,star3):
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

    ### Hessian ###
    hessian = Gal.hessian(t)

    for i in range(3):
        for j in range(3):
            yp[0:3] += hessian[i,j][0]*(dot(n[j],jv)*cross(ev,n[i])-5*dot(n[j],ev)*cross(jv,n[i])+int(i==j)*cross(jv,ev))
            yp[3:6] += hessian[i,j][0]*(dot(n[j],jv)*cross(jv,n[i])-5*dot(n[j],ev)*cross(ev,n[i]))
            yp[7:10] += hessian[i,j][0]*(dot(n[j],Jv)*cross(Ev,n[i])-5*dot(n[j],Ev)*cross(Jv,n[i])+int(i==j)*cross(Jv,Ev))
            yp[10:13] += hessian[i,j][0]*(dot(n[j],Jv)*cross(Jv,n[i])-5*dot(n[j],Ev)*cross(Ev,n[i]))

    yp[0:6] *= a**(3/2)/2/np.sqrt(G*m12)
    yp[7:13] *= A**(3/2)/2/np.sqrt(G*m123)

    return yp