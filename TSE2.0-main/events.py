'''
Define terminating events during the integration
'''

from numba import jit
import numpy as np
from numpy.linalg import norm
from units import *
import globals

# Returns aout and eout
#@jit(nopython=True,fastmath=True,cache=True)
def Kelperout(y,m1,m2,m3):
    m123 = m1+m2+m3
    routv = y[7:10]
    voutv = y[10:13]
    eoutv = np.cross(voutv,np.cross(routv,voutv))/G/m123-routv/norm(routv) 
    eout = norm(eoutv)
    Loutv = np.cross(routv,voutv)
    aout = norm(Loutv)**2/G/m123/(1.-eout**2)
    return aout,eout

# Roche lobe factor
#@jit(nopython=True,cache=True,fastmath=True)
def Roche(m1,m2):
    p = (m1/m2)**(1./3.)
    return .49*p*p/(.6*p*p+np.log(1.+p))

# DCO merger
#@jit(nopython=True,fastmath=True,cache=True)
def merger(t,y,limit=1e-3*AU):
    return y[6]-limit
merger.terminal = True
#merger.direction = -1.

# Dynamical instablility
def stability(t,y):
    return -1
    m1 = globals.star1.get_m(t/Myr)*Msol
    m2 = globals.star2.get_m(t/Myr)*Msol
    m3 = globals.star3.get_m(t/Myr)*Msol
    m12 = m1+m2
    m123 = m12+m3
    a = y[6]
    aout,eout = Kelperout(y,m1,m2,m3)
    if(eout>1.):print("Error in stability crit.")
    MA=2.8*((1.+m3/m12)*(1.+eout)/np.sqrt(1.-eout))**(2./5.) #Mardling+Aarseth
    return aout*(1.-eout)-MA*a
stability.terminal = False
stability.direction = -1.

# Primary RLO
def PRLO(t,y,filling_factor=1.):
    m1 = globals.star1.get_m(t/Myr)*Msol
    m2 = globals.star2.get_m(t/Myr)*Msol
    r1 = globals.star1.get_r(t/Myr)*Rsol
    m12 = m1+m2
    ev = y[0:3]
    e = norm(ev)
    a = y[6]
    if(r1>=filling_factor*a*(1.-e)*Roche(m1,m2)):print(r1/AU,filling_factor*a*(1.-e)*Roche(m1,m2)/AU)
    return r1-filling_factor*a*(1.-e)*Roche(m1,m2)
PRLO.terminal = True
#PRLO.direction = -1.

# Secondary RLO
def SRLO(t,y,filling_factor=1.):
    m1 = globals.star1.get_m(t/Myr)*Msol
    m2 = globals.star2.get_m(t/Myr)*Msol
    r2 = globals.star2.get_r(t/Myr)*Rsol
    m12 = m1+m2
    ev = y[0:3]
    e = norm(ev)
    a = y[6]
    return r2-filling_factor*a*(1.-e)*Roche(m2,m1)
SRLO.terminal = True
#SRLO.direction = -1.

# Tertiary RLO
def TRLO(t,y,filling_factor=1.):
    m1 = globals.star1.get_m(t/Myr)*Msol
    m2 = globals.star2.get_m(t/Myr)*Msol
    m3 = globals.star3.get_m(t/Myr)*Msol
    r3 = globals.star3.get_r(t/Myr)*Rsol
    m12 = m1+m2
    r3 = globals.r3
    aout,eout = Kelperout(y,m1,m2,m3)
    return r3-filling_factor*aout*(1.-eout)*Roche(m3,m12)
TRLO.terminal = True
#TRLO.direction = -1.

# Primary SN
def PSN(t,y):
    return globals.star1.tco*Myr-t
PSN.terminal = True
PSN.direction = -1.

# Secondary SN
def SSN(t,y):
    return globals.star2.tco*Myr-t
SSN.terminal = True
SSN.direction = -1.

# Tertiary SN
def TSN(t,y):
    return globals.star3.tco*Myr-t
TSN.terminal = True
TSN.direction = -1.

# Unbound outer orbit
def TUnbound(t,y,factor=10.):
    m1 = globals.star1.get_m(t/Myr)*Msol
    m2 = globals.star2.get_m(t/Myr)*Msol
    m3 = globals.star3.get_m(t/Myr)*Msol
    m123 = m1+m2+m3
    routv = y[7:10]
    rout = norm(routv)
    voutv = y[10:13]
    vout = norm(voutv)
    return factor*G*m123/rout-vout**2/2
TUnbound.terminal = False
TUnbound.direction = -1.