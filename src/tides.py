import numpy as np

def f1(e):
    return (1-e**2)**(-13/2)*(1+15/4*e**2+15/8*e**4+5/64*e**6)
def f2(e):
    return (1-e**2)**(-5)*(1+3/2*e**2+1/8*e**4)
def f3(e):
    return (1-e**2)**(-5)*(1+9/2*e**2+5/8*e**4)
def f4(e):
    return (1-e**2)**(-13/2)*(1+15/2*e**2+45/8*e**4+5/16*e**6)
def f5(e):
    return (1-e**2)**(-5)*(3+5*e**2)
def f6(e):
    return (1-e**2)**(-8)*(1+31/2*e**2+255/8*e**4+185/16*e**6+25/64*e**8)

def tidal_evolve(t,y,star1,star2,star3):
    yp = np.zeros_like(y)

    # Get masses
    m1 = star1.interpolators['m'](t)
    m2 = star2.interpolators['m'](t)

    # Get radii
    r1 = 10**star1.interpolators['logr'](t)
    r2 = 10**star2.interpolators['logr'](t)

    # Get stellar types
    k1 = star1.interpolators['k'](t)
    k2 = star2.interpolators['k'](t)

    # Arrays
    m_prim = np.array([m1,m2])
    r_prim = np.array([r1,r2])
    k_prim = np.array([k1,k2])
    m_comp = np.flip(m_prim)

    # Semi-major axis
    a = y[6]
    ev = y[0:3]
    e = np.linalg.norm(ev)

    for i in range(2):



    return yp