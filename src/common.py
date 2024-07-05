import numpy as np

####################################################################################################
# Useful functions
####################################################################################################
def dot(a,b): 
    return np.sum(a*b)

def cross(a,b):
    return np.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])

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