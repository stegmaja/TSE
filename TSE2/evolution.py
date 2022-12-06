from numpy import zeros_like,sqrt,sum
from numpy.linalg import norm
import numpy as np
from numba import jit
import globals
from units import *

@jit(nopython=True,fastmath=True,cache=True)
def dot(a,b): 
    return np.sum(a*b)

@jit(nopython=True,fastmath=True,cache=True)
def cross(a,b):
    return np.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])

#@jit(nopython=True,cache=True)
def evolve(t,y,CO1=False,CO2=False,CO3=False):

    m1,m2,m3 = globals.m1,globals.m2,globals.m3
    r1,r2,r3 = globals.r1,globals.r2,globals.r3

    m12 = m1+m2
    mu = m1*m2/m12
    m123 = m12+m3
    muout = m12*m3/m123
    
    yp = zeros_like(y)

    ev = y[0:3]
    e = norm(ev)
    e = np.clip(e,1e-7,1-1e-7)
    
    jv = y[3:6]
    j = norm(jv)
    
    a = y[6]
    
    #f = sqrt(G*m123/aout**3)*t
    #y[7:10] = r*(u1*np.cos(f+om)+u2*np.sin(f+om))

    routv = y[7:10]
    rout = norm(routv)
    
    voutv = y[10:13]

    S1v = y[13:16]
    S1 = norm(S1v)
    
    S2v = y[16:19]
    S2 = norm(S2v)

    S0v = (1.+m2/m1)*S1v+(1.+m1/m2)*S2v
    Seffv = (1.+3.*m2/4./m1)*S1v+(1.+3.*m1/4./m2)*S2v

    Om1v = y[19:22]
    Om2v = y[22:25]

    n = sqrt(G*m12/a**3)
      
    tLK = m12*rout**3/n/m3/a**3
    
    eps = (m1-m2)*a/m12/rout
    
    ### Quad ###
    
    yp[0:3] += 3.*(5.*dot(ev,routv/rout)*cross(jv,routv/rout)
                   -dot(jv,routv/rout)*cross(ev,routv/rout)
                   -2.*cross(jv,ev))/2./tLK
    yp[3:6] += 3.*(5.*dot(ev,routv/rout)*cross(ev,routv/rout)
                   -dot(jv,routv/rout)*cross(jv,routv/rout))/2./tLK
    
    ### Oct ###
    
    yp[0:3] += 15.*eps*(16.*dot(ev,routv/rout)*cross(jv,ev/e)
                       -(1.-8.*e**2)*cross(jv,routv/rout)
                       +5.*dot(jv,routv/rout)**2*cross(jv,routv/rout)
                       -35.*dot(ev,routv/rout)**2*cross(jv,routv/rout)
                       +10.*dot(jv,routv/rout)*dot(ev,routv/rout)*cross(ev,routv/rout))/16./tLK
    yp[3:6] += 15.*eps*(10.*dot(jv,routv/rout)*dot(ev,routv/rout)*cross(jv,routv/rout)
                        -(1.-8.*e**2)*cross(ev,routv/rout)
                        +5.*dot(jv,routv/rout)**2*cross(ev,routv/rout)
                        -35.*dot(ev,routv/rout)**2*cross(ev,routv/rout))/16./tLK
    
    ### rout ###
    
    yp[7:10] += voutv
    yp[10:13] += -G*m123*routv/rout**3 # Keplerian

    Phi0 = G*m3*a**2/rout**3

    yp[10:13] += -mu*Phi0/4./muout*(6.*dot(routv,jv)*jv/rout**2-6.*dot(routv,jv)**2*routv/rout**4
                                  -30.*dot(routv,ev)*ev/rout**2+30.*dot(routv,ev)**2*routv/rout**4) # Quad
    yp[10:13] += -5.*mu*Phi0*eps/16./muout*((3.-24.*e**2)*(ev/rout-routv*dot(ev,routv)/rout**3)
                                  +45.*routv*dot(ev,routv)*dot(jv,routv)**2/rout**5-30.*jv*dot(ev,routv)*dot(jv,routv)/rout**3-15.*ev*dot(jv,routv)**2/rout**3
                                  +115.*ev*dot(ev,routv)**2/rout**3-115.*routv*dot(ev,routv)**3/rout**5) # Oct
    
    ### Mass loss ###

    yp[6] += a*np.abs(globals.dm1)/m12
    yp[6] += a*np.abs(globals.dm2)/m12
    
    ### GW ###
    
    yp[0:3] += -304.*ev*G**3*m1*m2*m12*(1.+121.*e**2/304.)/15./c**5/a**4/(1.-e**2)**2.5
    yp[3:6] += -32.*G**3.5*mu**2*m12**2.5*(1.+7.*e**2/8.)*jv/5./c**5/a**3.5/(1.-e**2)**2/j
    yp[6] += -64.*G**3*mu*m12**2*(1.+73.*e**2/24.+37.*e**4/96.)/5./c**5/a**3/(1.-e**2)**3.5

    ### GR ###
    
    yp[0:3] += 3.*G*n*m12*cross(jv,ev)/c**2/a/(1.-e**2)**1.5

    return yp

    ### SO ###

    if(CO1 and CO2): 
        yp[0:3] += 2.*G*cross((Seffv-3.*dot(Seffv,jv/j)*jv/j),ev)/c**2/a**3/(1.-e**2)**1.5
        yp[3:6] += 2.*G*cross(Seffv,jv)/c**2/a**3/(1.-e**2)**1.5
        yp[13:16] += 2.*G*1.5*mu*sqrt(m12)*(1.+3.*m2/4./m1)*cross(jv,S1v)/c**2/a**2.5/(1.-e**2)**1.5
        yp[16:19] += 2.*G*1.5*mu*sqrt(m12)*(1.+3.*m1/4./m2)*cross(jv,S2v)/c**2/a**2.5/(1.-e**2)**1.5

    ### SS ###

        yp[0:3] += 3.*sqrt(G)*cross(5.*dot(S0v,jv/j)**2*jv/j-2.*dot(S0v,jv/j)*S0v-dot(S0v,S0v)*jv/j,ev)/4./c**2/m12**1.5/a**3.5/(1.-e**2)**2
        yp[3:6] += -3.*sqrt(G)*dot(S0v,jv/j)*cross(S0v,jv)/2./c**2/m12**1.5/a**3.5/(1.-e**2)**2
        yp[13:16] += G*mu*(1.+m2/m1)*cross(S0v-3.*dot(S0v,jv/j)*jv/j,S1v)/2./c**2/m12/a**3/(1.-e**2)**1.5
        yp[16:19] += G*mu*(1.+m1/m2)*cross(S0v-3.*dot(S0v,jv/j)*jv/j,S2v)/2./c**2/m12/a**3/(1.-e**2)**1.5

    ### Tides and rotation ###

    uv = ev/e
    hv = jv/j
    qv = cross(hv,uv)

    Ome = np.array([dot(Om1v,uv),dot(Om2v,uv)])
    Omh = np.array([dot(Om1v,hv),dot(Om2v,hv)])
    Omq = np.array([dot(Om1v,qv),dot(Om2v,qv)])

    k = np.array([k1,k2])
    kL = np.array([kL1,kL2])
    m = np.array([m1,m2])
    r = np.array([r1,r2])
    I = k*m*r**2
    omega = np.sqrt(G*m12/a**3)
    h = np.sqrt(a*G*m12*(1.-e**2))
    f1 = 1.+3./2.*e**2+1./8.*e**4
    A = r**5*kL/(1.-kL)
    tV = 3.*(1.+1./kL)*r**3/G/m12/Dt
    tF_eq = tV*(a/r)**8*m**2*(1.-kL)**2/9./m12/np.flip(m)
    tF_dyn = (a/r)**9*np.sqrt(a**3/G/m)*m/np.flip(m)*(1.+np.flip(m)/m)**(11/6)/(9.*1.592e-9)*(m/Msol)**2.84

    one_over_tF = np.zeros(2)

    if(tF_eq[0]>tF_dyn[0]):print("Shout")
    if(tF_eq[1]>tF_dyn[1]):print("Shout")

    # Check if equilibrium or dynamical tides apply
    if(((globals.k1==1) & (m1>=1.25)) | (globals.k1==4) | (globals.k1==7)): 
        one_over_tF[0] = 1./tF_dyn[0]
    elif((globals.k1==0) | (globals.k1==2) | (globals.k1==3) | (globals.k1==5) | (globals.k1==6) | (globals.k1==8) | (globals.k1==9)):
        one_over_tF[0] = 1./tF_eq[0]
    if(((globals.k2==1) & (m2>=1.25)) | (globals.k2==4) | (globals.k2==7)): 
        one_over_tF[1] = 1./tF_dyn[1]
    elif((globals.k2==0) | (globals.k2==2) | (globals.k2==3) | (globals.k2==5) | (globals.k2==6) | (globals.k2==8) | (globals.k2==9)):
        one_over_tF[1] = 1./tF_eq[1]


    X = -np.flip(m)*A*Omh*Ome/2./mu/omega/a**5/(1.-e**2)**2\
    -Omq*(1.+9./2.*e**2+5./8.*e**4)/2./omega*one_over_tF/(1.-e**2)**5

    Y = -np.flip(m)*A*Omh*Omq/2./mu/omega/a**5/(1.-e**2)**2\
        +Ome*f1/2./omega*one_over_tF/(1.-e**2)**5

    Z = np.flip(m)*A/2./mu/omega/a**5*((2.*Omh**2-Ome**2-Omq**2)/2./(1.-e**2)**2\
        +15.*np.flip(m)*f1/a**3/(1.-e**2)**5)
        
    V = 9.*one_over_tF*((1.+15./4.*e**2+15./8.*e**4+5./64.*e**6)/(1.-e**2)**(13./2.)\
        -11.*Omh*f1/18./omega/(1.-e**2)**5)

    W = 1.*one_over_tF*((1.+15./2.*e**2+45./8.*e**4+5./16.*e**6)/(1.-e**2)**(13./2.)\
        -Omh*(1.+3.*e**2+3./8.*e**4)/18./omega/(1.-e**2)**5)

    K = np.array([sum(X),sum(Y),sum(Z)])

    #print(norm(K),1./tLK)

    yp[0:3] += -sum(V)*ev+cross(K,ev)
    yp[3:6] += -(e/j)**2*sum(V)*jv+cross(K,jv)
    yp[6] += -2.*a*(sum(W)+sum(V)*(e/j)**2)
    yp[19:22] += mu*h/I[0]*(-Y[0]*uv+X[0]*qv+W[0]*hv)
    yp[22:25] += mu*h/I[1]*(-Y[1]*uv+X[1]*qv+W[1]*hv)
    return yp
