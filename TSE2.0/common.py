import numpy as np
import my_lib as f
from numba import jit

#Myr = np.sqrt(3.942e13)
Myr = 6278534.860936905
yr = Myr/1e6
yeardy = 365.24
Rsol = 1./214.95
kms = 210805./Myr
c = 6.32e10/Myr
c2 = c**2
c5 = c**5
#z = 0.02
z = 0.0002
zpars = np.zeros(20)
y = np.zeros(20)
kL1 = .028
kL2 = .028
k1 = .1
k2 = .1
Dt1 = 1.996e-7 # 1 Sec lag time
Dt2 = 1.996e-7 # 1 Sec lag time
Dt = np.array([Dt1,Dt2])
kL = np.array([kL1,kL2])
k = np.array([k1,k2])
tphysf = 100.

label = ['INITIAL','KW_CHNGE','BEG_RCHE','END_RCHE','CONTACT ','COELESCE','COMENV  ','GNTAGE  ','NO_REMNT','MAX_TIME','DISRUPT ','BEG_SYMB','END_SYMB','BEG_BSS']
labels = lambda i : label[int(i-1)]

def reset():
    f.value1.neta,f.value1.bwind,f.value1.hewind,f.value2.alpha1,f.value2.lambda1 = 0.5, 0.0, 1.0, 3.0, 0.1
    f.flags.ceflag,f.flags.tflag,f.flags.ifflag,f.flags.wdflag,f.value4.bhflag,f.flags.nsflag,f.flags.piflag,f.value1.mxns,f.value3.idum = 0, 1, 0, 2, 2, 1, 1, 3.0, np.random.randint(678555)
    f.value4.bhflag = 1
    f.points.pts1,f.points.pts2,f.points.pts3 = 0.05,0.01,0.02
    f.value4.sigma1,f.value4.sigma2,f.value5.beta,f.value5.xi,f.value5.acc2,f.value5.epsnov,f.value5.eddfac,f.value5.gamma = 265.0, 265.0, 0.125, 1.0, 1.5, 0.001, 1.0, -1.0

    f.types.ktype = np.zeros_like(f.types.ktype)
    f.load.saveflag = np.zeros_like(f.load.saveflag)
    f.kicksn.ffb = np.zeros_like(f.kicksn.ffb)
    f.kicksn.directcollapse = np.zeros_like(f.kicksn.directcollapse)
    f.kicksn.ecs = np.zeros_like(f.kicksn.ecs)
    f.kicksn.mfin = np.zeros_like(f.kicksn.mfin)
    f.rand3.idum2 = np.zeros_like(f.rand3.idum2)
    f.rand3.iy = np.zeros_like(f.rand3.iy)
    f.rand3.ir = np.zeros_like(f.rand3.ir)
    f.mscff.msp = np.zeros_like(f.mscff.msp)
    f.gbcff.gbp = np.zeros_like(f.gbcff.gbp)
    f.corer.coef = np.zeros_like(f.corer.coef)
    f.binary.bcm = np.zeros_like(f.binary.bcm)
    f.binary.bpp = np.zeros_like(f.binary.bpp)
    
    if(f.value3.idum>0): f.value3.idum = -int(f.value3.idum)

@jit(nopython=True,cache=True,fastmath=True)
def Roche(m1,m2):
    p = (m1/m2)**(1./3.)
    return .49*p*p/(.6*p*p+np.log(1.+p))

@jit(nopython=True,cache=True,fastmath=True)
def stability(m12,m3,a1,a2,e2):
    stable = True
    MA=2.8*((1.+m3/m12)*(1.+e2)/np.sqrt(1.-e2))**(2./5.) #Mardling+Aarseth
    if(a2*(1.-e2)<MA*a1): stable=False
    return stable

#@jit(nopython=True,cache=True,fastmath=True)
def RLO(m12,m3,a2,e2,r3,filling_factor=1.01):
    fill = False
    if(r3>filling_factor*a2*(1.-e2)*Roche(m3,m12)): fill=True
    return fill

@jit(nopython=True,cache=True,fastmath=True)
def cross(a,b):
    return np.array([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])

@jit(nopython=True,cache=True,fastmath=True)
def power_law(alpha,low,high,harv=np.random.uniform()):
    N=(alpha+1.)/(high**(alpha+1.)-low**(alpha+1.)) # Normalise
    return (harv*(alpha+1.)/N+low**(alpha+1.))**(1./(alpha+1.))

#
########################################################################
#
#     Eq. 5
#
@jit(nopython=True,cache=True,fastmath=True)
def Ftwin_func(m1,logPrd):
    Ftwinle1=0.3-0.15*np.log10(m1)
    if(m1>6.5):
        logPtwin=1.5
    else:
        logPtwin=8.-m1
    
    if(logPrd<1.):
        Ftwin_func=Ftwinle1
    elif(logPrd<logPtwin):
        Ftwin_func=Ftwinle1*(1.-(logPrd-1.)/(logPtwin-1.))
    else:
        Ftwin_func=0.
    
    return Ftwin_func
#
########################################################################
#
#     Eq. 11, valid if M>6
#
@jit(nopython=True,cache=True,fastmath=True)
def gamma_lq_func(m1,logPrd):

    if(m1<=6.):
        print('Error: primary too light. Eqs. are only valid if M>6')
    if(logPrd<1.):
        gamma_lq_func=-0.5
    elif(logPrd<2.):
        gamma_lq_func=-0.5-0.9*(logPrd-1.)
    elif(logPrd<4.):
        gamma_lq_func=-1.4-0.3*(logPrd-2.)
    else:
        gamma_lq_func=-2.

    return gamma_lq_func
#
########################################################################
#
#     Eq. 15, valid if M>6
#
@jit(nopython=True,cache=True,fastmath=True)
def gamma_sq_func(m1,logPrd):

    if(logPrd<1.):
        gamma_sq_func=0.1
    elif(logPrd<3.):
        gamma_sq_func=0.1-0.15*(logPrd-1.)
    elif(logPrd<5.6):
        gamma_sq_func=-0.2-0.5*(logPrd-3.)
    else:
        gamma_sq_func=-1.5
      
    return gamma_sq_func
#
########################################################################
#
#     Eq. 3
#
@jit(nopython=True,cache=True,fastmath=True)
def emax_func(Prd):
    return 1.-(Prd/2.)**(-2./3.)
#
########################################################################
#
#     Eq. 18, valid if M>7
#
@jit(nopython=True,cache=True,fastmath=True)
def eta_func(m1,logPrd):
    # Note that eta is not well constrained beyond logP>5, see Sec. 9.2.
    return 0.9-0.2/(logPrd-0.5)