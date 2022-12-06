import numpy as np
import my_lib as f
from scipy import interpolate
from stepinterp import StepInterp
from common import *

class IsolatedStar:   
    def __init__(self,z,kstar,mass0,mass,rad=0.,lum=0.,massc=0.,radc=0.,menv=0.,renv=0.,ospin=0.,epoch=0.,tphys=0.,tphysf=100.):
        reset()

        if(kstar>=12):
            self.js = True
        else:
            self.js = False

        ### Place a BH companion at extremely large distance s.t. the star evolves effectively single. ###
        kstar = (14,int(kstar)) # BH companion
        mass0 = (1.,mass0)
        mass = (1.,mass)
        rad = (0.,rad)
        lum = (0.,lum)
        massc = (0.,massc)
        radc = (0.,radc)
        menv = (0.,menv)
        renv = (0.,renv)
        ospin = (0.,ospin)
        epoch = (0.,epoch)
        tms = (0.,0.)
        dtp = 0.
        tb = 1e9 # extremely large distance
        ecc = 0.

        zpars = np.zeros(20)
        f.zcnsts(z,zpars)
        f.instar()

        f.evolve(kstar,mass0,mass,rad,lum,massc,radc,menv,renv,ospin,epoch,tms,tphys,tphysf,dtp,z,zpars,tb,ecc)

        self.bcm,self.bpp = f.binary.bcm,f.binary.bpp
        self.ffb,self.directcollapse,self.ecs,self.mfin = f.kicksn.ffb,f.kicksn.directcollapse,f.kicksn.ecs,f.kicksn.mfin 

        self.print()

        t = self.bcm[:,0]
        k = self.bcm[:,1+14]
        m = self.bcm[:,3+14]
            
        i = 0
        self.mpre = m[i]
        while((k[i]<13) and (t[i]!=-1)):
            self.mpre = m[i]
            i = i+1
        if(k[i]>=13):
            self.mpost = m[i]
            self.tco = t[i]
            self.kco = k[i]
        else:
            self.mpost = -1.
            self.tco = -1.
            self.kco = -1

        self.cont = np.argwhere(t==-1)[0][0]
        self.interpolations = {}
        self.interpolate()
        
    def print(self):

        print("###  Evolve single star  ###")
        print("     Time Kstar     Mass    Mass0   LogRad")
        i = 0
        print("%9.3f" % self.bcm[i][0],
                "%5i" % self.bcm[i][1+14],
                "%8.3f" % self.bcm[i][3+14],
                "%8.3f" % self.bcm[i][2+14],
                "%8.3f" % self.bcm[i][5+14])
        while(self.bpp[i][0]!=-1):
            print("%9.3f" % self.bpp[i][0],
                "%5i" % self.bpp[i][1+14],
                "%8.3f" % self.bpp[i][3+14],
                "%8.3f" % self.bpp[i][2+14],
                "%8.3f" % self.bpp[i][5+14])
            i = i+1
        print("")

    def interpolate(self):

        t = self.bcm[:,0]

        k = self.bcm[:,1+14]
        m = self.bcm[:,3+14]
        r = 10.**self.bcm[:,5+14]

        mass0 = self.bcm[:,2+14]
        lum = 10.**self.bcm[:,4+14]
        massc = self.bcm[:,7+14]
        radc = self.bcm[:,8+14]
        menv = self.bcm[:,9+14]
        renv = self.bcm[:,10+14]
        ospin = self.bcm[:,12+14]
        epoch = self.bcm[:,11+14]

        self.interpolations = {
            #"k":        StepInterp(t[:self.cont],k[:self.cont]),
            "k":        interpolate.interp1d(t[:self.cont],k[:self.cont]),
            "m":        interpolate.interp1d(t[:self.cont],m[:self.cont]),
            "rad":      interpolate.interp1d(t[:self.cont],r[:self.cont]),
            "mass0":    interpolate.interp1d(t[:self.cont],mass0[:self.cont]),
            "lum":      interpolate.interp1d(t[:self.cont],lum[:self.cont]),
            "massc":    interpolate.interp1d(t[:self.cont],massc[:self.cont]),
            "radc":     interpolate.interp1d(t[:self.cont],menv[:self.cont]),
            "menv":     interpolate.interp1d(t[:self.cont],menv[:self.cont]),
            "renv":     interpolate.interp1d(t[:self.cont],renv[:self.cont]),
            "ospin":    interpolate.interp1d(t[:self.cont],ospin[:self.cont]),
            "epoch":    interpolate.interp1d(t[:self.cont],epoch[:self.cont])
        }

        t = t[:self.cont]
        m = m[:self.cont]
        t,ind = np.unique(t,return_index=True)
        m = m[ind]

        self.interpolations["dm"] = interpolate.interp1d(t[:-1],-np.diff(m)/(np.diff(t)),bounds_error=False,fill_value=0.)

    def get_basics(self,t):
        return int(self.interpolations["k"](t)),\
                self.interpolations["m"](t),\
                self.interpolations["rad"](t),\
                self.interpolations["dm"](t)

    def get_k(self,t):
        return int(self.interpolations["k"](t))

    def get_m(self,t):
        return self.interpolations["m"](t)

    def get_dm(self,t):
        return self.interpolations["dm"](t)
    
    def get_r(self,t):
        return self.interpolations["rad"](t)

    def get_full(self,t):
        return self.interpolations["mass0"](t),\
                self.interpolations["lum"](t),\
                self.interpolations["massc"](t),\
                self.interpolations["radc"](t),\
                self.interpolations["menv"](t),\
                self.interpolations["renv"](t),\
                self.interpolations["ospin"](t),\
                self.interpolations["epoch"](t)


class IsolatedCO:   
    def __init__(self,kstar,m):
        self.m = m
        self.k = kstar
        self.tco = -1.
        self.r = 1e-9*Rsol

    def get_k(self,t):
        return self.k

    def get_m(self,t):
        return self.m

    def get_dm(self,t):
        return 0.
    
    def get_r(self,t):
        return self.r

    def get_basics(self,t):
        return self.k,self.m,self.r,0.

    def get_full(self,t):
        return self.m,0.,0.,0.,0.,0.,0.,0.

    '''While stellar spins are described by y[]...'''
    def spin_conversion(Omv,jv,alignment="aligned",magnitude="unity"):
        if(alignment=="aligned"):
            uv = jv/np.linalg.norm(jv)
        elif(alignment=="random"):
            uv = np.random.uniform(size=3)
            uv /= np.linalg.norm(uv)
        else:
            print("Error: Which case are you looking for?")
        
        if(magnitude=="unity"):
            return uv
        else:
            print("Error: Which case are you looking for?")
        
