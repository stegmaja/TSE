import numpy as np
import my_lib as f
from scipy import interpolate
from stepinterp import StepInterp
from common import *
from copy import copy

class BinaryStar:
    def __init__(self,z,kstar,mass0,mass,rad=(0,0),lum=(0,0),massc=(0,0),radc=(0,0),menv=(0,0),renv=(0,0),ospin=(0,0),epoch=(0,0),tms=(0,0),tphys=0.,tphysf=100.,dtp=0.,tb=1e3,ecc=0.):
        reset()

        zpars = np.zeros(20)
        f.zcnsts(z,zpars)
        f.instar()

        f.evolve(kstar,mass0,mass,rad,lum,massc,radc,menv,renv,ospin,epoch,tms,tphys,tphysf,dtp,z,zpars,tb,ecc)

        self.bcm,self.bpp = copy(f.binary.bcm),copy(f.binary.bpp)
        self.ffb,self.directcollapse,self.ecs,self.mfin = copy(f.kicksn.ffb),copy(f.kicksn.directcollapse),copy(f.kicksn.ecs),copy(f.kicksn.mfin)

        self.print()

        t = self.bcm[:,0]
        k = self.bcm[:,1]
        m = self.bcm[:,3]
        '''    
        i = 0
        self.mpre1 = m[i]
        while((k[i]<13) and (t[i]!=-1)):
            self.mpre1 = m[i]
            i = i+1
        if(k[i]>=13):
            self.mpost1 = m[i]
        else:
            self.mpost1 = -1.

        k = self.bcm[:,15]
        m = self.bcm[:,17]
            
        i = 0
        self.mpre2 = m[i]
        while((k[i]<13) and (t[i]!=-1)):
            self.mpre2 = m[i]
            i = i+1
        if(k[i]>=13):
            self.mpost2 = m[i]
        else:
            self.mpost2 = -1.
        '''
        try:
            self.cont = np.argwhere(t<=0)[0][0]
        except:
            self.cont = len(t)
        self.interpolations = {}
        self.interpolate()
        
    def print(self):
        print("###  Evolve binary star  ###")
        print("     Time Kstar1 Kstar2    Mass1   Mass10   Mass2    Mass20  LogRad1  LogRad2          Sep      Ecc   Event")
        i = 0
        if(self.bcm[0][0]>0.):
            print("%9.3f" % self.bcm[i][0],
                "%6i" % self.bcm[i][1],
                "%6i" % self.bcm[i][15],
                "%8.3f" % self.bcm[i][3],
                "%8.3f" % self.bcm[i][2],
                "%8.3f" % self.bcm[i][17],
                "%8.3f" % self.bcm[i][16],
                "%8.3f" % self.bcm[i][14],
                "%8.3f" % self.bcm[i][28],
                "%12.3f" % self.bcm[i][30],
                "%8.3f" % self.bcm[i][31],
                "  "+labels(1))
        while(self.bpp[i][0]!=0.):
            print("%9.3f" % self.bpp[i][0],
                "%6i" % self.bpp[i][1],
                "%6i" % self.bpp[i][15],
                "%8.3f" % self.bpp[i][3],
                "%8.3f" % self.bpp[i][2],
                "%8.3f" % self.bpp[i][17],
                "%8.3f" % self.bpp[i][16],
                "%8.3f" % self.bpp[i][5],
                "%8.3f" % self.bpp[i][19],
                "%12.3f" % self.bpp[i][30],
                "%8.3f" % self.bpp[i][31],
                "  "+labels(self.bpp[i][32]))
            i = i+1
        print("")    

    def interpolate(self):

        t = self.bcm[:,0]

        k1 = self.bcm[:,1]
        m1 = self.bcm[:,3]
        r1 = 10.**self.bcm[:,5]

        mass01 = self.bcm[:,2]
        lum1 = 10.**self.bcm[:,4]
        massc1 = self.bcm[:,7]
        radc1 = self.bcm[:,8]
        menv1 = self.bcm[:,9]
        renv1 = self.bcm[:,10]
        ospin1 = self.bcm[:,12]
        epoch1 = self.bcm[:,11]

        k2 = self.bcm[:,15]
        m2 = self.bcm[:,17]
        r2 = 10.**self.bcm[:,19]

        mass02 = self.bcm[:,16]
        lum2 = 10.**self.bcm[:,18]
        massc2 = self.bcm[:,21]
        radc2 = self.bcm[:,22]
        menv2 = self.bcm[:,23]
        renv2 = self.bcm[:,24]
        ospin2 = self.bcm[:,26]
        epoch2 = self.bcm[:,25]

        self.interpolations = {
            "k1":        StepInterp(t[:self.cont],k1[:self.cont]),
            "m1":        interpolate.interp1d(t[:self.cont],m1[:self.cont]),
            "rad1":      interpolate.interp1d(t[:self.cont],r1[:self.cont]),
            "mass01":    interpolate.interp1d(t[:self.cont],mass01[:self.cont]),
            "lum1":      interpolate.interp1d(t[:self.cont],lum1[:self.cont]),
            "massc1":    interpolate.interp1d(t[:self.cont],massc1[:self.cont]),
            "radc1":     interpolate.interp1d(t[:self.cont],menv1[:self.cont]),
            "menv1":     interpolate.interp1d(t[:self.cont],menv1[:self.cont]),
            "renv1":     interpolate.interp1d(t[:self.cont],renv1[:self.cont]),
            "ospin1":    interpolate.interp1d(t[:self.cont],ospin1[:self.cont]),
            "epoch1":    interpolate.interp1d(t[:self.cont],epoch1[:self.cont]),
            "k2":        StepInterp(t[:self.cont],k2[:self.cont]),
            "m2":        interpolate.interp1d(t[:self.cont],m2[:self.cont]),
            "rad2":      interpolate.interp1d(t[:self.cont],r2[:self.cont]),
            "mass02":    interpolate.interp1d(t[:self.cont],mass02[:self.cont]),
            "lum2":      interpolate.interp1d(t[:self.cont],lum2[:self.cont]),
            "massc2":    interpolate.interp1d(t[:self.cont],massc2[:self.cont]),
            "radc2":     interpolate.interp1d(t[:self.cont],menv2[:self.cont]),
            "menv2":     interpolate.interp1d(t[:self.cont],menv2[:self.cont]),
            "renv2":     interpolate.interp1d(t[:self.cont],renv2[:self.cont]),
            "ospin2":    interpolate.interp1d(t[:self.cont],ospin2[:self.cont]),
            "epoch2":    interpolate.interp1d(t[:self.cont],epoch2[:self.cont])
        }
        """
        t = t[:self.cont]

        m1 = m1[:self.cont]
        t1,ind = np.unique(t,return_index=True)
        m1 = m1[ind]

        self.interpolations["dm1"] = interpolate.interp1d(t1[:-1],-np.diff(m1)/(np.diff(t1)),bounds_error=False,fill_value=0.)

        m2 = m2[:self.cont]
        t2,ind = np.unique(t,return_index=True)
        m2 = m2[ind]

        self.interpolations["dm2"] = interpolate.interp1d(t2[:-1],-np.diff(m2)/(np.diff(t2)),bounds_error=False,fill_value=0.)
        """
    
    def get_m1(self,t):
        return self.interpolations["m1"](t)
    def get_m2(self,t):
        return self.interpolations["m2"](t)
    
    def get_star1(self,t,bpp=False,i=0):
        if(bpp):
            return self.bpp[i,1],\
                self.bpp[i,3],\
                10.**self.bpp[i,5],\
                self.bpp[i,2],\
                10.**self.bpp[i,4],\
                self.bpp[i,7],\
                self.bpp[i,8],\
                self.bpp[i,9],\
                self.bpp[i,10],\
                self.bpp[i,12],\
                self.bpp[i,11]
        else:
            return self.interpolations["k1"].interp1d(t),\
                    self.interpolations["m1"](t),\
                    self.interpolations["rad1"](t),\
                    self.interpolations["mass01"](t),\
                    self.interpolations["lum1"](t),\
                    self.interpolations["massc1"](t),\
                    self.interpolations["radc1"](t),\
                    self.interpolations["menv1"](t),\
                    self.interpolations["renv1"](t),\
                    self.interpolations["ospin1"](t),\
                    self.interpolations["epoch1"](t)

    def get_star2(self,t,bpp=False,i=0):
        if(bpp):
            return self.bpp[i,15],\
                self.bpp[i,17],\
                10.**self.bpp[i,19],\
                self.bpp[i,16],\
                10.**self.bpp[i,18],\
                self.bpp[i,21],\
                self.bpp[i,22],\
                self.bpp[i,23],\
                self.bpp[i,24],\
                self.bpp[i,26],\
                self.bpp[i,25]
        else:
            return self.interpolations["k2"].interp1d(t),\
                    self.interpolations["m2"](t),\
                    self.interpolations["rad2"](t),\
                    self.interpolations["mass02"](t),\
                    self.interpolations["lum2"](t),\
                    self.interpolations["massc2"](t),\
                    self.interpolations["radc2"](t),\
                    self.interpolations["menv2"](t),\
                    self.interpolations["renv2"](t),\
                    self.interpolations["ospin2"](t),\
                    self.interpolations["epoch2"](t)

    def get_outcome(self,interaction=True,i=0):
        # Only for RLO interaction
        if(interaction):
            if((10.**self.bpp[i,5]<=self.bpp[i,30]*(1.-self.bpp[i,31])*Roche(self.bpp[i,3],self.bpp[i,17]))
            and (10.**self.bpp[i,19]<=self.bpp[i,30]*(1.-self.bpp[i,31])*Roche(self.bpp[i,17],self.bpp[i,3]))):
                print("BSE did not find anything.\n")
            else:
                while(np.isin(self.bpp[i,32],[4,6,10,11],assume_unique=True,invert=True)\
                    and (max(self.bpp[i,1],self.bpp[i,15])!=15)): # Terminating events
                    i = i+1
                if(self.bpp[i,32]==4):
                    print("RLO stops at","%6.3f" % self.bpp[i,0],"Myr.\n")
                elif((self.bpp[i,32]==6)\
                    or (max(self.bpp[i,1],self.bpp[i,15])==15)):
                    print("Binary interaction leads to coalescence at","%6.3f" % self.bpp[i,0],"Myr.\n")
                elif(self.bpp[i,32]==11):
                    print("Binary disrupted during interaction at","%6.3f" % self.bpp[i,0],"Myr.\n")
                elif(self.bpp[i,32]==10):
                    print("Maximum integration time reached at","%6.3f" % self.bpp[i,0],"Myr.\n")
                else:
                    print("Error. Which case?\n")
            
            t = self.bpp[i,0]
            return t,labels(self.bpp[i,32]),i
        
        else:
            t = tphysf
            DCO = lambda i : np.isin(self.bpp[i,1],[13,14]) and np.isin(self.bpp[i,15],[13,14])\
                and (self.bpp[i,30]>0.) and (self.bpp[i,31]>=0.) and (self.bpp[i,31]<1.)
            if(np.isin(11,self.bpp[:,32])):
                i = np.where(self.bpp[:,32]==11)
                t = self.bpp[i,0]
                print("Binary disrupted at","%6.3f" % t,"Myr.\n")
            elif(np.isin(6,self.bpp[:,32])):
                i = np.where(self.bpp[:,32]==6)
                t = self.bpp[i,0]
                print("Binary coalescence at","%6.3f" % t,"Myr.\n")
            elif((np.isin(15,self.bpp[:,1])) or (np.isin(15,self.bpp[:,15]))):
                i = np.where(np.maximum(self.bpp[:,1],self.bpp[:,15])==15)[0][0]
                t = self.bpp[i,0]
                print("Binary coalescence at","%6.3f" % t,"Myr.\n")
            elif(DCO(np.where(self.bpp[:,32]==10))):
                i = 0
                while(not(DCO(i))):
                    i = i+1
                t = self.bpp[i,0]
                print("DCO formed at","%6.3f" % t,"Myr.\n")
            else:
                print("Maximum integration time reached at","%6.3f" % self.bpp[np.where(self.bpp[:,32]==10),0],"Myr.\n")
                i = np.where(self.bpp[:,32]==10)
                t = self.bpp[i,0]

            return t,labels(self.bpp[i,32]),i