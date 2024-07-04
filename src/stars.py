import numpy as np
from scipy import interpolate
import os
import sys
from initialconditions import ic

class SingleStar:   
    def __init__(self,
                 tphys=0.,
                 mass0_1=1.,
                 mass0_2=0.001,
                 mass_1=1.,
                 mass_2=0.001,
                 epoch_1=0.,
                 epoch_2=0.,
                 ospin_1=0.,
                 ospin_2=0.,
                 max_time=15000.,
                 period=1e10,
                 type1=1,
                 type2=14,
                 Z=ic.Z,
                 ecc=0.,
                 neta=ic.neta,
                 bwind=ic.bwind,
                 hewind=ic.hewind,
                 alpha1=ic.alpha1,
                 lambda_=ic.lamb,
                 ceflag=ic.ceflag,
                 tflag=ic.tflag,
                 ifflag=ic.ifflag,
                 wdflag=ic.wdflag,
                 bhflag=ic.bhflag,
                 nsflag=ic.nsflag,
                 piflag=ic.piflag,
                 mxns=ic.mxns,
                 idum=ic.idum,
                 pts1=ic.pts1,
                 pts2=ic.pts2,
                 pts3=ic.pts3,
                 sigma1=ic.sigma1,
                 sigma2=ic.sigma2,
                 beta=ic.beta,
                 xi=ic.xi,
                 acc2=ic.acc2,
                 epsnov=ic.epsnov,
                 eddfac=ic.eddfac,
                 gamma=ic.gamma):

        # Write input file
        with open(ic.MOBSE_DIR+'/../input/binary.input','w') as f:
            f.write(f'{mass0_1} {mass0_2} {max_time} {period} {type1} {type2} {Z} {ecc}\n')
            f.write(f'{neta} {bwind} {hewind} {alpha1} {lambda_}\n')
            f.write(f'{ceflag} {tflag} {ifflag} {wdflag} {bhflag} {nsflag} {piflag} {mxns} {idum}\n')
            f.write(f'{pts1} {pts2} {pts3}\n')
            f.write(f'{sigma1} {sigma2} {beta} {xi} {acc2} {epsnov} {eddfac} {gamma}\n')

            f.write(f'{tphys}\n')
            f.write(f'{epoch_1} {mass_1} {ospin_1}\n')
            f.write(f'{epoch_2} {mass_2} {ospin_2}\n')

            f.write('\n')

            # Write comments
            f.write('# m_zams1  m_zams2  max_time  period  type1  type2  Z  ecc\n')
            f.write('# neta  bwind  hewind  alpha1  lambda\n')
            f.write('# ceflag  tflag  ifflag  wdflag  bhflag  nsflag piflag mxns  idum\n')
            f.write('# pts1  pts2  pts3\n')
            f.write('# sigma1  sigma2  beta  xi  acc2  epsnov  eddfac  gamma\n')

        # Run MOBSE
        os.chdir(ic.MOBSE_DIR)
        os.system('./mobse.x')
        os.chdir(ic.SRC_DIR)

        # Read output
        t,k,m0,m,logr,epoch,ospin,ffb = np.loadtxt(ic.MOBSE_DIR+'/../mobse.out',unpack=True,usecols=(0,1,2,3,5,11,12,-1))

        # Only keep values where t is unique
        t, idx = np.unique(t,return_index=True)
        k = k[idx]
        m0 = m0[idx]
        m = m[idx]
        logr = logr[idx]
        epoch = epoch[idx]
        ospin = ospin[idx]
        ffb = ffb[idx]

        # Record fallback fraction
        self.ffb = ffb[0]

        # Write interpolators in time
        self.interpolators = {}
        self.interpolators['k'] = interpolate.interp1d(t,k,fill_value=(k[0],k[-1]),bounds_error=False,kind='previous')
        self.interpolators['m0'] = interpolate.interp1d(t,m0,fill_value=(m0[0],m0[-1]),bounds_error=False)
        self.interpolators['m'] = interpolate.interp1d(t,m,fill_value=(m[0],m[-1]),bounds_error=False)
        self.interpolators['logr'] = interpolate.interp1d(t,logr,fill_value=(logr[0],logr[-1]),bounds_error=False)
        self.interpolators['epoch'] = interpolate.interp1d(t,epoch,fill_value=(epoch[0],epoch[-1]),bounds_error=False)
        self.interpolators['ospin'] = interpolate.interp1d(t,ospin,fill_value=(ospin[0],ospin[-1]),bounds_error=False)
        self.interpolators['dm'] = interpolate.interp1d(t[:-1],np.diff(m)/(np.diff(t)),bounds_error=False,fill_value=0.)
        
        # Test if star is CO initially
        if k[0]>=13:
            self.initial_CO = True
        else:
            self.initial_CO = False
        
        # Get time of CO formation
        if np.any(k>=13):
            self.tco = t[np.where(k>=13)[0][0]]
            self.m_SN = m[(k<=12) & (m>0)][-1]
            self.dm_SN = m[(k<=12) & (m>0)][-1] - m[k>=13][0]
        else:
            self.tco = np.inf

        # Record the sigma with which it was initiated
        self.sigma = sigma1

        # Record the final stellar type
        self.kfinal = k[k>0][-1]

class InteractingBinaryStar:
    def __init__(self,
                 tphys=0.,
                 mass0_1=1.,
                 mass0_2=1.,
                 mass_1=1.,
                 mass_2=1.,
                 epoch_1=0.,
                 epoch_2=0.,
                 ospin_1=0.,
                 ospin_2=0.,
                 max_time=15000.,
                 period=100,
                 type1=1,
                 type2=1,
                 Z=ic.Z,
                 ecc=0.,
                 neta=ic.neta,
                 bwind=ic.bwind,
                 hewind=ic.hewind,
                 alpha1=ic.alpha1,
                 lambda_=ic.lamb,
                 ceflag=ic.ceflag,
                 tflag=ic.tflag,
                 ifflag=ic.ifflag,
                 wdflag=ic.wdflag,
                 bhflag=ic.bhflag,
                 nsflag=ic.nsflag,
                 piflag=ic.piflag,
                 mxns=ic.mxns,
                 idum=ic.idum,
                 pts1=ic.pts1,
                 pts2=ic.pts2,
                 pts3=ic.pts3,
                 sigma1=ic.sigma1,
                 sigma2=ic.sigma2,
                 beta=ic.beta,
                 xi=ic.xi,
                 acc2=ic.acc2,
                 epsnov=ic.epsnov,
                 eddfac=ic.eddfac,
                 gamma=ic.gamma):

        # Write input file
        with open(ic.MOBSE_DIR+'/../input/binary.input','w') as f:
            f.write(f'{mass0_1} {mass0_2} {max_time} {period} {type1} {type2} {Z} {ecc}\n')
            f.write(f'{neta} {bwind} {hewind} {alpha1} {lambda_}\n')
            f.write(f'{ceflag} {tflag} {ifflag} {wdflag} {bhflag} {nsflag} {piflag} {mxns} {idum}\n')
            f.write(f'{pts1} {pts2} {pts3}\n')
            f.write(f'{sigma1} {sigma2} {beta} {xi} {acc2} {epsnov} {eddfac} {gamma}\n')

            f.write(f'{tphys}\n')
            f.write(f'{epoch_1} {mass_1} {ospin_1}\n')
            f.write(f'{epoch_2} {mass_2} {ospin_2}\n')

            f.write('\n')

            # Write comments
            f.write('# m_zams1  m_zams2  max_time  period  type1  type2  Z  ecc\n')
            f.write('# neta  bwind  hewind  alpha1  lambda\n')
            f.write('# ceflag  tflag  ifflag  wdflag  bhflag  nsflag piflag mxns  idum\n')
            f.write('# pts1  pts2  pts3\n')
            f.write('# sigma1  sigma2  beta  xi  acc2  epsnov  eddfac  gamma\n')

        # Run MOBSE
        os.chdir(ic.MOBSE_DIR)
        os.system('./mobse.x')
        os.chdir(ic.SRC_DIR)

        # Read output
        t,k1,m0_1,m1,epoch1,ospin1,RL1,k2,m0_2,m2,epoch2,ospin2,RL2,sep,ecc = np.loadtxt(ic.MOBSE_DIR+'/../mobse.out',unpack=True,usecols=(0,1,2,3,11,12,14,15,16,17,25,26,28,30,31))

        # Detect when RLO ends again
        idx = 0
        RL1_value = RL1[0]
        RL2_value = RL2[0]
        for i in range(len(t)):
            if RL1_value >= 1 and RL1[i] < 1:
                idx = i
                break
            if RL2_value >= 1 and RL2[i] < 1:
                idx = i
                break
            if k1[i] >= 15 or k2[i] >= 15:
                idx = i
                break
            RL1_value = RL1[i]
            RL2_value = RL2[i]

        # Print empty line
        print('',end='\n\n')

        if idx == 0:
            print('MOBSE did not find a Roche lobe overflow in the bcm array. Search in bpp array...')
            t,k1,m0_1,m1,epoch1,ospin1,RL1,k2,m0_2,m2,epoch2,ospin2,RL2,sep,ecc,kw = np.loadtxt(ic.MOBSE_DIR+'/../mobse-bpp.out',unpack=True,usecols=(0,1,2,3,11,12,14,15,16,17,25,26,28,30,31,32))
            # Search for kw=4
            idx = np.where(kw==4)[0]
            if len(idx)==0:
                print('MOBSE did not find a Roche lobe overflow in the bpp array')
                self.event_status = -1
                self.t = np.nan
                self.sep = np.nan
                self.ecc = np.nan
                self.m12 = np.nan
                self.star1 = None
                self.star2 = None
                return
            else:
                idx = idx[0]
                print('MOBSE found a Roche lobe overflow in the bpp array')

        # Get final post-interaction properties
        t = t[idx]
        k1 = k1[idx]
        m0_1 = m0_1[idx]
        m1 = m1[idx]
        epoch1 = epoch1[idx]
        ospin1 = ospin1[idx]
        k2 = k2[idx]
        m0_2 = m0_2[idx]
        m2 = m2[idx]
        epoch2 = epoch2[idx]
        ospin2 = ospin2[idx]
        sep = sep[idx]
        ecc = ecc[idx]

        # Print post-interaction properties
        print('Roche lobe ends at',t,end='\n\n')
        print('Star 1:',k1,m0_1,m1,epoch1,ospin1)
        print('Star 2:',k2,m0_2,m2,epoch2,ospin2)
        print('Separation:',sep)
        print('Eccentricity:',ecc,end='\n\n')

        # Test if either star is k=15
        if k1>=15 or k2>=15:
            print('MOBSE found a merger')
            self.event_status = -1
            self.t = np.nan
            self.sep = np.nan
            self.ecc = np.nan
            self.m12 = np.nan
            self.star1 = None
            self.star2 = None
            return
            
        # Test if orbit is intact
        elif sep<=0 or ecc>=1 or ecc<0 or np.isnan(ecc) or np.isnan(sep):
            print('MOBSE found that orbit is not intact')
            self.event_status = -1
            self.t = np.nan
            self.sep = np.nan
            self.ecc = np.nan
            self.m12 = np.nan
            self.star1 = None
            self.star2 = None
            return 

        # Initiate two single stars to continue evolution
        else:
            self.event_status = 1 # Successful interaction completed

        # Prepare output
        print('Evolve single stars (ignore entries for dummy secondary)',end='\n\n')

        self.star1 = SingleStar(tphys=t,
                           mass0_1=m0_1,
                           mass_1=m1,
                           epoch_1=epoch1,
                           ospin_1=ospin1,
                           max_time=max_time,
                           type1=-int(k1),
                           Z=Z,
                           ecc=ecc,
                           neta=neta,
                           bwind=bwind,
                           hewind=hewind,
                           alpha1=alpha1,
                           lambda_=lambda_,
                           ceflag=ceflag,
                           tflag=tflag,
                           ifflag=ifflag,
                           wdflag=wdflag,
                           bhflag=bhflag,
                           nsflag=nsflag,
                           piflag=piflag,
                           mxns=mxns,
                           idum=idum,
                           pts1=pts1,
                           pts2=pts2,
                           pts3=pts3,
                           sigma1=sigma1,
                           sigma2=sigma2,
                           beta=beta,
                           xi=xi,
                           acc2=acc2,
                           epsnov=epsnov,
                           eddfac=eddfac,
                           gamma=gamma)
        
        self.star2 = SingleStar(tphys=t,
                           mass0_1=m0_2,
                           mass_1=m2,
                           epoch_1=epoch2,
                           ospin_1=ospin2,
                           max_time=max_time,
                           type1=-int(k2),
                           Z=Z,
                           ecc=ecc,
                           neta=neta,
                           bwind=bwind,
                           hewind=hewind,
                           alpha1=alpha1,
                           lambda_=lambda_,
                           ceflag=ceflag,
                           tflag=tflag,
                           ifflag=ifflag,
                           wdflag=wdflag,
                           bhflag=bhflag,
                           nsflag=nsflag,
                           piflag=piflag,
                           mxns=mxns,
                           idum=idum,
                           pts1=pts1,
                           pts2=pts2,
                           pts3=pts3,
                           sigma1=sigma1,
                           sigma2=sigma2,
                           beta=beta,
                           xi=xi,
                           acc2=acc2,
                           epsnov=epsnov,
                           eddfac=eddfac,
                           gamma=gamma)
        
        self.t = t
        self.sep = sep
        self.ecc = max(ecc,1e-3)
        self.m12 = m1+m2