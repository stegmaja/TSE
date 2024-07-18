import numpy as np
from scipy import interpolate
import os
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
        with open(ic.MOBSE_DIR+'/'+ic.mobse_input,'w') as f:
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
        os.system('./mobse.x '+ic.mobse_input+' '+ic.mobse_output+' '+ic.mobse_log)
        os.chdir(ic.SRC_DIR)

        # Read output
        t,k,m0,m,loglum,logr,massc,radc,menv,renv,epoch,ospin,dm,ffb = np.loadtxt(ic.MOBSE_DIR+'/'+ic.mobse_output,unpack=True,usecols=(0,1,2,3,4,5,7,8,9,10,11,12,13,-1))

        # Only keep values where t is unique
        t, idx = np.unique(t,return_index=True)
        k = k[idx]
        m0 = m0[idx]
        m = m[idx]
        loglum = loglum[idx]
        logr = logr[idx]
        massc = massc[idx]
        radc = radc[idx]
        menv = menv[idx]
        renv = renv[idx]
        epoch = epoch[idx]
        ospin = ospin[idx]
        dm = dm[idx]
        self.ffb = ffb[0]

        # Write interpolators in time
        self.interpolators = {}
        self.interpolators['k'] = interpolate.interp1d(t,k,fill_value=(k[0],k[-1]),bounds_error=False,kind='previous')
        self.interpolators['m0'] = interpolate.interp1d(t,m0,fill_value=(m0[0],m0[-1]),bounds_error=False)
        self.interpolators['m'] = interpolate.interp1d(t,m,fill_value=(m[0],m[-1]),bounds_error=False)
        self.interpolators['logr'] = interpolate.interp1d(t,logr,fill_value=(logr[0],logr[-1]),bounds_error=False)
        self.interpolators['epoch'] = interpolate.interp1d(t,epoch,fill_value=(epoch[0],epoch[-1]),bounds_error=False)
        self.interpolators['ospin'] = interpolate.interp1d(t,ospin,fill_value=(ospin[0],ospin[-1]),bounds_error=False)
        #self.interpolators['dm'] = interpolate.interp1d(t[:-1],np.diff(m)/(np.diff(t)),bounds_error=False,fill_value=0.)
        self.interpolators['dm'] = interpolate.interp1d(t,dm,fill_value=(dm[0],0),bounds_error=False)
        self.interpolators['loglum'] = interpolate.interp1d(t,loglum,fill_value=(loglum[0],loglum[-1]),bounds_error=False)
        self.interpolators['menv'] = interpolate.interp1d(t,menv,fill_value=(menv[0],menv[-1]),bounds_error=False)
        self.interpolators['renv'] = interpolate.interp1d(t,renv,fill_value=(renv[0],renv[-1]),bounds_error=False)
        self.interpolators['massc'] = interpolate.interp1d(t,massc,fill_value=(massc[0],massc[-1]),bounds_error=False)
        self.interpolators['radc'] = interpolate.interp1d(t,radc,fill_value=(radc[0],radc[-1]),bounds_error=False)
        
        # Test if star is CO initially
        if k[0]>=13 or (k[0]==0 and k[1]>=13):
            self.initial_CO = True
        if k[0]>=13:
            self.initial_CO = True
        else:
            self.initial_CO = False
        
        # Get time of CO formation
        if np.any(k>=13) and np.any(k<=12):
            self.tco = t[np.where(k>=13)[0][0]]
            self.m_SN = m[(k<=12)][-1]
            self.dm_SN = m[(k<=12)][-1] - m[k>=13][0]
        else:
            self.tco = np.inf

        if self.initial_CO:
            self.tco = np.inf


        # Record the sigma with which it was initiated
        self.sigma = sigma1

        # Record the final stellar type
        self.kfinal = k[k>=0][-1]

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
                 gamma=ic.gamma,
                 just_print=False):

        # Write input file
        with open(ic.MOBSE_DIR+'/'+ic.mobse_input,'w') as f:
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
        os.system('./mobse.x '+ic.mobse_input+' '+ic.mobse_output+' '+ic.mobse_log)
        os.chdir(ic.SRC_DIR)

        # If you don't want to do anything but just let mobse print the output (intendend for comparison with isolated binaries)
        if just_print:
            os.chdir(ic.MOBSE_DIR)
            os.system('cp '+ic.mobse_output+' '+(ic.mobse_output).replace('.out','_isolated.out'))
            os.system('cp '+ic.mobse_log+' '+(ic.mobse_log).replace('.out','_isolated.out'))
            os.chdir(ic.SRC_DIR)
            return

        # Read output
        t,k1,m0_1,m1,epoch1,ospin1,RL1,k2,m0_2,m2,epoch2,ospin2,RL2,sep,ecc = np.loadtxt(ic.MOBSE_DIR+'/'+ic.mobse_output,unpack=True,usecols=(0,1,2,3,11,12,14,15,16,17,25,26,28,30,31))
        t_bpp,k1_bpp,m0_1_bpp,m1_bpp,epoch1_bpp,ospin1_bpp,RL1_bpp,k2_bpp,m0_2_bpp,m2_bpp,epoch2_bpp,ospin2_bpp,RL2_bpp,sep_bpp,ecc_bpp,kw_bpp = np.loadtxt(ic.MOBSE_DIR+'/'+ic.mobse_log,unpack=True,usecols=(0,1,2,3,11,12,14,15,16,17,25,26,28,30,31,32))

        # Test if RLO within first Myr, and determine indices
        RL1_roll = np.roll(RL1,1)
        idx1_start = np.where((RL1>1) & (t-t[0]<1))[0]
        idx1_end = np.where((RL1<1) & (RL1_roll>=1) & (t>=0))[0]

        RL2_roll = np.roll(RL2,1)
        idx2_start = np.where((RL2>1) & (t-t[0]<1))[0]
        idx2_end = np.where((RL2<1) & (RL2_roll>=1) & (t>=0))[0]

        idx1_start_bpp = np.where((kw_bpp==3) & (t_bpp-t_bpp[0]<1))[0]
        idx1_end_bpp = np.where((kw_bpp==4) & (t_bpp>=0))[0]

        idx2_start_bpp = np.where((kw_bpp==3) & (t_bpp-t_bpp[0]<1))[0]
        idx2_end_bpp = np.where((kw_bpp==4) & (t_bpp>=0))[0]

        if len(idx1_start)>0 or len(idx1_end)>0:
            # Primary RLO detected
            idx = idx1_end[0]
        elif len(idx2_start)>0 or len(idx2_end)>0:
            # Secondary RLO detected
            idx = idx2_end[0]
        elif len(idx1_start_bpp)>0 or len(idx1_end_bpp)>0:
            print('MOBSE did not find a Roche lobe overflow in the bcm array. Search in bpp array...')
            idx = idx1_end_bpp[0]
            t,k1,m0_1,m1,epoch1,ospin1,RL1,k2,m0_2,m2,epoch2,ospin2,RL2,sep,ecc = t_bpp,k1_bpp,m0_1_bpp,m1_bpp,epoch1_bpp,ospin1_bpp,RL1_bpp,k2_bpp,m0_2_bpp,m2_bpp,epoch2_bpp,ospin2_bpp,RL2_bpp,sep_bpp,ecc_bpp
        elif len(idx2_start_bpp)>0 or len(idx2_end_bpp)>0:
            print('MOBSE did not find a Roche lobe overflow in the bcm array. Search in bpp array...')
            idx = idx2_end_bpp[0]
            t,k1,m0_1,m1,epoch1,ospin1,RL1,k2,m0_2,m2,epoch2,ospin2,RL2,sep,ecc = t_bpp,k1_bpp,m0_1_bpp,m1_bpp,epoch1_bpp,ospin1_bpp,RL1_bpp,k2_bpp,m0_2_bpp,m2_bpp,epoch2_bpp,ospin2_bpp,RL2_bpp,sep_bpp,ecc_bpp
        else:
            print('MOBSE did not find a Roche lobe overflow in the bcm or bpp array')
            self.event_status = -1

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

        # Check if the binary underwent a CE
        if 7 in kw_bpp[:idx+1]:
            ic.CE = True
            print('MOBSE found a common envelope phase')
        else:
            print('MOBSE found no common envelope phase')

        # Test if either star is k=15
        if k1>=15 or k2>=15:
            print('MOBSE found a merger')
            self.event_status = -2
            
        # Test if orbit is intact
        elif sep<=0 or ecc>=1 or ecc<0 or np.isnan(ecc) or np.isnan(sep):
            print('MOBSE found that orbit is not intact')
            self.event_status = -3

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
                           max_time=ic.max_time,
                           type1=-int(k1),
                           Z=Z,
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
                           max_time=ic.max_time,
                           type1=-int(k2),
                           Z=Z,
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
        self.ospin1 = ospin1
        self.ospin2 = ospin2
        self.k1 = k1
        self.k2 = k2