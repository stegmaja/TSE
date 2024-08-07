***
      PROGRAM mobse
***
*
* Evolves a binary by calling evolv3.f 
* (see header of subroutine for algorithm description). 
*
* Required input is described below. 
***
*           Massive Objects
*           *************** 
* See Giacobbo et al. MNRAS, 2018, 474 for a description of the upgrades.
*** 
* See Tout et al., MNRAS, 1997, 291, 732 for a description of many of the
* processes in this code as well as the relevant references mentioned
* within the code.
* Updated reference is:
*           Hurley J.R., Tout C.A., & Pols O.R., 2002, MNRAS, 329, 897
***
* For single star evolution see Hurley, Pols & Tout, 2000, MNRAS, 315, 543.
* or Hurley, 2000, PhD Thesis, University of Cambridge (Chapter 2).
* The binary evolution algorithm is described in Chapter 3 of the thesis.
***
*
*           B I N A R Y
*           ***********
*
*       Roche lobe overflow.
*       --------------------
*
*       Developed by Jarrod Hurley, IOA, Cambridge.
*       .........................................................
*
*       Advice by Christopher Tout, Onno Pols & Sverre Aarseth.
*       ++++++++++++++++++++++++++++++++++++++++++++++++++
***
      implicit none
*
      INCLUDE '../input/const_mobse.h'
*
      integer kw,kw2,kstar(2),j,k
*
      real*8 mass0(2),mass(2),z,zpars(20)
      real*8 epoch(2),tms(2),tphys,tphysf,dtp,aj
      real*8 rad(2),lum(2),ospin(2)
      real*8 massc(2),radc(2),menv(2),renv(2)
      real*8 tb,ecc,yearsc
      PARAMETER(yearsc=3.1557d+07)
      CHARACTER*8 label(14)
      character(len=255) :: arg1, arg2, arg3

      real*8 ffb
      COMMON /KICKSN/ ffb
*
************************************************************************
* Input:
*
* mass is in solar units from 0.1 -> 150 Msun.
* tphysf is the maximum evolution time in Myr.
* tb is the orbital period in days.
* kstar is the stellar type: 0 or 1 on the ZAMS - unless in evolved state. 
* z is metallicity in the range 0.0001 -> 0.03 where 0.02 is Population I.
* eccentricity can be anywhere in the range 0.0 -> 1.0.
*
* neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
* bwind is the binary enhanced mass loss parameter (inactive for single).
* hewind is a helium star mass loss factor (normally inactive).
* alpha1 is the common-envelope efficiency parameter (1.0).  
* lambda is the binding energy factor for common envelope evolution (0.1).
*
* ceflag > 0 activates spin-energy correction in common-envelope (0). #defunct#
* ceflag = 3 activates de Kool common-envelope model (0). 
* tflag > 0 activates tidal circularisation (1).
* ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
* wdflag > 0 uses modified-Mestel cooling for WDs (1). 
* bhflag > 0 allows velocity kick at BH formation (3). 
*     1 --> old case
*     2 --> fallback case: Fryer et al. 2012, ApJ, 749, 91
*     3 --> mass ejected case: Giacobbo & Mapelli 2020, ApJ, 891
*     4 --> full kick case
* nsflag > 0 takes NS/BH mass from (default 3): 
*     1 --> Belczynski et al. 2002, ApJ, 572, 407
*     2 --> Rapid supernova model Fryer et al. 2012, ApJ, 749, 91
*     3 --> Delayed supernova model Fryer et al. 2012, ApJ, 749, 91
* piflag > 0 activates the PPISNs and PISNe (1) Spera et al. 2017, MNRAS, 470, 4739.  
* mxns is the maximum NS mass (3.0).
* idum is the random number seed used by the kick routine. 
*
* Next come the parameters that determine the timesteps chosen in each
* evolution phase:
*                 pts1 - MS                  (0.05) 
*                 pts2 - GB, CHeB, AGB, HeGB (0.01)
*                 pts3 - HG, HeMS            (0.02)
* as decimal fractions of the time taken in that phase.
*
* sigma1 is the dispersion in the Maxwellian for the SN kick speed (265. km/s)
*            from Hobbs et al. 2005, ApJ, 591, 288. 
* sigma2 is the dispersion in the Maxwellian for the SN kick speed (15. km/s)
*            to considere the different mechanims involve in ECS. 
* beta is wind velocity factor: proportional to vwind**2 (1/8). 
* xi is the wind accretion efficiency factor (1.0). 
* acc2 is the Bondi-Hoyle wind accretion factor (3/2). 
* epsnov is the fraction of accreted matter retained in nova eruption (0.001). 
* eddfac is Eddington limit factor for mass transfer (1.0).
* gamma is the angular momentum factor for mass lost during Roche (-1.0). 
*
* If you enter a negative kstar then parameters for an evolved star are
* required in the order of:
* current age, initial mass and spin rate, 
* otherwise the star will start on the ZAMS.
*
* Get filenames
      call get_command_argument(1, arg1)
      call get_command_argument(2, arg2)
      call get_command_argument(3, arg3)
* Trim the filenames
      arg1 = trim(adjustl(arg1))
      arg2 = trim(adjustl(arg2))
      arg3 = trim(adjustl(arg3))
*
*      OPEN(22,file='../input/binary.input', status='old')
      OPEN(22,file=arg1, status='old')
      READ(22,*)mass0(1),mass0(2),tphysf,tb,kstar(1),kstar(2),z,ecc
      READ(22,*)neta,bwind,hewind,alpha1,lambda
      READ(22,*)ceflag,tflag,ifflag,wdflag,bhflag,nsflag,piflag,
     &           mxns,idum
      READ(22,*)pts1,pts2,pts3
      READ(22,*)sigma1,sigma2,beta,xi,acc2,epsnov,eddfac,gamma
      if(kstar(1).lt.0.or.kstar(2).lt.0)then
         READ(22,*)tphys
*         READ(22,*)aj,mass(1),ospin(1)
*         epoch(1) = tphys - aj
         READ(22,*)epoch(1),mass(1),ospin(1)
         kstar(1) = ABS(kstar(1))
*         READ(22,*)aj,mass(2),ospin(2)
*         epoch(2) = tphys - aj
         READ(22,*)epoch(2),mass(2),ospin(2)
         kstar(2) = ABS(kstar(2))
      else
*
* Initialize the parameters.
* Set the initial spin of the stars. If ospin is zero (actually < 0.001)
* at time zero then evolv3 will set an appropriate ZAMS spin. If 
* ospin is greater than zero then it will start with that spin regardless
* of the time. If you want to start at time zero with negligible spin 
* then I suggest using a negligible value (but greater than 0.001).
* If ospin is negative then the stars will be in co-rotation with the orbit.
*
         tphys = 0.d0
         mass(1) = mass0(1)
         epoch(1) = 0.d0
         ospin(1) = 0.d0
         mass(2) = mass0(2)
         epoch(2) = 0.d0
         ospin(2) = 0.d0
      endif
      if(idum.gt.0) idum = -idum
      CLOSE(22)
      WRITE(*,*)
*
* Note that this routine can be used to evolve a single star if you 
* simply set mass(2) = 0.0 or tb = 0.0 (setting both is advised as  
* well as some dummy value for ecc). 
*
************************************************************************
*
* Set parameters which depend on the metallicity 
*
      CALL zcnsts(z,zpars)
*
* Set the collision matrix.
*
      CALL instar
*
* Set the array with the labels.
*
      label(1) = 'INITIAL '
      label(2) = 'KW_CHNGE'
      label(3) = 'BEG_RCHE'
      label(4) = 'END_RCHE'
      label(5) = 'CONTACT '
      label(6) = 'COELESCE'
      label(7) = 'COMENV  '
      label(8) = 'GNTAGE  '
      label(9) = 'NO_REMNT'
      label(10) = 'MAX_TIME'
      label(11) = 'DISRUPT '
      label(12) = 'BEG_SYMB'
      label(13) = 'END_SYMB'
      label(14) = 'BEG_BSS'
*
* Set the data-save parameter. If dtp is zero then the parameters of the 
* star will be stored in the bcm array at each timestep otherwise they 
* will be stored at intervals of dtp. Setting dtp equal to tphysf will 
* store data only at the start and end while a value of dtp greater than 
* tphysf will mean that no data is stored.
*
      dtp = (tphysf-tphys)/45000.d0
      dtp = 0.d0
*
* Evolve the binary.
* 
      CALL evolve(kstar,mass0,mass,rad,lum,massc,radc,
     &            menv,renv,ospin,epoch,tms,
     &            tphys,tphysf,dtp,z,zpars,tb,ecc)
*
************************************************************************
* Output:
* First check that bcm is not empty.
*
      if(bcm(1,1).lt.0.0) goto 50
* The bcm array stores the stellar and orbital parameters at the 
* specified output times. The parameters are (in order of storage in bcm):
*
*    Time, 
*    [stellar type, initial mass, current mass, log10(L), log10(r),
*    log10(T), core mass, core radius, mass of any convective 
*    envelope, radius of the envelope, epoch, spin, mass loss rate and 
*    ratio of radius to roche lobe radius (repeated for secondary)],
*    period[year], separation, eccentricity.
*
*      OPEN(24,file='../mobse.out',status='unknown')
      OPEN(24,file=arg2,status='unknown')
      j = 0
 30   j = j + 1
      kw = INT(bcm(j,2))
      kw2 = INT(bcm(j,16))
      WRITE(24,*)bcm(j,1),bcm(j,2),bcm(j,3),bcm(j,4),bcm(j,5),
     &          bcm(j,6),bcm(j,7),bcm(j,8),bcm(j,9),bcm(j,10),
     &          bcm(j,11),bcm(j,12),bcm(j,13),bcm(j,14),bcm(j,15),
     &          bcm(j,16),bcm(j,17),bcm(j,18),bcm(j,19),bcm(j,20),
     &          bcm(j,21),bcm(j,22),bcm(j,23),bcm(j,24),bcm(j,25),
     &          bcm(j,26),bcm(j,27),bcm(j,28),bcm(j,29),bcm(j,30),
     &          bcm(j,31),bcm(j,32),bcm(j,33),bcm(j,34),ffb

      if(bcm(j,1).ge.0.0) goto 30
      CLOSE(24)
* 99   FORMAT(f10.4,i3,5f10.4,i3,5f10.4,2f12.4,f7.3,a8)

*      OPEN(25,file='../mobse-bpp.out',status='unknown')
      OPEN(25,file=arg3,status='unknown')
      j = 0
 31   j = j + 1
      kw = INT(bpp(j,2))
      kw2 = INT(bpp(j,16))
      WRITE(25,*)bpp(j,1),bpp(j,2),bpp(j,3),bpp(j,4),bpp(j,5),
     &          bpp(j,6),bpp(j,7),bpp(j,8),bpp(j,9),bpp(j,10),
     &          bpp(j,11),bpp(j,12),bpp(j,13),bpp(j,14),bpp(j,15),
     &          bpp(j,16),bpp(j,17),bpp(j,18),bpp(j,19),bpp(j,20),
     &          bpp(j,21),bpp(j,22),bpp(j,23),bpp(j,24),bpp(j,25),
     &          bpp(j,26),bpp(j,27),bpp(j,28),bpp(j,29),bpp(j,30),
     &          bpp(j,31),bpp(j,32),bpp(j,33),ffb

      if(bpp(j,1).ge.0.0) goto 31
      CLOSE(25)
 99   FORMAT(f10.4,i3,5f10.4,i3,5f10.4,2f12.4,f7.3,a8)
*
* The bpp array acts as a log, storing parameters at each change
* of evolution stage (it has the same information of bcm plus the 
* labels at the end).
*
 50   j = 0
      WRITE(*,*)'     TIME      M1       M2   K1 K2        SEP    ECC',
     &          '  R1/ROL1 R2/ROL2  TYPE'
 52   j = j + 1
      if(bpp(j,1).lt.0.0) goto 60
         kstar(1) = INT(bpp(j,2))
         kstar(2) = INT(bpp(j,16))
         kw = INT(bpp(j,33))
         WRITE(*,100)bpp(j,1),bpp(j,4),bpp(j,18),kstar,bpp(j,31),
     & bpp(j,32),bpp(j,15),bpp(j,29),label(kw)

      goto 52

 60   continue
 100   FORMAT(f11.4,2f11.5,2i3,f15.5,f6.2,2f8.3,2x,a8)
      WRITE(*,*)'Fallback:',ffb
*
************************************************************************
*
      STOP
      END
***
