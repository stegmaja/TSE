***   
      subroutine bse(m10,m20,mass1,mass2,spin1,spin2,k1,k2,z,tb,tphys,
     &     e1,a2,a1,m3,ep1,ep2,true,col,dm3,tmax,CE,KCE1,KCE2,RLO)
***   
*     
*     Evolves a binary by calling evolv2.f 
*     (see header of subroutine for algorithm description). 
*     
*     Required input is described below. 
***   
*     See Tout et al., MNRAS, 1997, 291, 732 for a description of many of the
*     processes in this code as well as the relevant references mentioned
*     within the code.
*     Updated reference is:
*     Hurley J.R., Tout C.A., & Pols O.R., 2002, MNRAS, 329, 897
*     (please use this one).
***   
*     For single star evolution see Hurley, Pols & Tout, 2000, MNRAS, 315, 543.
*     or Hurley, 2000, PhD Thesis, University of Cambridge (Chapter 2).
*     The binary evolution algorithm is described in Chapter 3 of the thesis.
***   
*     
*     B I N A R Y
*     ***********
*     
*     Roche lobe overflow.
*     --------------------
*     
*     Developed by Jarrod Hurley, IOA, Cambridge.
*     .........................................................
*     
*     Advice by Christopher Tout, Onno Pols & Sverre Aarseth.
*     ++++++++++++++++++++++++++++++++++++++++++++++++++
***   
      use my_subs
      implicit none
*     
      INCLUDE 'const_bse.h'
*     
      integer kw,kw2,kstar(2),j,k,time,k1,k2,true,true2,col,CE
      integer KCE1,KCE2,jRLO,RLO
*     
      real*8 mass0(2),mass(2),z,zpars(20),m10,m20,dm3
      real*8 epoch(2),tms(2),tphys,tphysf,dtp,aj,a1
      real*8 rad(2),lum(2),ospin(2),a2,dm1,dm2,tfin,tmax
      real*8 massc(2),radc(2),menv(2),renv(2),m3,tfin0,e1
      real*8 tb,ecc,yearsc,spin1,spin2,mass1,mass2,ep1,ep2
      PARAMETER(yearsc=3.1557d+07)
      CHARACTER*8 label(14)
      
*     
************************************************************************
*     Input:
*     
*     mass is in solar units.
*     tphysf is the maximum evolution time in Myr.
*     tb is the orbital period in days.
*     kstar is the stellar type: 0 or 1 on the ZAMS - unless in evolved state. 
*     z is metallicity in the range 0.0001 -> 0.03 where 0.02 is Population I.
*     eccentricity can be anywhere in the range 0.0 -> 1.0.
*     
*     neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
*     bwind is the binary enhanced mass loss parameter (inactive for single).
*     hewind is a helium star mass loss factor (1.0 normally).
*     alpha1 is the common-envelope efficiency parameter (1.0).  
*     lambda is the binding energy factor for common envelope evolution (0.5).
*     
*     ceflag > 0 activates spin-energy correction in common-envelope (0). #defunct#
*     ceflag = 3 activates de Kool common-envelope model (0). 
*     tflag > 0 activates tidal circularisation (1).
*     ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
*     wdflag > 0 uses modified-Mestel cooling for WDs (0). 
*     bhflag > 0 allows velocity kick at BH formation (0). 
*     nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1). 
*     mxns is the maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1). 
*     idum is the random number seed used by the kick routine. 
*     
*     Next come the parameters that determine the timesteps chosen in each
*     evolution phase:
*     pts1 - MS                  (0.05) 
*     pts2 - GB, CHeB, AGB, HeGB (0.01)
*     pts3 - HG, HeMS            (0.02)
*     as decimal fractions of the time taken in that phase.
*     
*     sigma is the dispersion in the Maxwellian for the SN kick speed (190 km/s). 
*     beta is wind velocity factor: proportional to vwind**2 (1/8). 
*     xi is the wind accretion efficiency factor (1.0). 
*     acc2 is the Bondi-Hoyle wind accretion factor (3/2). 
*     epsnov is the fraction of accreted matter retained in nova eruption (0.001). 
*     eddfac is Eddington limit factor for mass transfer (1.0).
*     gamma is the angular momentum factor for mass lost during Roche (-1.0). 
*     
*     If you enter a negative kstar then parameters for an evolved star are
*     required in the order of:
*     current age, initial mass and spin rate, 
*     otherwise the star will start on the ZAMS.
*
      true=0
      true2=0
      tfin0=tphys
      col=0
      tphysf=tmax
      mass0(1)=m10
      mass0(2)=m20
      kstar(1)=k1
      kstar(2)=k2
      neta=0.5
      bwind=0.0
      hewind=0.5
*      alpha1=1.0
      lambda=0.5
      ceflag=0
      tflag=1
      ifflag=0
      wdflag=1
      mxns=3.
      idum=29769      
      pts1=0.01
      pts2=0.01
      pts3=0.02
      beta=0.125
      xi=1.0
      acc2=1.5
      epsnov=0.01
      eddfac=1.
      gamma=-1.0
      
      epoch(1)=ep1
      epoch(2)=ep2
      
      mass(1)=mass1
      mass(2)=mass2
      ospin(1) = spin1
      ospin(2) = spin2
      ecc=e1
      
c     tphys=0
*     
*     Initialize the parameters.
*     Set the initial spin of the stars. If ospin is zero (actually < 0.001)
*     at time zero then evolv2 will set an appropriate ZAMS spin. If 
*     ospin is greater than zero then it will start with that spin regardless
*     of the time. If you want to start at time zero with negligible spin 
*     then I suggest using a negligible value (but greater than 0.001).
*     If ospin is negative then the stars will be in co-rotation with the orbit.
*     
      
      if(idum.gt.0) idum = -idum
*      CLOSE(22)
c      WRITE(*,*) 
*     
*     Note that this routine can be used to evolve a single star if you 
*     simply set mass(2) = 0.0 or tb = 0.0 (setting both is advised as  
*     well as some dummy value for ecc). 
*     
************************************************************************
*     
*     Set parameters which depend on the metallicity 
*     
      CALL zcnsts(z,zpars)
*     
*     Set the collision matrix.
*     
      CALL instar
*     
      label(1) = 'INITIAL '
      label(2) = 'KW CHNGE'
      label(3) = 'BEG RCHE'
      label(4) = 'END RCHE'
      label(5) = 'CONTACT '
      label(6) = 'COALESCE'
      label(7) = 'COMENV  '
      label(8) = 'GNTAGE  '
      label(9) = 'NO REMNT'
      label(10) = 'MAX TIME'
      label(11) = 'DISRUPT '
      label(12) = 'BEG SYMB'
      label(13) = 'END SYMB'
      label(14) = 'BEG BSS'
*     
*     Set the data-save parameter. If dtp is zero then the parameters of the 
*     star will be stored in the bcm array at each timestep otherwise they 
*     will be stored at intervals of dtp. Setting dtp equal to tphysf will 
*     store data only at the start and end while a value of dtp greater than 
*     tphysf will mean that no data is stored.
*     
      dtp=0
*     
*     Evolve the binary.
*      
      CALL evolv2(kstar,mass0,mass,rad,lum,massc,radc,
     &     menv,renv,ospin,epoch,tms,
     &     tphys,tphysf,dtp,z,zpars,tb,ecc)
*
      tphys=tfin0
      
************************************************************************
*     Output:
*     First check that bcm is not empty and that RLO occurred.
*     
      if(bcm(1,1).lt.0.0) goto 50
*     
*     The bpp array acts as a log, storing parameters at each change
*     of evolution stage.
*     
 50   j = 0
      WRITE(*,*)'     TIME      M1       M2   K1 K2        SEP    ECC',  
     &     '  R1     R2      TYPE'
 52   j = j + 1 
      if(bpp(j,1).lt.0.0) goto 60
      kstar(1) = INT(bpp(j,4))
      kstar(2) = INT(bpp(j,5))
      kw = INT(bpp(j,10))     
      WRITE(*,100)(bpp(j,k),k=1,3),kstar,(bpp(j,k),k=6,9),label(kw)

      if(bpp(j,2).lt.0.or.bpp(j,3).lt.0)kw=6

      if(a2.eq.-1.d0)then ! BSE is called when there is no longer a tertiary
            if(kw.eq.3)then ! RLO actually starts
                  if((kstar(1).eq.13.or.kstar(1).eq.14).and.
     &                  (kstar(2).eq.13.or.kstar(2).eq.14))then ! CO-CO merger
                        col=3
                        tfin=bpp(j-1,1)
                        goto 60
                  endif
                  true = 1
                  jRLO = j
                  RLO = 1
            endif
c            if(kw.eq.4)then ! RLO finishes
c                  tfin=bpp(j,1)
c                  goto 60
c            endif
            if(kw.eq.7)then ! CE occurs
                  CE=1 
                  jRLO=INT(MIN(jRLO,j))
                  if(bpp(jRLO,8).gt.bpp(jRLO,6)*(1.d0-bpp(jRLO,7))
     &            *Roche(bpp(jRLO,2),bpp(jRLO,3)))then ! Primary was donor
                        if(KCE1.ne.2)KCE1=INT(bpp(jRLO,4)) ! ensure that previous HG donor is not overwritten
                        print*,'Unstable donor type',INT(bpp(jRLO,4))
                  else
                        if(KCE2.ne.2)KCE2=INT(bpp(jRLO,5)) ! ensure that previous HG donor is not overwritten
                        print*,'Unstable donor type',INT(bpp(jRLO,5))
                  endif
                  if(kstar(1).eq.15.or.kstar(2).eq.15)then ! ! Stellar coaelscence during CE
                        col=1
                        tfin=bpp(j,1)
                        goto 60
                  endif
            endif
            if(kw.eq.5.or.kw.eq.6)then ! Stellar coaelscence
                  col=1
                  tfin=bpp(j,1)
                  goto 60
            endif 
            if(kw.eq.11)then ! Disruption
                  col=2
                  tfin=bpp(j,1)
                  goto 60
            endif 
            if(kw.eq.10)then
                  if((kstar(1).eq.13.or.kstar(1).eq.14).and.
     &                  (kstar(2).eq.13.or.kstar(2).eq.14).and.
     &                  bpp(j,6)*0.00465d0.lt.1.d-3) col=3 ! CO-CO merger
                  tfin=bpp(j,1)
            endif
      else
            if(kw.eq.3)then ! RLO actually starts
                  if((kstar(1).eq.13.or.kstar(1).eq.14).and.
     &                  (kstar(2).eq.13.or.kstar(2).eq.14))then ! CO-CO merger
                        col=3
                        tfin=bpp(j-1,1)
                        goto 60
                  endif
                  true = 1
                  jRLO = j
                  RLO = 1
            endif
            if(kw.eq.4)then ! RLO finishes
                  tfin=bpp(j,1)
                  goto 60
            endif
            if(kw.eq.7)then ! CE occurs
                  CE=1 
                  jRLO=INT(MIN(jRLO,j))
                  if(bpp(jRLO,8).gt.bpp(jRLO,6)*(1.d0-bpp(jRLO,7))
     &            *Roche(bpp(jRLO,2),bpp(jRLO,3)))then ! Primary was donor
                        if(KCE1.ne.2)KCE1=INT(bpp(jRLO,4)) ! ensure that previous HG donor is not overwritten
                        print*,'Unstable donor type',INT(bpp(jRLO,4))
                  else
                        if(KCE2.ne.2)KCE2=INT(bpp(jRLO,5)) ! ensure that previous HG donor is not overwritten
                        print*,'Unstable donor type',INT(bpp(jRLO,5))
                  endif
                  if(kstar(1).eq.15.or.kstar(2).eq.15)then ! ! Stellar coaelscence during CE
                        col=1
                        tfin=bpp(j,1)
                        goto 60
                  endif
            endif
            if(kw.eq.5.or.kw.eq.6)then ! Stellar coaelscence
                  col=1
                  tfin=bpp(j,1)
                  goto 60
            endif 
            if(kw.eq.11)then ! Disruption
                  col=2
                  tfin=bpp(j,1)
                  goto 60
            endif 
      endif

      goto 52
 60   continue

*     The bcm array stores the stellar and orbital parameters at the 
*     specified output times. The parameters are (in order of storage):
*     
*     Time, 
*     [stellar type, initial mass, current mass, log10(L), log10(r),
*     log10(Teff), core mass, core radius, mass of any convective 
*     envelope, radius of the envelope, epoch, spin, mass loss rate and 
*     ratio of radius to roche lobe radius (repeated for secondary)],
*     period, separation, eccentricity.
*
      if(true.eq.1.or.a2.eq.-1.d0)then
        
c     OPEN(23,file='binary.dat', status='unknown')
         j = 0
 30      j = j + 1
c         open(44,file="evolution.dat",action='write',access='append')
c            write(44,*)bcm(j,1),bcm(j,4),bcm(j,18),m3,e1,'NaN',a1,a2,
c     &bcm(j,15)*0.00465d0,bcm(j,29)*0.00465d0
c            close(44)
         if(bcm(j,1).lt.0.0)then
            bcm(j-1,1) = bcm(j,1)
            j = j - 1
         endif
         kw = INT(bcm(j,2))
         kw2 = INT(bcm(j,16))

c         WRITE(23,99)bcm(j,1),kw,kw2,bcm(j,4),bcm(j,18),
c     &        bcm(j,8),bcm(j,22), 
c     &        bcm(j,6),bcm(j,20),bcm(j,15),bcm(j,29),
c     &        bcm(j,5),bcm(j,19),bcm(j,13),bcm(j,27),
c     &        bcm(j,26),bcm(j,12),
c     &        bcm(j,31),bcm(j,32)
         dm1=(bcm(j-1,4)-bcm(j,4))
         dm2=(bcm(j-1,18)-bcm(j,18))
         if(bcm(j-1,4).eq.0)dm1=0.d0
         if(bcm(j-1,18).eq.0)dm2=0.d0
         
         if(bcm(j-1,4).ge.0.d0.and.bcm(j-1,18).ge.0.d0.and.
     &      bcm(j,1).ge.0.0d0.and.bcm(j,1).le.tfin)then
            if(bcm(j,32).lt.0)then
              print*,'Ecc time',bcm(j,1),bcm(j,32),bcm(j,31)
              bcm(j,32)=10.d0
            endif
            if(a2.ge.0.0)a2=a2+a2*(dm1+dm2+dm3)/(bcm(j,4)+bcm(j,18)+m3)
            a1=bcm(j,31)*0.00465d0 ! in AU
            e1=max(bcm(j,32),1.d-3)           
            m10=bcm(j,3)
            m20=bcm(j,17)
            mass1=bcm(j,4)
            mass2=bcm(j,18)
            spin1=bcm(j,13)
            spin2=bcm(j,27)
            tphys=bcm(j,1)
            if(tphys.lt.0)tphys=tphysf            
            ep1=bcm(j,12)
            ep2=bcm(j,26)
            k1=kw
            k2=kw2
            goto 30
         endif
c         if(bcm(j,1).ge.0.0.and.bcm(j,1).le.tfin) goto 30
c     CLOSE(23)
      endif
      print*,'T M1 M2 SEP ECC after BSE.F ',bcm(j,1),mass1,mass2,a1,e1
 99   FORMAT(f10.4,2i3,10f10.4,5e12.4,f7.3)
 999  FORMAT(f10.4,2f10.4,1p,2e12.4)

 100      FORMAT(f11.4,2f9.3,2i3,f13.3,f6.2,2f9.3,2x,a8)
c      WRITE(*,*)
*     
************************************************************************
*
      RETURN
      END SUBROUTINE bse
***   
      
