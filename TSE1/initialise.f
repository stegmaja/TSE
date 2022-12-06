*
************************************************************************
*
*     Main Code
*     Variable step-size Runge-Kutta-Fehlberg method 7(8) is used
*     Kozai with stellar spin evolution  
*
*     Authors: Jakob Stegmann and Fabio Antonini, Cardiff 2020
*       
      program spinTriple     
      use my_subs
      use commonV
      use tables
      implicit none
      external sse,vrotf,bse
      integer nq,j,h,N,iterator
      integer bhflag,kw1,kw2,true,col,CE,KCE1,KCE2,RLO
      integer kw
      integer,dimension(8) :: time_values
      parameter (nq=20,N=1000000)
      real*8 vrotf,dms1,dms2,dms3,vs(3),jorb
      real*8 vstore1,vstore2,vstore3
      real*8 y(nq),tol(nq),S1(3),S2(3),Rcr
      real*8 e1,e2,a1,a2,t,dt
      real*8 DT0,z,sigma
      real*8 tmax,tmax1,tmax2,tmax3
      real*8 j1,j2,Om1,Om2,Sp1,Sp2
      real*8 SevalDouble,MA
      real*8 fb1,fb2,fb3
      real*8 tRmax3,tRmax
      real*8 harv
      real*8 i1,i2
      real*8 an1,an2
      real*8 o1,o2
      real*8 Pol_S1G1,Azi_S1G1,Pol_S2G1,Azi_S2G1
      real*8 Prd1,Prd2
      real*8 logPrd1,logPrd2
      real*8 q,ecc
      real*8 G,low_mass
      real*8 power_law,logp
      real scm(1000000,14),spp(20,3),scm1(1000000,14),scm2(1000000,14)
      real scm3(1000000,14)
      real*8 m1,m2,m3,R1,R2,R3
      real*8 dtp,mco0
      integer ceflag,tflag,ifflag,nsflag,wdflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag
      real*8 alpha1,lambda
      COMMON /VALUE2/ alpha1,lambda
      COMMON /VALUE4/ sigma,bhflag
      COMMON /SINGLE/ scm,spp
      data tol/nq*1d-12/
      print*,'--------------------------------------------------------'
*
************************************************************************
*
*     Randomise initial seed, see sub.f
*
      call init_random_seed()
*
************************************************************************
*
*     Constants and units
*       
      pi=2.d0*ASIN(1.d0)
      G=2.959d-4 ! in AU^3/Solar mass/day^2
      c=1.d4                    !units M_Sun, AU, G=1
      ut=58.d0/365.d0           !to convert time in years
      RSun=0.00465d0            !R_Sun in AU
      low_mass=8.d0
      dtp=0.01d0
*
************************************************************************
*
*     Parameters from input file
*  
      open(unit=99,file="input")
      read(99,*)
      read(99,*)
      read(99,*)sigma,tmax,tmax1,tmax2,tmax3
      read(99,*)
      read(99,*)K1,K2,KL1,KL2,Kq1,Kq2
      close(99)
*
************************************************************************
*
*     Metallicity
*  
      CALL GETARG(1,arg1)
      read(arg1,*)z

      print*,z
*
************************************************************************
*
*     iterator
*  
      CALL GETARG(2,arg5)
      read(arg5,*)iterator
      print*,'IC line',iterator
      h=0
*
************************************************************************
*
*     Primary mass
*
 1    continue
      M1=power_law(-2.3d0,low_mass,100.d0)
      logPrd1 = logp(m1,0.2d0,5.5d0)
      Prd1 = 10.d0**logPrd1
      m2 = m1*q(m1,logPrd1)
      e1 = ecc(m1,logPrd1)
      if(M1.lt.low_mass.or.M2.lt.low_mass)then
            goto 1
      else
            a1=(Prd1*sqrt(G*(M1+M2))/(2.d0*pi))**(2.d0/3.d0) ! in AU
      endif
*
************************************************************************
*
*           Compute stellar evolution using sse
*    
      print*,'evolve star 1'
      call sse(M1,N,tabt1,tabm1,tabr1,tabdm1,tabmenv1,tabrenv1,
     &tablumin1,kwv1,j,tco1,dms1,Rmax1,fb1,z,tmax1,kw,tRmax,mco0,dtp)
      kwfin1 = kw
      scm1=scm
      cont1=j-1
c      CALL FMMsplineDouble(tabt1,tabdm1,cont1)
      print*,'evolve star 2'
      call sse(M2,N,tabt2,tabm2,tabr2,tabdm2,tabmenv2,tabrenv2,
     &tablumin2,kwv2,j,tco2,dms2,Rmax2,fb2,z,tmax2,kw,tRmax,mco0,dtp)
      kwfin2 = kw
      scm2=scm      
      cont2=j-1     
c      CALL FMMsplineDouble(tabt2,tabr2,cont2)

      R1=Rsun*10.d0**SevalDouble(0.0d0,tabt1,tabr1,cont1) 
      R2=Rsun*10.d0**SevalDouble(0.0d0,tabt2,tabr2,cont2) 

      if(R1.gt.a1*(1.d0-e1)*Roche(m1,m2))then   
            print*,'Initial condition unrealistic (Roche comp 1)'
            print*,R1,a1,m1
            goto 1
      elseif(R2.gt.a1*(1.d0-e1)*Roche(m2,m1))then  
            print*,'Initial condition unrealistic (Roche comp 2)'
            print*,R2,a1,m1
            goto 1
      endif

      call random_number(harv)
      if(harv.ge.1.d0/3.d0)then ! Massive inner binary member
 2          continue
            logPrd2 = logp(m1,0.2d0,8.d0)
            Prd2 = 10.d0**logPrd2
            m3 = m1*q(m1,logPrd2)
            e2 = ecc(m1,logPrd2)
            a2=(Prd2*sqrt(G*(M1+M2+M3))/(2.d0*pi))**(2.d0/3.d0)

            MA=2.8d0*((1.d0+M3/(M1+M2))*(1.d0+e2)/sqrt(1.d0-e2))
     &      **(2.d0/5.d0) !Mardling+Aarseth   
            if(e1.gt.1.d0.or.e2.gt.1.d0.or.a1.lt.0.d0
     &      .or.a2.lt.0.d0)then
               print*,'Initial condition unrealistic (apart)'
               goto 2
            elseif(a2*(1.d0-e2).lt.MA*a1)then
               print*,'Initial condition unrealistic (unstable)'
               goto 2
            endif
*
************************************************************************
*
*           Compute stellar evolution using sse
*       
            print*,'evolve Tertiary'
            call sse(M3,N,tabt3,tabm3,tabr3,tabdm3,tabmenv3,tabrenv3,
     &      tablumin3,kwv3,j,tco3,dms3,Rmax3,fb3,z,tmax3,kw,tRmax,mco0,
     &      dtp)
            tRmax3 = tRmax
            kwfin3 = kw
            scm3=scm 
            cont3=j-1
c            CALL FMMsplineDouble(tabt3,tabr3,cont3)

            R3=Rsun*10.d0**SevalDouble(0.0d0,tabt3,tabr3,cont3) 
            
            if(R3.gt.a2*(1.d0-e2)*Roche(m3,m1+m2))then 
               print*,'Initial condition unrealistic (Roche comp 3)'
               goto 2
            endif

      else ! massive outer companion
 3          continue
            call random_number(harv)
            m3 = m1+harv*(100.d0-m1)
            logPrd2 = logp(m3,0.2d0,8.d0)
            Prd2 = 10.d0**logPrd2
            e2 = ecc(m3,logPrd2)
            a2=(Prd2*sqrt(G*(M1+M2+M3))/(2.d0*pi))**(2.d0/3.d0)

            MA=2.8d0*((1.d0+M3/(M1+M2))*(1.d0+e2)/sqrt(1.d0-e2))
     &      **(2.d0/5.d0) !Mardling+Aarseth   
            if(e1.gt.1.d0.or.e2.gt.1.d0.or.a1.lt.0.d0
     &      .or.a2.lt.0.d0)then
               print*,'Initial condition unrealistic (apart)'
               goto 3
            elseif(a2*(1.d0-e2).lt.MA*a1)then
               print*,'Initial condition unrealistic (unstable)'
               goto 3
            endif
*
************************************************************************
*
*           Compute stellar evolution using sse
*       
            print*,'evolve Tertiary'
            call sse(M3,N,tabt3,tabm3,tabr3,tabdm3,tabmenv3,tabrenv3,
     &      tablumin3,kwv3,j,tco3,dms3,Rmax3,fb3,z,tmax3,kw,tRmax,mco0,
     &      dtp)
            tRmax3 = tRmax
            kwfin3 = kw
            scm3=scm 
            cont3=j-1
c            CALL FMMsplineDouble(tabt3,tabr3,cont3)

            R3=Rsun*10.d0**SevalDouble(0.0d0,tabt3,tabr3,cont3) 
            
            if(R3.gt.a2*(1.d0-e2)*Roche(m3,m1+m2))then 
               print*,'Initial condition unrealistic (Roche comp 3)'
               goto 3
            endif

      endif

!     Avoid divergences in evolve.f
      e1=max(e1,1.d-10)  
      e2=max(e2,1.d-10)  

!     Inner orbit inclination
      call random_number(harv)
      harv=2.d0*harv-1.d0
      i1=acos(harv)
      
!     Outer orbit inclination
      i2=0.d0
      
!     Ascending node of the inner orbit
      an1=0.d0
      
!     Ascending node of the outer orbit
      an2=pi-an1
      
!     Argument of the inner orbit periapse
      call random_number(harv)
      o1=harv*2.d0*pi
      
!     Argument of the outer orbit periapse
      call random_number(harv)
      o2=harv*2.d0*pi

!     Polar angle between S1 and G1
c      call random_number(harv)
c      harv=2.d0*harv-1.d0
c      Pol_S1G1=acos(harv)

!     Azimuthal angles of S1 around G1
c      call random_number(harv)
c      Azi_S1G1=harv*2.d0*pi

!     Polar angle between S2 and G1
c      call random_number(harv)
c      harv=2.d0*harv-1.d0
c      Pol_S2G1=acos(harv)

!     Azimuthal angles of S2 around G1
c      call random_number(harv)
c      Azi_S2G1=harv*2.d0*pi

!     Optional alignment of spins and G1
      Pol_S1G1=i1
      Azi_S1G1=an1
      Pol_S2G1=i1
      Azi_S2G1=an1

!     initial vectors
      y(1)=(cos(an1)*cos(o1)-sin(an1)*sin(o1)*cos(i1))*e1 !e1x
      y(2)=(sin(an1)*cos(o1)+cos(an1)*sin(o1)*cos(i1))*e1 !e1y
      y(3)=(sin(o1)*sin(i1))*e1 !e1z

      j1=sqrt(1.d0-e1**2)
      y(4)=(sin(an1)*sin(i1))*j1 !j1x
      y(5)=(-cos(an1)*sin(i1))*j1 !j1y
      y(6)=(cos(i1))*j1         !j1z
      
      y(7)=(cos(an2)*cos(o2)-sin(an2)*sin(o2)*cos(i2))*e2 !e2x
      y(8)=(sin(an2)*cos(o2)+cos(an2)*sin(o2)*cos(i2))*e2 !e2y
      y(9)=(sin(o2)*sin(i2))*e2 !e2z      

      j2=sqrt(1.d0-e2**2)
      y(10)=(sin(an2)*sin(i2))*j2 !j2x
      y(11)=(-cos(an2)*sin(i2))*j2 !j2y
      y(12)=(cos(i2))*j2        !j2z      
      
!     initial misalignment             
      y(13)=(sin(Azi_S1G1)*sin(Pol_S1G1)) !S1x/S1
      y(14)=(-cos(Azi_S1G1)*sin(Pol_S1G1)) !S1y/S1
      y(15)=(cos(Pol_S1G1)) !S1z/S1
      
      y(16)=(sin(Azi_S2G1)*sin(Pol_S2G1)) !S2x/S2
      y(17)=(-cos(Azi_S2G1)*sin(Pol_S2G1)) !S2y/S2
      y(18)=(cos(Pol_S2G1)) !S2z/S2

      Om1=45.35d0*vrotf(M1)/(R1/Rsun)*ut !from Lang 1992, and Hurley et al. 2000
      Om2=45.35d0*vrotf(M2)/(R2/Rsun)*ut   
      Sp1=Om1*K1*M1*R1**2       !spin AM = moment of inertia*stellar spin frequency
      Sp2=Om2*K2*M2*R2**2
                 
      y(13)=y(13)*Sp1 !S1x
      y(14)=y(14)*Sp1 !S1y
      y(15)=y(15)*Sp1 !S1z
      
      y(16)=y(16)*Sp2 !S2x
      y(17)=y(17)*Sp2 !S2y
      y(18)=y(18)*Sp2 !S2z
      
      y(19)=a1
      y(20)=a2

      print*,'      m1=',m1
      print*,'      m2=',m2
      print*,'      m3=',m3
      do h=1,20
         print*,'      y(',h,')=',y(h)
      enddo
      print*,'Initial masses (3), semi-major axes (2), and ecc (2):'
      print*,m1,m2,m3
      print*,a1,a2
      print*,e1,e2

*
************************************************************************
*
*     Record initial set-up
*            
      open(61,file="./stable-IC.dat",action='write',access='append')
      write(61,*)m1,m2,m3,y
      close(61)

      h=h+1
      if(h.lt.iterator)goto 1
 999  continue
      end program spinTriple


*
************************************************************************
*
*     Eq. 5
*
      real*8 function Ftwin_func(M1,logPrd)
        implicit none
        real*8 M1,logPrd,Ftwinle1,logPtwin

        Ftwinle1=0.3d0-0.15d0*LOG10(M1)
        if(M1.gt.6.5d0)then
            logPtwin=1.5d0
        else
            logPtwin=8.d0-M1
        endif
        if(logPrd.lt.1.d0)then
            Ftwin_func=Ftwinle1
        elseif(logPrd.lt.logPtwin)then
            Ftwin_func=Ftwinle1*(1.d0-(logPrd-1.d0)/(logPtwin-1.d0))
        else
            Ftwin_func=0.d0
        endif
      end function Ftwin_func
*
************************************************************************
*
*     Eq. 11, valid if M>6
*
      real*8 function gamma_lq_func(M1,logPrd)
      implicit none
      real*8 M1,logPrd

      if(M1.le.6.d0)then
        print*,'Error: primary too light. Eqs. are only valid if M>6'
      endif

      if(logPrd.lt.1.d0)then
        gamma_lq_func=-0.5d0
      elseif(logPrd.lt.2.d0)then
        gamma_lq_func=-0.5d0-0.9d0*(logPrd-1.d0)
      elseif(logPrd.lt.4.d0)then
        gamma_lq_func=-1.4d0-0.3d0*(logPrd-2.d0)
      else
        gamma_lq_func=-2.d0
      endif
      end function gamma_lq_func
*
************************************************************************
*
*     Eq. 15, valid if M>6
*
      real*8 function gamma_sq_func(M1,logPrd)
      implicit none
      real*8 M1,logPrd

      if(logPrd.lt.1.d0)then
        gamma_sq_func=0.1d0
      elseif(logPrd.lt.3.d0)then
        gamma_sq_func=0.1d0-0.15d0*(logPrd-1.d0)
      elseif(logPrd.lt.5.6d0)then
        gamma_sq_func=-0.2d0-0.5d0*(logPrd-3.d0)
      else
        gamma_sq_func=-1.5d0
      endif
      end function gamma_sq_func
*
************************************************************************
*
*     Eq. 3
*
      real*8 function emax_func(Prd)
      implicit none
      real*8 Prd

      emax_func=1.d0-(Prd/2.d0)**(-2.d0/3.d0)
      end function emax_func
*
************************************************************************
*
*     Eq. 18, valid if M>7
*
      real*8 function eta_func(M1,logPrd)
      implicit none
      real*8 M1,logPrd

      ! Note that eta is not well constrained beyond logP>5, see Sec. 9.2.
      eta_func=0.9d0-0.2d0/(logPrd-0.5d0)
      end function eta_func
*
************************************************************************
*
*     Power-law
*
      real*8 function power_law(alpha,low,high)
      implicit none
      real*8 alpha,low,high,N,harv

      call random_number(harv)
      N=(alpha+1.d0)/(high**(alpha+1.d0)-low**(alpha+1.d0)) ! Normalise
      power_law=(harv*(alpha+1.d0)/N+low**(alpha+1.d0))
     &**(1.d0/(alpha+1.d0))

      end function power_law
*
************************************************************************
*
*     Period 
*
      real*8 function logp(m1,lower_bound,upper_bound)
      implicit none
      real*8 m1
      real*8 alpha,low,high,harv,DlogP,dP
      real*8 CDF(1000000),PDF(1000000)
      real*8 flogPle1,flogPeq2p7,flogPeq5p5
      real*8 lower_bound,upper_bound
      integer h,loc
      real*8 P_iterator

      alpha=0.018d0 
      DlogP=0.7d0  
      dP=(upper_bound-lower_bound)/1000000
      P_iterator=lower_bound
      flogPle1   = 0.020d0+0.04d0*LOG10(M1)+0.07d0*LOG10(M1)**2.d0   ! Eqn. 20
      flogPeq2p7 = 0.039d0+0.07d0*LOG10(M1)+0.01d0*LOG10(M1)**2.d0   ! Eqn. 21
      flogPeq5p5 = 0.078d0-0.05d0*LOG10(M1)+0.04d0*LOG10(M1)**2.d0   ! Eqn. 22
      do h=1,1000000
        if(P_iterator.lt.1.d0)then
            PDF(h)=flogPle1
        elseif(P_iterator.lt.2.7d0-DlogP)then
            PDF(h)=flogPle1+(P_iterator-1.d0)/(1.7d0-DlogP)*(flogPeq2p7
     &             -flogPle1-alpha*DlogP)
        elseif(P_iterator.lt.2.7d0+DlogP)then
            PDF(h)=flogPeq2p7+alpha*(P_iterator-2.7d0)
        elseif(P_iterator.lt.5.5d0)then
            PDF(h)=flogPeq2p7+alpha*DlogP+(P_iterator-2.7d0-DlogP)
     &             /(2.8d0-DlogP)*(flogPeq5p5-flogPeq2p7-alpha*DlogP)
        else
            PDF(h)=flogPeq5p5*EXP(-0.3d0*(P_iterator-5.5d0))
        endif
!        if(P_iterator.lt.1.5d0)then
!            PDF(h)=PDF(h)*(1.d0)
!        else
!            PDF(h)=PDF(h)*MAX(0.d0,
!     &      1.d0-0.11d0*(logPrd2-1.5d0)**1.43d0*(M1/10.d0)**0.56d0)
!        endif
        PDF(h)=PDF(h)*dP
        P_iterator=P_iterator+dP
        if(h.eq.1)then 
            CDF(h)=PDF(h)
        else 
            CDF(h)=CDF(h-1)+PDF(h)
        endif
      end do

      CDF=CDF/SUM(PDF) ! Normalise

      call random_number(harv)

      loc=MINLOC(ABS(CDF-harv),1)
      if(harv.lt.CDF(loc))then
        logp=lower_bound+(loc-1)*dP+(harv-CDF(loc-1))
     &          /(CDF(loc)-CDF(loc-1))*dP
      else
        logp=lower_bound+loc*dP+(harv-CDF(loc))
     &          /(CDF(loc+1)-CDF(loc))*dP
      endif

      end function logp
*
************************************************************************
*
*     Mass ratio
*
      real*8 function q(m1,logPrd1)
      implicit none
      real*8 M1,logPrd1
      real*8 Ftwin,gamma_sq,gamma_lq,red,blue,white
      real*8 Ftwin_func,gamma_lq_func,gamma_sq_func
      real*8 CDF(1000000),PDF(1000000)
      real*8 lower_bound,upper_bound,inter_bound
      real*8 P_iterator,dP,q_iterator,dq
      real*8 harv
      integer h,loc

      Ftwin=Ftwin_func(M1,logPrd1)
      gamma_lq=gamma_lq_func(M1,logPrd1)
      gamma_sq=gamma_sq_func(M1,logPrd1)
      lower_bound=0.1d0
      inter_bound=0.3d0
      upper_bound=1.d0
      blue=(upper_bound**(1.d0+gamma_lq)-inter_bound**(1.d0+gamma_lq))
     &      /(1.d0+gamma_lq)
      white=(inter_bound**(1.d0+gamma_sq)-lower_bound**(1.d0+gamma_sq))
     &      /(1.d0+gamma_sq)
      red=Ftwin*blue/(1.d0-Ftwin)
      
      CDF=0.d0
      PDF=0.d0
      dq=(upper_bound-lower_bound)/1000000
      q_iterator=lower_bound
      do h=1,1000000
        if(q_iterator.lt.inter_bound)then
            PDF(h)=q_iterator**gamma_sq
        elseif(q_iterator.lt.0.95d0)then
            PDF(h)=q_iterator**gamma_lq
        else
            PDF(h)=q_iterator**gamma_lq+red/0.05d0
        endif
        PDF(h)=PDF(h)*dq
        q_iterator=q_iterator+dq
        if(h.eq.1)then 
            CDF(h)=PDF(h)
        else 
            CDF(h)=CDF(h-1)+PDF(h)
        endif
      end do

      CDF=CDF/SUM(PDF) ! Normalise

      call random_number(harv)

      loc=MINLOC(ABS(CDF-harv),1)
      if(harv.lt.CDF(loc))then
        q=lower_bound+(loc-1)*dq+(harv-CDF(loc-1))
     &          /(CDF(loc)-CDF(loc-1))*dq
      else
        q=lower_bound+loc*dq+(harv-CDF(loc))
     &          /(CDF(loc+1)-CDF(loc))*dq
      endif

      end function q
*
************************************************************************
*
*     Eccentricity
*
      real*8 function ecc(m1,logPrd1)
      implicit none
      real*8 M1,logPrd1,e1,eta,eta_func,power_law

      if(logPrd1.gt.23.d0/38.d0)then ! Period limit ensures convergence of the CDF
        eta=eta_func(M1,logPrd1)
        e1=power_law(eta,0.d0,1.d0)
      else
        e1=0.d0
      endif  
      ecc = e1
      end function ecc
