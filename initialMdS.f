********************************************************************
C     SET UP INITIAL CONDITIONS                                    C
********************************************************************

!      program initial_MdS(M1,M2,M3,e1,e2,a1,a2,i1,i2,an1,an2,o1,o2,
!     &     Prd1,Prd2,Pol_S1G1,Azi_S1G1,Pol_S2G1,Azi_S2G1,y)
      subroutine initial(M1,M2,M3,y)
      implicit none
      integer h,loc
      integer kw,nq
      real*8 pi,G
      real*8 m1,m2,m3,mp,ma
      real*8 alpha,lower_bound,upper_bound,inter_bound,Norm,harv,DlogP
      real*8 Prd1,Prd2
      real*8 logPrd1,logPrd2
      real*8 a1,a2
      real*8 flogPle1,flogPeq2p7,flogPeq5p5
      real*8 CDF(1000000),PDF(1000000)
      real*8 P_iterator,dP,q_iterator,dq
      real*8 Ftwin,gamma_sq,gamma_lq,red,blue,white
      real*8 Ftwin_func,gamma_lq_func,gamma_sq_func
      real*8 q_in,q_out
      real*8 emax_func,eta,eta_func,e1,e2
      real*8 P_hier
      real*8 j1,j2
      real*8 i1,i2
      real*8 an1,an2
      real*8 o1,o2
      real*8 Pol_S1G1,Azi_S1G1,Pol_S2G1,Azi_S2G1
      real*8 power_law
      real*8 low_mass
      parameter (nq=20)
      real*8 y(nq)
*
************************************************************************
*
*   
*     
      low_mass=8.d0
      pi=2.d0*ASIN(1.d0)
      G=2.959d-4 ! in AU^3/Solar mass/day^2
*
************************************************************************
*
*     Primary mass
*
 1    continue
      M1=power_law(-2.3d0,low_mass,1000000.d0)
*
************************************************************************
*
*     Inner Period
*
      alpha=0.018d0 
      DlogP=0.7d0  
      lower_bound=0.2d0
      upper_bound=5.5d0
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
        logPrd1=lower_bound+(loc-1)*dP+(harv-CDF(loc-1))
     &          /(CDF(loc)-CDF(loc-1))*dP
        Prd1=10.d0**logPrd1
      else
        logPrd1=lower_bound+loc*dP+(harv-CDF(loc))
     &          /(CDF(loc+1)-CDF(loc))*dP
        Prd1=10.d0**logPrd1
      endif
*
************************************************************************
*
*     Inner mass ratio
*
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
        q_in=lower_bound+(loc-1)*dq+(harv-CDF(loc-1))
     &          /(CDF(loc)-CDF(loc-1))*dq
      else
        q_in=lower_bound+loc*dq+(harv-CDF(loc))
     &          /(CDF(loc+1)-CDF(loc))*dq
      endif
      M2=q_in*M1
*
************************************************************************
*
*     Ensure massive inner binary
*
      if(M1.lt.low_mass.or.M2.lt.low_mass)then
            goto 1
      else
            a1=(Prd1*sqrt(G*(M1+M2))/(2.d0*pi))**(2.d0/3.d0) ! in AU
      endif
*
************************************************************************
*
*     Inner eccentricity
*
      if(logPrd1.gt.23.d0/38.d0)then ! Period limit ensures convergence of the CDF
        eta=eta_func(M1,logPrd1)
        e1=power_law(eta,0.d0,1.d0)
      else
        e1=0.d0
      endif  
*
************************************************************************
*
*     Outer mass ratio, agnostic
*
      lower_bound=0.1d0
      upper_bound=1.d0
      q_out=power_law(0.0d0,lower_bound,upper_bound)
      M3=q_out*(M1+M2)
      MP=MAX(M3,M1+M2) ! Determine primary mass
*
************************************************************************
*
*     Outer Period
*
      alpha=0.018d0 
      DlogP=0.7d0  
      lower_bound=0.2d0
      upper_bound=8.d0
      dP=(upper_bound-lower_bound)/1000000
      P_iterator=lower_bound
      flogPle1   = 0.020d0+0.04d0*LOG10(MP)+0.07d0*LOG10(MP)**2.d0   ! Eqn. 20
      flogPeq2p7 = 0.039d0+0.07d0*LOG10(MP)+0.01d0*LOG10(MP)**2.d0   ! Eqn. 21
      flogPeq5p5 = 0.078d0-0.05d0*LOG10(MP)+0.04d0*LOG10(MP)**2.d0   ! Eqn. 22
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
!            PDF(h)=PDF(h)*(1.d0-1.d0)
!        else
!            PDF(h)=PDF(h)*(1.d0-MAX(0.d0,
!     &      1.d0-0.11d0*(logPrd2-1.5d0)**1.43d0*(MP/10.d0)**0.56d0))
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
        logPrd2=lower_bound+(loc-1)*dP+(harv-CDF(loc-1))
     &          /(CDF(loc)-CDF(loc-1))*dP
        Prd2=10.d0**logPrd2
      else
        logPrd2=lower_bound+loc*dP+(harv-CDF(loc))
     &          /(CDF(loc+1)-CDF(loc))*dP
        Prd2=10.d0**logPrd2
      endif
      a2=(Prd2*sqrt(G*(M1+M2+M3))/(2.d0*pi))**(2.d0/3.d0) ! in AU
*
************************************************************************
*
*     Outer eccentricity
*
      if(logPrd2.gt.23.d0/38.d0)then ! Period limit ensures convergence of the CDF
        eta=eta_func(MP,logPrd2)
        e2=power_law(eta,0.d0,1.d0)
      else
        e2=0.d0
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
      call random_number(harv)
      harv=2.d0*harv-1.d0
      Pol_S1G1=acos(harv)

!     Azimuthal angles of S1 around G1
      call random_number(harv)
      Azi_S1G1=harv*2.d0*pi

!     Polar angle between S2 and G1
      call random_number(harv)
      harv=2.d0*harv-1.d0
      Pol_S2G1=acos(harv)

!     Azimuthal angles of S2 around G1
      call random_number(harv)
      Azi_S2G1=harv*2.d0*pi

!     Optional alignment of spins and G1
      Pol_S1G1=i1
      Azi_S1G1=an1
      Pol_S2G1=i1
      Azi_S2G1=an1

!     initial vectors
      y(1)=(cos(an1)*cos(o1)-sin(an1)*sin(o1)*cos(i1))*e1 !e1x
      y(2)=(sin(an1)*cos(o1)+cos(an1)*sin(o1)*cos(i1))*e1 !e1y
      y(3)=(sin(o1)*sin(i1))*e1 !e1z

      j1=sqrt(1.-e1**2)
      y(4)=(sin(an1)*sin(i1))*j1 !j1x
      y(5)=(-cos(an1)*sin(i1))*j1 !j1y
      y(6)=(cos(i1))*j1         !j1z
      
      y(7)=(cos(an2)*cos(o2)-sin(an2)*sin(o2)*cos(i2))*e2 !e2x
      y(8)=(sin(an2)*cos(o2)+cos(an2)*sin(o2)*cos(i2))*e2 !e2y
      y(9)=(sin(o2)*sin(i2))*e2 !e2z      

      j2=sqrt(1.-e2**2)
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
      
      y(19)=a1
      y(20)=a2

 98   continue
      RETURN
      end
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