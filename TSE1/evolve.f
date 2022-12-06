***********************************************************************
*     CALCULATE TIME DERIVATIVES                                      *
***********************************************************************
      subroutine evolution(t,y,yp)
      use my_subs
      use commonV
      use tables
      implicit real*8 (a-z)
      integer n,k,h,i(1)
      parameter  (n=20)
      real*8 y(n),yp(n)
      real*8 j(3),edot(3),jdot(3),e2dot(3),j2dot(3)
      real*8 S1dot(3),S2dot(3),S1(3),S2(3),S0(3)
      real*8 ev(3),jv(3),ev2(3),jv2(3)
      real*8 SevalDouble
      real*8 cross_jv_ev(3)
      real*8 cross_jv_jv2_j2t(3)
      real*8 cross_ev_jv2_j2t(3)
      real*8 cross_jv_ev2_e2(3)
      real*8 cross_ev_ev2_e2(3)
      real*8 cross_ev2_jv(3)
      real*8 cross_ev2_ev(3)
      real*8 cross_jv2_j2t_ev2(3)
      real*8 cross_jv2_j2t_jv(3)
      real*8 cross_jv2_j2t_ev(3)
      real*8 cross_jv_j1_ev(3)
      real*8 cross_ev2_jv2_j2t(3)
      real*8 cross_ev2_e2_ev(3)
      real*8 cross_ev2_e2_jv(3)
      real*8 cross_jv_j1_S1(3)
      real*8 cross_jv_j1_S2(3)
      real*8 jv_j1(3)
      real*8 ev_e1(3)
      real*8 jv2_j2(3)
      real*8 ev2_e2(3)

      c2 = c**2
      c5 = c**5
      
      if(t.lt.tco1)then
            m1=SevalDouble(t,tabt1,tabm1,cont1)
            dm1=-10.d0**SevalDouble(t,tabt1,tabdm1,cont1) 
      else
            m1=mco1
            dm1=0.d0
      endif
      if(t.lt.tco2)then
            m2=SevalDouble(t,tabt2,tabm2,cont2)
            dm2=-10.d0**SevalDouble(t,tabt2,tabdm2,cont2)
      else
            m2=mco2
            dm2=0.d0
      endif 
      if(t.lt.tco3)then
            m3=SevalDouble(t,tabt3,tabm3,cont3)
            dm3=-10.d0**SevalDouble(t,tabt3,tabdm3,cont3)
      else
            m3=mco3
            dm3=0.d0
      endif 
      m12 = m1+m2
      m123 = m12+m3
      reduc1=m1*m2/(m1+m2)
      reduc2=(m1+m2)*m3/(m1+m2+m3)
      R1=Rsun*10.d0**SevalDouble(t,tabt1,tabr1,cont1)
      if(R1.gt.Rmax1)R1=Rmax1
      R2=Rsun*10.d0**SevalDouble(t,tabt2,tabr2,cont2)
      if(R2.gt.Rmax2)R2=Rmax2

!      R3=Rsun*10.d0**SevalDouble(t,tabt3,tabr3,cont3)
!      if(R3.gt.Rmax3)R3=Rmax3   !avoid overpredicting R 

      yp=0.d0 
      do h=1,3
         ev(h)=y(h)    
         jv(h)=y(h+3)
         ev2(h)=y(h+6)     
         jv2(h)=y(h+9)
         S1(h)=y(12+h)
         S2(h)=y(15+h)
         S0(h)=(1.d0+m2/m1)*S1(h)+(1.d0+m1/m2)*S2(h)
      end do
**      
!     define various quantities      
      a1=y(19)
      a2=y(20)
      e1=sqrt(y(1)**2+y(2)**2+y(3)**2)  
      e2=sqrt(y(7)**2+y(8)**2+y(9)**2)
      j1=sqrt(y(4)**2+y(5)**2+y(6)**2)
      j2t=sqrt(y(10)**2+y(11)**2+y(12)**2)


c     Avoid divergences 
      e1=max(e1,1.d-10)
      e2=max(e2,1.d-10)
      j1=min(j1,sqrt(1.-(1.d-10)**2)) 
      j2t=min(j2t,sqrt(1.-(1.d-10)**2)) 

      do h=1,3
            y(h)=y(h)*e1/sqrt(dotp(ev,ev))
            y(h+3)=y(h+3)*j1/sqrt(dotp(jv,jv))
            y(h+6)=y(h+6)*e2/sqrt(dotp(ev2,ev2))
            y(h+9)=y(h+9)*j2t/sqrt(dotp(jv2,jv2))
      end do
      do h=1,3
            ev(h)=y(h)    
            jv(h)=y(h+3)
            ev2(h)=y(h+6)     
            jv2(h)=y(h+9)
      enddo
      
c      if(e1.gt.1.d0.or.e1.lt.0.d0)print*,'e1:',e1,t*(58.d0/365.d0)/1.d6
      
      mm=sqrt(m12/a1**3)
      
      L1=reduc1*sqrt(m12*a1)
      L2=reduc2*sqrt(m123*a2)

      e1p2 = e1**2
      e1p4 = e1**4
      e1p6 = e1**6
     
!     eccentricity dependent coefficients
      f1=1.d0+31.d0*e1p2/2.d0+255.d0*e1p4/8.d0+
     &   185.d0*e1p6/16.d0+25.d0*e1**8/64.d0
      f2=1.d0+15.d0*e1p2/2.d0+45.d0*e1p4/8.d0+5.d0*e1p6/16.d0
      f3=1.d0+15.d0*e1p2/4.d0+15.d0*e1p4/8.d0+5.d0*e1p6/64.d0
      f4=1.d0+3.d0*e1p2/2.d0+e1p4/8.d0
      f5=1.d0+3.d0*e1p2+3.d0*e1p4/8.d0

!     spin modules
      S1t=sqrt(S1(1)**2+S1(2)**2+S1(3)**2)
      S2t=sqrt(S2(1)**2+S2(2)**2+S2(3)**2)
      
!     Some dot-prod

      jv_j1 = jv/j1
      ev_e1 = ev/e1
      jv2_j2 = jv2/j2t
      ev2_e2 = ev2/e2

      e1j2=dotp(ev,jv2_j2)
      j1j2=dotp(jv,jv2_j2)
      e1e2=dotp(ev,ev2_e2)
      j1e2=dotp(jv,ev2_e2)


      cross_jv_j1_ev=cross(jv_j1,ev)
      cross_jv_j1_S1=cross(jv_j1,S1)
      cross_jv_j1_S2=cross(jv_j1,S2)
*
************************************************************************
*
*     Lidov-Kozai
*

!      if(m3.lt.1d-3)goto 100
      TLK=m12*(a2*j2t/a1)**3/mm/m3*4.d0/3.d0
      Om_LK = 1.d0/TLK 
      oct=(m1-m2)/m12*a1/a2*e2/j2t**2*75.d0/64.d0*(4.d0/3.d0)

      cross_jv_ev=cross(jv,ev)
      cross_jv_jv2_j2t=cross(jv,jv2_j2)
      cross_ev_jv2_j2t=cross(ev,jv2_j2)
      cross_jv_ev2_e2=cross(jv,ev2_e2)
      cross_ev_ev2_e2=cross(ev,ev2_e2)
      cross_ev2_jv=cross(ev2,jv)
      cross_ev2_ev=cross(ev2,ev)
      cross_jv2_j2t_ev2=cross(jv2_j2,ev2)
      cross_jv2_j2t_jv=cross(jv2_j2,jv)
      cross_jv2_j2t_ev=cross(jv2_j2,ev)

      cross_ev2_jv2_j2t=cross(ev2,jv2_j2)
      cross_ev2_e2_ev=cross(ev2_e2,ev)
      cross_ev2_e2_jv=cross(ev2_e2,jv)
           
      edot=(2.d0*cross_jv_ev
     &     -5.d0*e1j2*cross_jv_jv2_j2t
     &     +j1j2*cross_ev_jv2_j2t)*Om_LK 
     &     -(2.d0*e1j2*j1j2*cross_ev_ev2_e2
     &     +(8.d0/5.d0*e1p2-1.d0/5.d0-7.d0*e1j2**2+j1j2**2)
     &     *cross_jv_ev2_e2
     &     +2.d0*(e1e2*j1j2+e1j2*j1e2)*cross_ev_jv2_j2t
     &     +2.d0*(j1j2*j1e2-7.d0*e1j2*e1e2)*cross_jv_jv2_j2t
     &     +16.d0/5.d0*e1e2*cross_jv_ev)
     &     *oct*Om_LK
      
      jdot=-(2.d0*(e1e2*j1j2+e1j2*j1e2)*cross_jv_jv2_j2t
     &     +2.d0*(j1e2*j1j2-7.d0*e1e2*e1j2)*cross_ev_jv2_j2t
     &     +2.d0*(e1j2*j1j2)*cross_jv_ev2_e2
     &     +(8.d0/5.d0*e1p2-1.d0/5.d0-7.d0*e1j2**2
     &     +j1j2**2)*cross_ev_ev2_e2)
     &     *oct*Om_LK
     &     +(j1j2*cross_jv_jv2_j2t
     &     -5.d0*e1j2*cross_ev_jv2_j2t)*Om_LK !djx/dt

      yp(1)=yp(1)+edot(1)*LK
      yp(2)=yp(2)+edot(2)*LK
      yp(3)=yp(3)+edot(3)*LK       
      yp(4)=yp(4)+jdot(1)*LK  
      yp(5)=yp(5)+jdot(2)*LK
      yp(6)=yp(6)+jdot(3)*LK
      
      e2dot=(j1j2*cross_ev2_jv
     &     -5.d0*e1j2*cross_ev2_ev
     &     -(1.d0/2.d0-3.d0*e1p2+25.d0/2.d0*e1j2**2
     &     -5.d0/2.d0*j1j2**2)*cross_jv2_j2t_ev2)
     &     *Om_LK*L1/L2/j2t
     &     -(2.d0*(e1j2*j1e2*cross_ev2_jv
     &     +j1j2*e1e2*cross_ev2_jv
     &     +j2t**2/e2*e1j2*j1j2*cross_jv2_j2t_jv)
     &     +2.d0*j1e2*j1j2*cross_ev2_ev
     &     -14.d0*e1e2*e1j2*cross_ev2_ev
     &     +j2t**2/e2*(8.d0/5.d0*e1p2-1.d0/5.d0
     &     -7.d0*e1j2**2+j1j2**2)*cross_jv2_j2t_ev
     &     -2.d0*(1.d0/5.d0-8.d0/5.d0*e1p2)*e1e2*cross_ev2_jv2_j2t
     &     -14.d0*e1j2*j1e2*j1j2*cross_ev2_jv2_j2t
     &     -7.d0*e1e2*(8.d0/5.d0*e1p2-1.d0/5.d0-7.d0*e1j2**2
     &     +j1j2**2)*cross_ev2_jv2_j2t)
     &     *oct*Om_LK*L1/L2/j2t

      j2dot=(j1j2*cross_jv2_j2t_jv
     &     -5.d0*e1j2*cross_jv2_j2t_ev)*Om_LK*L1/L2
     &     -(2.d0*(e1j2*j1e2*cross_jv2_j2t_jv
     &     +e1e2*j1j2*cross_jv2_j2t_jv
     &     +e1j2*j1j2*cross_ev2_e2_jv)
     &     +2.d0*j1e2*j1j2*cross_jv2_j2t_ev
     &     -14.d0*e1e2*e1j2*cross_jv2_j2t_ev
     &     +(8.d0/5.d0*e1p2-1.d0/5.d0-7.d0*e1j2**2+j1j2**2)
     &     *cross_ev2_e2_ev)
     &     *oct*Om_LK*L1/L2

      yp(7)=yp(7)+e2dot(1)*LK
      yp(8)=yp(8)+e2dot(2)*LK
      yp(9)=yp(9)+e2dot(3)*LK
      yp(10)=yp(10)+j2dot(1)*LK
      yp(11)=yp(11)+j2dot(2)*LK
      yp(12)=yp(12)+j2dot(3)*LK


*
************************************************************************
*
*     Schwarzschild precession
*
 100  continue
      om_gr=3.d0*mm*m12/c2/(a1*j1**2)

      edot=om_gr*cross_jv_j1_ev
      
      yp(1)=yp(1)+edot(1)*GR
      yp(2)=yp(2)+edot(2)*GR
      yp(3)=yp(3)+edot(3)*GR   
*
************************************************************************
*
*     Gravitational wave emission
*      
      e_GW=-304.d0*m1*m2*m12*e1*(1.d0+e1p2*121.d0/304.d0)
     &     /(15.d0*c5*a1**4*((1.d0-e1p2)**(5.d0/2.d0)))
      l_GW=-e_GW*e1/sqrt(1.d0-e1p2)
      a_GW=(-64.d0*m1*m2*m12/(5.d0*(c5)*(a1**3)
     &     *(1.d0-e1p2)**3.5d0))
     &     *(1.d0+(e1p2)*73.d0/24.d0+(e1p4)*37.d0/96.d0)
      
      edot=e_GW*ev_e1        
      jdot=l_GW*jv_j1 

      yp(1)=yp(1)+edot(1)*GW
      yp(2)=yp(2)+edot(2)*GW
      yp(3)=yp(3)+edot(3)*GW       
      yp(4)=yp(4)+jdot(1)*GW      
      yp(5)=yp(5)+jdot(2)*GW
      yp(6)=yp(6)+jdot(3)*GW
      yp(19)=yp(19)+a_GW*GW
*
************************************************************************
*
*     de Sitter precession
* 
      Om_SO_1=2.d0*reduc1*mm*(1.d0+3.d0*m2/(4.d0*m1))/(c2*a1*j1**2)
      Om_SO_2=2.d0*reduc1*mm*(1.d0+3.d0*m1/(4.d0*m2))/(c2*a1*j1**2)
      
      S1dot=Om_SO_1*cross_jv_j1_S1
      S2dot=Om_SO_2*cross_jv_j1_S2

      yp(13)=yp(13)+S1dot(1)*Sitter
      yp(14)=yp(14)+S1dot(2)*Sitter
      yp(15)=yp(15)+S1dot(3)*Sitter     
      yp(16)=yp(16)+S2dot(1)*Sitter
      yp(17)=yp(17)+S2dot(2)*Sitter
      yp(18)=yp(18)+S2dot(3)*Sitter
*
************************************************************************
*
*     Lense-Thirring precession
* 
      edot=2.d0*(1.d0+3.d0*m2/(4.d0*m1))/(c2*a1**3*j1**3)
     &      *cross(S1-3.d0*dotp(S1,jv_j1),ev)
     &    +2.d0*(1.d0+3.d0*m1/(4.d0*m2))/(c2*a1**3*j1**3)
     &      *cross(S2-3.d0*dotp(S2,jv_j1),ev)

      jdot=2.d0*(1.d0+3.d0*m2/(4.d0*m1))/(c2*a1**3*j1**3)
     &      *cross(S1,jv)
     &    +2.d0*(1.d0+3.d0*m1/(4.d0*m2))/(c2*a1**3*j1**3)
     &      *cross(S2,jv)

      yp(1)=yp(1)+edot(1)*Thirring      
      yp(2)=yp(2)+edot(2)*Thirring       
      yp(3)=yp(3)+edot(3)*Thirring     
      yp(4)=yp(4)+jdot(1)*Thirring       
      yp(5)=yp(5)+jdot(2)*Thirring 
      yp(6)=yp(6)+jdot(3)*Thirring 

!     EVOLUTION OF SPINS due to spin-spin interaction, Rodriguez & Antonini (2018): Eq. 20
!      Om_SS_1=reduc1*(1.d0+m2/m1)/(2.d0*c2*m12*a1**3*j1**3)
!      Om_SS_2=reduc1*(1.d0+m1/m2)/(2.d0*c2*m12*a1**3*j1**3)

!      S1dot=Om_SS_1*cross(S0-3.d0*dotp(S0,jv_j1)*jv_j1,S1) 
!      S2dot=Om_SS_2*cross(S0-3.d0*dotp(S0,jv_j1)*jv_j1,S2) 

!      yp(13)=yp(13)+S1dot(1)
!      yp(14)=yp(14)+S1dot(2)
!      yp(15)=yp(15)+S1dot(3)
!      yp(16)=yp(16)+S2dot(1)
!      yp(17)=yp(17)+S2dot(2)
!      yp(18)=yp(18)+S2dot(3)    
*
************************************************************************
*
*     If primary is CO skip its remaining stellar physics
*
      if(js1.eq.1) goto 11
      !dm1=-abs(SevalDouble(t,tabt1,tabdm1,cont1))
!      dm1=-10.d0**SevalDouble(t,tabt1,tabdm1,cont1) 
*
************************************************************************
*
*     Tides, Correia et al. (2011)
*     star 1
*

*      Dt1=1.996d0/1.d7 !lag time (1 sec); 1sec/58days->1.996d0/1.d7

!      R1=Rsun*10.d0**SevalDouble(t,tabt1,tabr1,cont1)
      if(R1.gt.Rmax1)R1=Rmax1   !avoid overpredicting R
      Om1=S1t/(k1*m1*R1**2)     !rotation rate
      Om1_n=Om1/sqrt(m1/R1**3)  !rotation rate in units of breakup freq  
c      i=MINLOC(ABS(tabt1-t))
c      kw=kwv1(i(1))
c      if(((kw.eq.1).and.m1.ge.1.25d0).or.kw.eq.4.or.kw.eq.7)then
c            TidE2 = 1.592d-09*(m1**2.84d0)
c            kT = 1.9782d+04*m1*(R1/Rsun)**2/((a1/Rsun)**5)*(1.d0+m2/m1)
c     &      **(5.d0/6.d0)*TidE2*(58.d0/365.d0)
c      elseif(((kw.eq.1).and.m1.lt.1.25d0).or.kw.eq.0.or.kw.eq.2.or.
c     &      kw.eq.3.or.kw.eq.5.or.kw.eq.6.or.kw.eq.8.or.kw.eq.9)then
c            menv1 = SevalDouble(t,tabt1,tabmenv1,p10,p11,p12,cont1)
c            renv1 = SevalDouble(t,tabt1,tabrenv1,p13,p14,p15,cont1)
c            menv1 = MAX(menv1,1d-10)
c            renv1 = MAX(renv1,1d-10)
c            lum1 = 10.d0**SevalDouble(t,tabt1,tablumin1,p16,p17,p18,
c     &            cont1)
c            tconv = 0.4311d0*(menv1*renv1*MAX(R1/Rsun-renv1/2.d0,1d-10)
c     &            /3.d0/lum1)**(1.d0/3.d0)*(365.d0/58.d0)
c            Ptid = 1.d0/(1.0d-10+ABS(mm-Om1)) 
c            fconv = MIN(1.d0,(Ptid/2.d0/tconv)**2)
c            kT = 2.d0*fconv*menv1/21.d0/tconv/m1
c            Dt1 = R1**3/m1/(kL1/2.d0)*kT
c            print*,'T1',t*58.d0/365.d0/1.d6,Dt1*58.d0*86400.d0,m1,
c     &tabm1(i(1)),
c     &R1/Rsun,10.d0**tabr1(i(1)),renv1/2.d0,tabrenv1(i(1))/2.d0
c      elseif(kw.ge.10.and.kw.le.13)then
c            lum1 = 10.d0**SevalDouble(t,tabt1,tablumin1,p16,p17,p18,
c     &            cont1)
c            kT = 2.564d-08*0.21d0*(lum1/m1)**(5.d0/7.d0)
c     &            *(58.d0/365.d0)
c      else
c            goto 11
c      endif
c      TidE2 = 1.592d-09*(m1**2.84d0)
c      kT = 1.9782d+04*m1*R1**2/(a1**5)*(1.d0+m2/m1)**(5.d0/6.d0)
c     &      *TidE2*(58.d0/365.d0)
c      Dt1 = R1**3/m1/(kL1/2.d0)*kT

      Dt1 = 1.996d-7 ! 1 Sec lag time

*      ta1=1.d0/(6.d0*KL1*Dt1*m2/m1*(R1/a1)**5*mm**2)

      adot_td=2.d0*(3.d0*kL1*Dt1*m2**2*R1**5)/(a1**7)/reduc1
     &      *(f2/j1**12*Om1/mm*dotp(S1/S1t,jv_j1)-f1/j1**15)

*     No Orbit-average
      edot=7.5d0*kL1*mm*m2/m1*(R1/a1)**5*f4/j1**10
     &      *cross_jv_j1_ev
     &      -(3.d0*kL1*Dt1*m2**2*R1**5)/(a1**8*reduc1)
     &      *(
     &      Om1/(2.d0*mm)*dotp(S1/S1t,ev)*f4/j1**10*jv_j1
     &      +(9.d0*f3/j1**13-11.d0/2.d0*Om1/mm*dotp(S1/S1t,jv_j1)
     &      *f4/j1**10)*ev
     &      )

*     Orbit-averaged
*      edot=-1.d0/(2.d0*ta1*j1**13)*
*     &     (j1**3*f4*Om1/(2.d0*mm)*dotp(ev,S1/S1t)*jv_j1
*     &     -(11.d0/2.d0*j1**3*f4*Om1/mm-9.d0*f3)*ev)
*     &     -1.d0/(2.d0*ta2*j1**13)*
*     &     (j1**3*f4*Om2/(2.d0*mm)*dotp(ev,S2/S2t)*jv_j1
*     &     -(11.d0/2.d0*j1**3*f4*Om2/mm-9.d0*f3)*ev)

*     No Orbit-average
*      S1dot=(3.d0*kL1*Dt1*m2**2*R1**5)/(a1**6*j1**9)
*     &      *(
*     &      Om1/(2.d0*mm)*dotp(S1/S1t,cross(jv_j1,ev_e1))*(f4-2.d0*f5)
*     &      *cross(jv_j1,ev_e1)
*     &      +(f2/j1**3-Om1/mm*dotp(S1/S1t,jv_j1)*f5)*jv_j1
*     &      -Om1/(2.d0*mm)*dotp(S1/S1t,ev_e1)*f4*ev_e1
*     &      )  

*     No Orbit-average
      S1dot=(3.d0*kL1*Dt1*m2**2*R1**5)/(a1**6)*mm
     &      *(
     &      f4/j1**9*Om1/(2.d0*mm)*(S1/S1t-dotp(S1/S1t,jv_j1)*jv_j1)
     &      -f5/j1**9*Om1/mm*S1/S1t
     &      +f2/j1**12*jv_j1
     &      +dotp(ev,S1/S1t)*(6.d0+e1p2)/(4.d0*j1**9)*Om1/mm*ev
     &      )

*     Orbit-averaged
*      S1dot=-(j1*L1)/(2.d0*ta1*j1**13)
*     &     *(j1**3*f5*Om1/(2.d0*mm)*(S1/S1t+dotp(S1/S1t,jv_j1)*jv_j1)
*     &     -f2*jv_j1)
*      S2dot=-(j1*L1)/(2.d0*ta2*j1**13)
*     &     *(j1**3*f5*Om2/(2.d0*mm)*(S2/S2t+dotp(S2/S2t,jv_j1)*jv_j1)
*     &     -f2*jv_j1)

      jdot=-S1dot/L1-adot_td/(2.d0*a1)*jv

      yp(1)=yp(1)+edot(1)*Tides      
      yp(2)=yp(2)+edot(2)*Tides         
      yp(3)=yp(3)+edot(3)*Tides       
      yp(4)=yp(4)+jdot(1)*Tides         
      yp(5)=yp(5)+jdot(2)*Tides   
      yp(6)=yp(6)+jdot(3)*Tides   
      yp(13)=yp(13)+S1dot(1)*Tides   
      yp(14)=yp(14)+S1dot(2)*Tides   
      yp(15)=yp(15)+S1dot(3)*Tides         
      yp(19)=yp(19)+adot_td*Tides    
*
************************************************************************
*
*     Rotational distortion star 1
*

      om_r_1=kL1/2.d0*mm/(j1**4*a1**2)*Om1_n**2*R1**2

      Ops_1=-kL1*m2*R1**3*Om1/(2.d0*k1*m1*a1**3*j1**3)
     &      *dotp(jv_j1,S1/S1t)

      S1dot=Ops_1*cross_jv_j1_S1 !spin-orbit precession ok

      jdot=-S1dot/L1 ! ok
      
      yp(13)=yp(13)+S1dot(1)*Rot
      yp(14)=yp(14)+S1dot(2)*Rot
      yp(15)=yp(15)+S1dot(3)*Rot  

      edot=(-om_r_1*(dotp(S1/S1t,jv_j1)*cross(S1/S1t,ev)
     &     +1.d0/2.d0*(1.d0-5.d0*dotp(S1/S1t,jv_j1)**2)
     &     *cross_jv_j1_ev))

      yp(1)=yp(1)+edot(1)*Rot      
      yp(2)=yp(2)+edot(2)*Rot      
      yp(3)=yp(3)+edot(3)*Rot    
      yp(4)=yp(4)+jdot(1)*Rot      
      yp(5)=yp(5)+jdot(2)*Rot
      yp(6)=yp(6)+jdot(3)*Rot    
*
************************************************************************
*
*     Mass-loss star 1 
* 
      S1dot=(K1*dm1*R1**2*Om1)*(S1/S1t)   

      yp(13)=yp(13)+S1dot(1)*Mass
      yp(14)=yp(14)+S1dot(2)*Mass
      yp(15)=yp(15)+S1dot(3)*Mass      
      yp(19)=yp(19)-dm1/m12*y(19)*Mass
      yp(20)=yp(20)-dm1/m123*y(20)*Mass 
*
************************************************************************
*
*     If secondary is CO skip its remaining stellar physics
*
 11   continue
      if(js2.eq.1) goto 22
      !dm2=-abs(SevalDouble(t,tabt2,tabdm2,cont2))
!      dm2=-10.d0**SevalDouble(t,tabt2,tabdm2,cont2)  
*
************************************************************************
*
*     Tides, Correia et al. (2011)
*     star 2
*
*      Dt2=1.996d0/1.d7

!      R2=Rsun*10.d0**SevalDouble(t,tabt2,tabr2,cont2)
      if(R2.gt.Rmax2)R2=Rmax2   !avoid overpredicting R
      Om2=S2t/(k2*m2*R2**2)
      Om2_n=Om2/sqrt(m2/R2**3)  
c      menv2 = Rsun*LinearInter(t,tabt2,tabmenv2,cont2)
c      renv2 = Rsun*LinearInter(t,tabt2,tabrenv2,cont2)
c      i=MINLOC(ABS(tabt2-t))
c      kw=kwv2(i(1))
c      if(kw.ge.10.and.kw.le.13)then
c            lum2 = 10.d0**SevalDouble(t,tabt2,tablumin2,d16,d17,d18,
c     &            cont2)
c            kT = 2.564d-08*0.21d0*(lum2/m2)**(5.d0/7.d0)
c     &            *(58.d0/365.d0)
c      elseif(((kw.eq.1).and.m2.ge.1.25d0).or.kw.eq.4.or.kw.eq.7.or.
c     &renv2.lt.1d-10.or.menv2.lt.1d-10)then
c            TidE2 = 1.592d-09*(m2**2.84d0)
c            kT = 1.9782d+04*m2*(R2/Rsun)**2/((a1/Rsun)**5)*(1.d0+m1/m2)
c     &      **(5.d0/6.d0)*TidE2*(58.d0/365.d0)
c            Dt2 = R2**3/m2/(kL2/2.d0)*kT
c      print*,'T2Dyn',t*58.d0/365.d0/1.d6,Dt2*58.d0*86400.d0
c      elseif(((kw.eq.1).and.m2.lt.1.25d0).or.kw.eq.0.or.kw.eq.2.or.
c     &      kw.eq.3.or.kw.eq.5.or.kw.eq.6.or.kw.eq.8.or.kw.eq.9)then
c            menv2 = SevalDouble(t,tabt2,tabmenv2,d10,d11,d12,cont2)
c            renv2 = SevalDouble(t,tabt2,tabrenv2,d13,d14,d15,cont2)
c            R2=Rsun*10.d0**LinearInter(t,tabt2,tabr2,cont2)
c            if(R2.gt.Rmax2)R2=Rmax2 
c            menv2 = MAX(menv2,1d-10) ! avoid underpredicting
c            renv2 = MAX(renv2,1d-10) ! avoid underpredicting
c            menv2 = MIN(menv2,m2) ! avoid overpredicting
c            renv2 = MIN(renv2,2.d0*R2-1.d-10) ! avoid overpredicting
c            lum2 = 10.d0**SevalDouble(t,tabt2,tablumin2,d16,d17,d18,
c     &            cont2)
c            tconv = 67.4d0*(menv2*renv2*(R2-renv2/2.d0)/lum2)
c     &            **(1.d0/3.d0)
c            Ptid = 1.d0/(1.0d-10+ABS(mm-Om2)) 
c            fconv = MIN(1.d0,(Ptid/2.d0/tconv)**2)
c            kT = 2.d0*fconv*menv2/21.d0/tconv/m2
c            Dt2 = R2**3/m2/(kL2/2.d0)*kT
c            Dt2=1.996d0/1.d7
c            print*,'T2',t*58.d0/365.d0/1.d6,Dt2*58.d0*86400.d0,
c     &1.d0/kT*58.d0/365.d0/1.d6,m2,menv2,
c     &R2,renv2,lum2
c      else
c            goto 22
c      endif
c      TidE2 = 1.592d-09*(m2**2.84d0)
c      kT = 1.9782d+04*m2*R2**2/(a1**5)*(1.d0+m1/m2)**(5.d0/6.d0)
c     &      *TidE2*(58.d0/365.d0)
c      Dt2 = R2**3/m2/(kL2/2.d0)*kT

      Dt2 = 1.996d-7 ! 1 Sec lag time
*      ta2=1.d0/(6.d0*KL2*Dt2*m2/m1*(R2/a1)**5*mm**2)

      adot_td=2.d0*(3.d0*kL2*Dt2*m1**2*R2**5)/(a1**7)/reduc1
     &      *(f2/j1**12*Om2/mm*dotp(S2/S2t,jv_j1)-f1/j1**15)

*     No Orbit-average
      edot=7.5d0*kL2*mm*m1/m2*(R2/a1)**5*f4/j1**10
     &      *cross_jv_j1_ev
     &      -(3.d0*kL2*Dt2*m1**2*R2**5)/(a1**8*reduc1)
     &      *(
     &      Om2/(2.d0*mm)*dotp(S2/S2t,ev)*f4/j1**10*jv_j1
     &      +(9.d0*f3/j1**13-11.d0/2.d0*Om2/mm*dotp(S2/S2t,jv_j1)
     &      *f4/j1**10)*ev
     &      )

*     Orbit-averaged
*      edot=-1.d0/(2.d0*ta1*j1**13)*
*     &     (j1**3*f4*Om1/(2.d0*mm)*dotp(ev,S1/S1t)*jv_j1
*     &     -(11.d0/2.d0*j1**3*f4*Om1/mm-9.d0*f3)*ev)
*     &     -1.d0/(2.d0*ta2*j1**13)*
*     &     (j1**3*f4*Om2/(2.d0*mm)*dotp(ev,S2/S2t)*jv_j1
*     &     -(11.d0/2.d0*j1**3*f4*Om2/mm-9.d0*f3)*ev)

*     No Orbit-average
*      S2dot=(3.d0*kL2*Dt2*m1**2*R2**5)/(a1**6*j1**9)
*     &      *(
*     &      Om2/(2.d0*mm)*dotp(S2/S2t,cross(jv_j1,ev_e1))*(f4-2.d0*f5)
*     &      *cross(jv_j1,ev_e1)
*     &      +(f2/j1**3-Om2/mm*dotp(S2/S2t,jv_j1)*f5)*jv_j1
*     &      -Om2/(2.d0*mm)*dotp(S2/S2t,ev_e1)*f4*ev_e1
*     &      )     

*     No Orbit-average
      S2dot=(3.d0*kL2*Dt2*m1**2*R2**5)/(a1**6)*mm
     &      *(
     &      f4/j1**9*Om2/(2.d0*mm)*(S2/S2t-dotp(S2/S2t,jv_j1)*jv_j1)
     &      -f5/j1**9*Om2/mm*S2/S2t
     &      +f2/j1**12*jv_j1
     &      +dotp(ev,S2/S2t)*(6.d0+e1p2)/(4.d0*j1**9)*Om2/mm*ev
     &      )

*     Orbit-averaged
*      S1dot=-(j1*L1)/(2.d0*ta1*j1**13)
*     &     *(j1**3*f5*Om1/(2.d0*mm)*(S1/S1t+dotp(S1/S1t,jv_j1)*jv_j1)
*     &     -f2*jv_j1)
*      S2dot=-(j1*L1)/(2.d0*ta2*j1**13)
*     &     *(j1**3*f5*Om2/(2.d0*mm)*(S2/S2t+dotp(S2/S2t,jv_j1)*jv_j1)
*     &     -f2*jv_j1)

      jdot=-S2dot/L1-adot_td/(2.d0*a1)*jv

      yp(1)=yp(1)+edot(1)*Tides      
      yp(2)=yp(2)+edot(2)*Tides          
      yp(3)=yp(3)+edot(3)*Tides        
      yp(4)=yp(4)+jdot(1)*Tides          
      yp(5)=yp(5)+jdot(2)*Tides    
      yp(6)=yp(6)+jdot(3)*Tides         
      yp(16)=yp(16)+S2dot(1)*Tides    
      yp(17)=yp(17)+S2dot(2)*Tides    
      yp(18)=yp(18)+S2dot(3)*Tides    
      yp(19)=yp(19)+adot_td*Tides   
*
************************************************************************
*
*     Rotational distortion star 2
*
      
      om_r_2=kL2/2.d0*mm/(j1**4*a1**2)*Om2_n**2*R2**2

      Ops_2=-kL2*m1*R2**3*Om2/(2.d0*k2*m2*a1**3*j1**3)
     &      *dotp(jv_j1,S2/S2t)

      S2dot=Ops_2*cross_jv_j1_S2 !spin-orbit precession ok

      jdot=-S2dot/L1 ! ok
        
      yp(16)=yp(16)+S2dot(1)
      yp(17)=yp(17)+S2dot(2)
      yp(18)=yp(18)+S2dot(3)

      edot=(-om_r_2*(dotp(S2/S2t,jv_j1)*cross(S2/S2t,ev)
     &     +1.d0/2.d0*(1.d0-5.d0*dotp(S2/S2t,jv_j1)**2)
     &     *cross_jv_j1_ev))

      yp(1)=yp(1)+edot(1)*Rot      
      yp(2)=yp(2)+edot(2)*Rot       
      yp(3)=yp(3)+edot(3)*Rot     
      yp(4)=yp(4)+jdot(1)*Rot       
      yp(5)=yp(5)+jdot(2)*Rot 
      yp(6)=yp(6)+jdot(3)*Rot        
*
************************************************************************
*
*     Mass-loss star 2
* 
      S2dot=(K2*dm2*R2**2*Om2)*(S2/S2t)      
 
      yp(16)=yp(16)+S2dot(1)*Mass
      yp(17)=yp(17)+S2dot(2)*Mass
      yp(18)=yp(18)+S2dot(3)*Mass
      yp(19)=yp(19)-dm2/m12*y(19)*Mass
      yp(20)=yp(20)-dm2/m123*y(20)*Mass
*
************************************************************************
*
*     Mass-loss star 3
* 
 22   continue
      yp(20)=yp(20)-dm3/m123*y(20)*Mass
      return
      end
