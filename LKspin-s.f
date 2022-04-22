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
      external evolution,sse,vrotf,bse,sserl
      integer IR1,nq,j,h,N,iterator
      integer bhflag,kw1,kw2,kw3,true,col,CE,KCE1,KCE2,RLO
      integer kw
      integer,dimension(8) :: time_values
      parameter (nq=20,N=1000000)
      real*8 vrotf,dms1,dms2,dms3,vs(3),jorb
      real*8 vstore1,vstore2,vstore3
      real*8 y(nq),tol(nq),S1(3),S2(3),Rcr
      real*8 jv(3),jv2(3),ev(3),ev2(3),et(3),tfin0,tfin
      real*8 e1,e2,a1,a2,is,t,dt
      real*8 DT0,z,sigma
      real*8 tmax,tmax1,tmax2,tmax3
      real*8 j1,j2,Om1,Om2,Sp1,Sp2
      real*8 SevalDouble,MA,tb,om10,om20,mp0,ms0,dbm3,mass10,mass20
      real*8 fb1,fb2,fb3,epoch1,epoch2
      real*8 m1,m2,m3,R1,R2,R3
      real*8 tRmax3,tRmax
      real :: start, finish
      real scm(1000000,14),spp(20,3),scm1(1000000,14),scm2(1000000,14)
      real scm3(1000000,14)
      real*8 dtp
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
      call cpu_time(start)
      call init_random_seed()
*
************************************************************************
*
*     Constants and units
*      
      dtp = 1.e-2 
 88   continue
      pi=2.d0*ASIN(1.d0)
      c=1.d4                    !units M_Sun, AU, G=1
      ut=58.d0/365.d0           !to convert time in years
      RSun=0.00465d0            !R_Sun in AU
      col=0
      DT0=1.d-4
      DT=DT0
      t=0.d0                   
      js1=0
      js2=0
      js3=0
      CE=0
      KCE1=0
      KCE2=0
      RLO=0
      
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
*
************************************************************************
*
*     bhflag (0=no kicks, 1 momentum, 2 fallback model); natal kick model
*  
      CALL GETARG(2,arg2)
      read(arg2,*)bhflag
*
************************************************************************
*
*     nsflag (0= old, 1= Belc '02, 2= startrack, 3=rapid, 4=delayed); CO formation model
*  
      CALL GETARG(3,arg3)
      read(arg3,*)nsflag
*
************************************************************************
*
*     alpha1
*  
      CALL GETARG(4,arg4)
      read(arg4,*)alpha1

      print*,z,bhflag,nsflag,alpha1
*
************************************************************************
*
*     iterator
*  
      CALL GETARG(5,arg5)
      read(arg5,*)iterator
      print*,'IC line',iterator

c      CALL GETARG(6,arg6)
c      CALL GETARG(7,arg7)

c      directory="Plane-iso/"//trim(arg1)//"-"//trim(arg2)//"-"
c     &//trim(arg3)//"-"//trim(arg4)//"-"//trim(arg6)//"-"//trim(arg7)
c      directory=trim(directory)
c      print*,arg7,directory

      directory=trim(arg1)//"-"//trim(arg2)//"-"//trim(arg3)
     &//"-"//trim(arg4)
      directory=trim(directory)


      if(iterator.eq.0)then  
*
************************************************************************
*
*     Sample initial triple configuration, y stores all orbital params
*     y(13-18) are normalised and will be scaled with the rotation rate
*     below
*  
c         CALL initial(M1,M2,M3,y)
          goto 999
      else
c         open(unit=99,file="Plane/"//trim(arg1)//"-"//trim(arg2)//"-"
c     &//trim(arg3)//"-"//trim(arg4)//"-"//trim(arg6)//"-"//trim(arg7)
c     &//"/stable-IC.dat")
         open(unit=99,file="./stable-IC.dat")
c         open(unit=99,file=trim(arg1)//"-"//trim(arg2)//"-"//trim(arg3)
c     &//"-"//trim(arg4)//"-"//trim(arg6)//"-"//trim(arg7)
c     &//"/stable-IC.dat")
c         open(unit=99,file="Table/MdS/"//trim(directory)//
c     &"/CO-triple.dat")
         do h=2,iterator
            read(99,*)
         enddo
         read(99,*)m1,m2,m3,y
c         read(99,*)dummy,
c     &           y,dummx,
c     &           m1,m2,m3,dummx,
c     &           dummy,dummy,dummy,
c     &           dummy,dummy,dummy,
c    &           dummy,dummy
        close(99)
      endif


c      CALL GETARG(6,arg6)
c      read(arg6,*)y(20)
c      y(20)=10.d0**y(20)

c      CALL GETARG(7,arg7)
c      read(arg7,*)m3
      print*,'      m1=',m1
      print*,'      m2=',m2
      print*,'      m3=',m3
      do h=1,20
         print*,'      y(',h,')=',y(h)
      enddo
      a1=y(19)
      a2=y(20)
      e1=sqrt(y(1)**2+y(2)**2+y(3)**2)
      e2=sqrt(y(7)**2+y(8)**2+y(9)**2)
      print*,'Initial masses (3), semi-major axes (2), and ecc (2):'
      print*,m1,m2,m3
      print*,a1,a2
      print*,e1,e2

c      y(20)=-1.d0
*
************************************************************************
*
*     Compute stellar evolution using sse
*       
      print*,'evolve star 1'
      call sse(M1,N,tabt1,tabm1,tabr1,tabdm1,tabmenv1,tabrenv1,
     &     tablumin1,kwv1,j,tco1,dms1,Rmax1,fb1,z,tmax1,kw,tRmax,mco1,
     &      dtp)
      kwfin1 = kw
      scm1=scm
      cont1=j-1
c      CALL FMMsplineDouble(tabt1,tabr1,cont1)
c      CALL FMMsplineDouble(tabt1,tabdm1,cont1)
c      CALL FMMsplineDouble(tabt1,tabmenv1,p10,p11,p12,cont1)
c      CALL FMMsplineDouble(tabt1,tabrenv1,p13,p14,p15,cont1)
c      CALL FMMsplineDouble(tabt1,tablumin1,p16,p17,p18,cont1)
      print*,'evolve star 2'
      call sse(M2,N,tabt2,tabm2,tabr2,tabdm2,tabmenv2,tabrenv2,
     &     tablumin2,kwv2,j,tco2,dms2,Rmax2,fb2,z,tmax2,kw,tRmax,mco2,
     &      dtp)
      kwfin2 = kw
      scm2=scm      
      cont2=j-1     
c      CALL FMMsplineDouble(tabt2,tabm2,cont2)
c      CALL FMMsplineDouble(tabt2,tabr2,cont2)
c      CALL FMMsplineDouble(tabt2,tabdm2,cont2)
c      CALL FMMsplineDouble(tabt2,tabmenv2,d10,d11,d12,cont2)
c      CALL FMMsplineDouble(tabt2,tabrenv2,d13,d14,d15,cont2)
c      CALL FMMsplineDouble(tabt2,tablumin2,d16,d17,d18,cont2)
      print*,'evolve Tertiary'
      call sse(M3,N,tabt3,tabm3,tabr3,tabdm3,tabmenv3,tabrenv3,
     &     tablumin3,kwv3,j,tco3,dms3,Rmax3,fb3,z,tmax3,kw,tRmax,mco3,
     &      dtp)
      tRmax3 = tRmax
      kwfin3 = kw
      scm3=scm 
      cont3=j-1
c      CALL FMMsplineDouble(tabt3,tabm3,cont3) !spline evolution of M,R, and DM 
c      CALL FMMsplineDouble(tabt3,tabr3,cont3)
c      CALL FMMsplineDouble(tabt3,tabdm3,cont3)
c      CALL FMMsplineDouble(tabt3,tabmenv3,t10,t11,t12,cont3)
c      CALL FMMsplineDouble(tabt3,tabrenv3,t13,t14,t15,cont3)
c      CALL FMMsplineDouble(tabt3,tablumin3,t16,t17,t18,cont3)
      print*,'CO masses',mco1,mco2,mco3
      print*,'CO times',tco1*ut/1.d6,tco2*ut/1.d6,tco3*ut/1.d6
*
************************************************************************
*
*     Initial radii and spins SSElike
*      
      R1=Rsun*10.d0**SevalDouble(0.0d0,tabt1,tabr1,cont1) 
      R2=Rsun*10.d0**SevalDouble(0.0d0,tabt2,tabr2,cont2) 
      Om1=45.35d0*vrotf(M1)/(R1/Rsun)*ut !from Lang 1992, and Hurley et al. 2000
      Om2=45.35d0*vrotf(M2)/(R2/Rsun)*ut   
      Sp1=Om1*K1*M1*R1**2       !spin AM = moment of inertia*stellar spin frequency
      Sp2=Om2*K2*M2*R2**2
                 
c      y(13)=y(13)*Sp1 !S1x
c      y(14)=y(14)*Sp1 !S1y
c      y(15)=y(15)*Sp1 !S1z
c      
c      y(16)=y(16)*Sp2 !S2x
c      y(17)=y(17)*Sp2 !S2y
c      y(18)=y(18)*Sp2 !S2z
*
************************************************************************
*
*     Record initial set-up, all-IC.dat
*            
      y0=y
      m10=m1
      m20=m2
      m30=m3
      R10=R1
      R20=R2



*
************************************************************************
*
*     Dynamical evolution
*     
*      do while((t*ut/1.d6).lt.tmax)
      print*,'Start Integration'
      do while(.true.)
            if(t.lt.tco1)then
               m1=SevalDouble(t,tabt1,tabm1,cont1)
            else
               m1=mco1
            endif
            if(t.lt.tco2)then
               m2=SevalDouble(t,tabt2,tabm2,cont2)
            else
               m2=mco2
            endif 
            if(t.lt.tco3)then
               m3=SevalDouble(t,tabt3,tabm3,cont3)
            else
               m3=mco3
            endif 
            R1=Rsun*10.d0**SevalDouble(t,tabt1,tabr1,cont1)
            if(R1.gt.Rmax1)R1=Rmax1
            R2=Rsun*10.d0**SevalDouble(t,tabt2,tabr2,cont2)
            if(R2.gt.Rmax2)R2=Rmax2
c            print*,t*ut/1d6,y(19),(m1+m2)/(m10+m20)
*
************************************************************************
*
*     Check for supernovae and update orbital parameters
*      
         if(t.ge.tco1.and.js1.eq.0)then
            js1=1
            if(kwfin1.eq.15)then
               print*,'Mass-less SN in inner binary (1) at (Myr)',
     &               t*ut/1.d6
               print*
               col=4
               tfin=t*ut/1.d6
               goto 841
            endif
            m1=SevalDouble(tco1,tabt1,tabm1,cont1)
            call Vkick(kwfin1,y,m1,m2,m3,dms1,jorb,vs,fb1,vstore1)
            print*,'SN in inner binary (1) at (Myr)',t*ut/1.d6,m1,y(19)
            print*
            DT=DT0
            do h=1,3
               ev(h)=y(h)           
               ev2(h)=y(h+6) 
            enddo
            e1=sqrt(dotp(ev,ev))
            e2=sqrt(dotp(ev2,ev2))
            if(e1.gt.1.d0.or.e1.lt.0.d0.or.y(19).lt.0.d0)goto 563
            if(e2.gt.1.d0.or.e2.lt.0.d0.or.y(20).lt.0.d0)goto 777
         end if                        
         if(t.ge.tco2.and.js2.eq.0)then
            js2=1
            if(kwfin2.eq.15)then
               print*,'Mass-less SN in inner binary (2) at (Myr)',
     &               t*ut/1.d6
               print*
               col=4
               tfin=t*ut/1.d6
               goto 841
            endif
            m2=SevalDouble(tco2,tabt2,tabm2,cont2)
            print*,'a1',y(19)
            call Vkick(kwfin2,y,m2,m1,m3,dms2,jorb,vs,fb2,vstore2)
            print*,'SN in inner binary (2) at (Myr)',t*ut/1.d6,m2,y(19)
            print*
            DT=DT0
            do h=1,3
               ev(h)=y(h)           
               ev2(h)=y(h+6) 
            enddo
            e1=sqrt(dotp(ev,ev))
            e2=sqrt(dotp(ev2,ev2))
            if(e1.gt.1.d0.or.e1.lt.0.d0.or.y(19).lt.0.d0)goto 563
            if(e2.gt.1.d0.or.e2.lt.0.d0.or.y(20).lt.0.d0)goto 777
         end if
         if(t.ge.tco3.and.js3.eq.0)then
            js3=1
            if(kwfin3.eq.15)then
              print*,'Mass-less SN in the tertiary at (Myr)',
     &              t*ut/1.d6
              print*
              e2=-1.d0
              goto 777
            endif
            m3=SevalDouble(tco3,tabt3,tabm3,cont3) 
            call Vkick3(kwfin3,y,m1,m2,m3,dms3,jorb,vs,fb3,vstore3)
            print*,'SN in the tertiary at (Myr)',t*ut/1.d6,m3,y(20)
            print*
            DT=DT0       
            do h=1,3
               ev(h)=y(h)           
               ev2(h)=y(h+6) 
            enddo
            e1=sqrt(dotp(ev,ev))
            e2=sqrt(dotp(ev2,ev2))
            if(e1.gt.1.d0.or.e1.lt.0.d0.or.y(19).lt.0.d0)goto 563
            if(e2.gt.1.d0.or.e2.lt.0.d0.or.y(20).lt.0.d0)goto 777
         end if   

         do h=1,3
            ev(h)=y(h)           
            ev2(h)=y(h+6)  
            jv(h)=y(h+3)        
            S1(h)=y(12+h)
            S2(h)=y(15+h)
         end do
         a1=y(19)
         a2=y(20)
         e1=sqrt(dotp(ev,ev))
         e2=sqrt(dotp(ev2,ev2))
         j1=sqrt(dotp(jv,jv))
         Sp1=sqrt(dotp(S1,S1))
         Sp2=sqrt(dotp(S2,S2))
*
************************************************************************
*
*     Check if inner binary is broken apart
*     If yes, terminate program and record result
*         
 563   continue
       if(e1.gt.1.d0.or.e1.lt.0.d0.or.y(19).lt.0.d0)then
            print*,'Inner binary broken apart at (Myr)',
     &           t*ut/1.d6,e1,e2,y(19),y(20)   
      open(61,file=trim(directory)//"/broken.dat",
     &action='write',access='append')
            j=0
            do while(scm1(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw1=scm1(j,2)
            j=0
            do while(scm2(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw2=scm2(j,2)
            j=0
            do while(scm3(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw3=scm3(j,2)
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
            close(61)
            goto 999
       end if 
*
************************************************************************
*
*     Check if outer binary is broken apart
*     If yes, finish inner binary with BSE
* 
 777   continue        
       if(e2.gt.1.d0.or.e2.lt.0.d0.or.y(20).lt.0.d0)then
            print*,'Outer binary broken apart at (Myr)',
     &           t*ut/1.d6,e1,e2,y(19),y(20)   
            print*,'Finish evolution of the inner binary with BSE'
            j=0
            do while(scm1(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            tb=2.d0*pi*sqrt(y(19)**3/(m1+m2))*58.d0
            mass10=scm1(j,3)
            mass20=scm2(j,3)
            kw1=scm1(j,2)
            kw2=scm2(j,2)
            om1=Sp1/(K1*M1*R1**2)/ut
            om2=Sp2/(K2*M2*R2**2)/ut
            epoch1=scm1(j,12)   
            tfin=t*ut/1.d6
            dbm3=-10.d0**SevalDouble(t,tabt3,tabdm3,cont3)
            a2=-1.d0 ! indicate that there is no longer a tertiary companion
            call bse(mass10,mass20,
     &           m1,m2,om1,om2,kw1,kw2,z,tb,tfin,
     &           e1,a2,a1,m3,epoch1,epoch2,true,
     &           col,dbm3,min(tmax1,tmax2),CE,KCE1,KCE2,RLO) !evolve the binary until BSE stops
            t=tfin/(ut/1.d6)
            kwfin1 = kw1
            kwfin2 = kw2
            do h=1,3
               y(h)=y(h)*e1/sqrt(dotp(ev,ev))
               y(h+3)=y(h+3)*sqrt(1.d0-e1**2)/j1  
               y(12+h)=-1.d0 ! no treatment for now
               y(15+h)=-1.d0 ! no treatment for now
            end do  
            y(19)=a1
            y(20)=-1.d0
            if(kw1.gt.10)js1=1
            if(kw2.gt.10)js2=1
            goto 841
       end if
*
************************************************************************
*
*     Check if system is broken apart
*     If yes, terminate program and record result
*         
*       if(e1.gt.1.d0.or.e2.gt.1.d0.or.e1.lt.0.d0.or.e2.lt.0.d0.
*     &        or.y(19).lt.0.d0.or.y(20).lt.0.d0)then
*            print*,'System broken apart at (Myr)',
*     &           t*ut/1.d6,e1,e2,y(19),y(20)   
*             open(61,file=trim(arg1)//"-"//trim(arg2)//"-"//trim(arg3)//"-"//trim(arg4)//"-"//trim(arg6)//"-"//trim(arg7)

*     &           action='write',access='append')
*            write(61,*)t*ut/1.d6,
*     &           y0,y,ySN1,ySN2,ySN3,
*     &           m10,m20,m30,m1,m2,m3,m1SN1,m2SN1,m3SN1,
*     &           m1SN2,m2SN2,m3SN2,m1SN3,m2SN3,m3SN3,
*     &           R10,R20,R1,R2,R1SN1,R2SN1,R1SN2,R2SN2,R1SN3,R2SN3,
*     &           kwfin1,kwfin2,kwfin3,tco1,tco2,tco3,js1,js2,js3,col,CE
*            close(61)
*            goto 999
*       end if     
         
!     Check for RLO, if yes evolve binary with BSE
         if(R1.gt.y(19)*(1.d0-e1)*Roche(m1,m2)
     &        .and.col.eq.0.and.js1.ne.1)then            
            print*,'RLO in the inner binary, comp.1 at (Myr)',t*ut/1.d6
            j=0
            do while(scm1(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            tb=2.d0*pi*sqrt(y(19)**3/(m1+m2))*58.d0
            mass10=scm1(j,3)
            mass20=scm2(j,3)
            kw1=scm1(j,2)
            kw2=scm2(j,2)
            om1=Sp1/(K1*M1*R1**2)/ut
            om2=Sp2/(K2*M2*R2**2)/ut
            epoch1=scm1(j,12)   
            tfin=t*ut/1.d6
            dbm3=-10.d0**SevalDouble(t,tabt3,tabdm3,cont3)

            call bse(mass10,mass20,
     &           m1,m2,om1,om2,kw1,kw2,z,tb,tfin,
     &           e1,a2,a1,m3,epoch1,epoch2,true,
     &           col,dbm3,min(tmax1,tmax2),CE,KCE1,KCE2,RLO) !evolve the binary until RLO stops
            if(true==0)goto 55  !bse did not find any RLO, so nothing is done
            kwfin1 = kw1
            kwfin2 = kw2
            if(col.eq.1.or.col.eq.2.or.col.eq.4)goto 841
            if(kw1.gt.10)js1=1
            if(kw2.gt.10)js2=1

!     continue the evolution of the stars as isolated
            tfin0=tfin
            om10=om1
            om20=om2
            mp0=m1
            ms0=m2          
            print*,'evolve star 1 after RLO'
            call sserl(m1,om1,tfin0,kw1,mass10, 
     &           N,tabt1,tabm1,tabr1,tabdm1,tabmenv1,tabrenv1,
     &           tablumin1,kwv1,j,tco1,dms1,Rmax1,fb1,z,epoch1,tmax1,
     &           mco1,dtp)
            kwfin1 = kw1
            scm1=scm
            cont1=j-1
c      c      CALL FMMsplineDouble(tabt1,tabm1,cont1) !spline evolution of M,R, and DM 
c      c      CALL FMMsplineDouble(tabt1,tabr1,cont1)
c      c      CALL FMMsplineDouble(tabt1,tabdm1,cont1)
c      c      CALL FMMsplineDouble(tabt1,tabmenv1,p10,p11,p12,cont1)
c      c      CALL FMMsplineDouble(tabt1,tabrenv1,p13,p14,p15,cont1)
c      c      CALL FMMsplineDouble(tabt1,tablumin1,p16,p17,p18,cont1)
            WRITE(*,*)
            print*,'evolve star 2 after RLO'
            call sserl(m2,om2,tfin0,kw2,mass20, 
     &           N,tabt2,tabm2,tabr2,tabdm2,tabmenv2,tabrenv2,
     &           tablumin2,kwv2,j,tco2,dms2,Rmax2,fb2,z,epoch2,tmax2,
     &           mco2,dtp)
            kwfin2 = kw2
            if(kwfin1.lt.13.or.kwfin2.lt.13)then
              print*,'No CO remnants in the inner binary'
              goto 841
             end if
            scm2=scm      
            cont2=j-1      
c      c      CALL FMMsplineDouble(tabt2,tabm2,cont2) !spline evolution of M,R, and DM 
c      c      CALL FMMsplineDouble(tabt2,tabr2,cont2)
c      c      CALL FMMsplineDouble(tabt2,tabdm2,cont2)
c      c      CALL FMMsplineDouble(tabt2,tabmenv2,d10,d11,d12,cont2)
c      c      CALL FMMsplineDouble(tabt2,tabrenv2,d13,d14,d15,cont2)
c      c      CALL FMMsplineDouble(tabt2,tablumin2,d16,d17,d18,cont2)
            t=tfin0/ut*1.d6     !time RLO ends
!     Change vectors             
            
            R1=Rsun*10.d0**SevalDouble(t,tabt1,tabr1,cont1)
            R2=Rsun*10.d0**SevalDouble(t,tabt2,tabr2,cont2)  
            R3=Rsun*10.d0**SevalDouble(t,tabt3,tabr3,cont3)
            y(19)=a1
            y(20)=a2
            m1=mp0
            m2=ms0
            m3=SevalDouble(t,tabt3,tabm3,cont3)           
            do h=1,3
               y(h)=y(h)*e1/sqrt(dotp(ev,ev))
               y(h+3)=y(h+3)*sqrt(1.d0-e1**2)/j1  
               y(12+h)=y(12+h)*om10*ut*(K1*M1*R1**2)/Sp1
               y(15+h)=y(15+h)*om20*ut*(K2*M2*R2**2)/Sp2
            end do        
            DT=DT0
            print*,'after RLO',a1,e1,
     &R1/y(19)*(1.d0-e1)*Roche(m1,m2)
 55         continue
         end if
         
         if(R2.gt.y(19)*(1.d0-e1)*Roche(m2,m1)
     &        .and.col.eq.0.and.js2.ne.1)then            
            print*,'RLO in the inner binary, comp. 2 at (Myr)',
     &            t*ut/1.d6
            j=0
            do while(scm1(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            tb=2.d0*pi*sqrt(y(19)**3/(m1+m2))*58.d0
            mass10=scm1(j,3)
            mass20=scm2(j,3)
            kw1=scm1(j,2)
            kw2=scm2(j,2)
            om1=Sp1/(K1*M1*R1**2)/ut
            om2=Sp2/(K2*M2*R2**2)/ut
            epoch1=scm1(j,12)   
            epoch2=scm2(j,12)   
            tfin=t*ut/1.d6
            dbm3=-10.d0**SevalDouble(t,tabt3,tabdm3,cont3)
            call bse(mass10,mass20,m1,m2,om1,om2,
     &           kw1,kw2,z,tb,tfin,
     &           e1,a2,a1,m3,epoch1,epoch2,true,
     &           col,dbm3,min(tmax1,tmax2),CE,KCE1,KCE2,RLO) !evolve the binary until RLO stops
            if(true==0)goto 66
            kwfin1 = kw1
            kwfin2 = kw2
            if(col.eq.1.or.col.eq.2.or.col.eq.4)goto 841
            if(kw1.gt.12)js1=1
            if(kw2.gt.12)js2=1
!     continue the evolution of the stars as isolated
            tfin0=tfin
            om10=om1
            om20=om2
            mp0=m1
            ms0=m2
            print*,'evolve star 1 after RLO'
            call sserl(m1,om1,tfin0,kw1,mass10, 
     &           N,tabt1,tabm1,tabr1,tabdm1,tabmenv1,tabrenv1,
     &           tablumin1,kwv1,j,tco1,dms1,Rmax1,fb1,z,epoch1,tmax1,
     &           mco1,dtp)    
            kwfin1 = kw1    
            scm1=scm
            cont1=j-1
c      c      CALL FMMsplineDouble(tabt1,tabm1,cont1) !spline evolution of M,R, and DM 
c      c      CALL FMMsplineDouble(tabt1,tabr1,cont1)
c      c      CALL FMMsplineDouble(tabt1,tabdm1,cont1)
c      c      CALL FMMsplineDouble(tabt1,tabmenv1,p10,p11,p12,cont1)
c      c      CALL FMMsplineDouble(tabt1,tabrenv1,p13,p14,p15,cont1)
c      c      CALL FMMsplineDouble(tabt1,tablumin1,p16,p17,p18,cont1)
            WRITE(*,*)
            print*,'evolve star 2 after RLO'
            call sserl(m2,om2,tfin0,kw2,mass20, 
     &           N,tabt2,tabm2,tabr2,tabdm2,tabmenv2,tabrenv2,
     &           tablumin2,kwv2,j,tco2,dms2,Rmax2,fb2,z,epoch2,tmax2,
     &           mco2,dtp)
            kwfin2 = kw2
            if(kwfin1.lt.13.or.kwfin2.lt.13)then
              print*,'No CO remnants in the inner binary'
              goto 841
             end if
            scm2=scm      
            cont2=j-1      
c      c      CALL FMMsplineDouble(tabt2,tabm2,cont2) !spline evolution of M,R, and DM 
c      c      CALL FMMsplineDouble(tabt2,tabr2,cont2)
c      c      CALL FMMsplineDouble(tabt2,tabdm2,cont2)
c      c      CALL FMMsplineDouble(tabt2,tabmenv2,d10,d11,d12,cont2)
c      c      CALL FMMsplineDouble(tabt2,tabrenv2,d13,d14,d15,cont2)
c      c      CALL FMMsplineDouble(tabt2,tablumin2,d16,d17,d18,cont2)
            t=tfin0/ut*1.d6     !time RLO ends
!     Change vectors             
            
            R1=Rsun*10.d0**SevalDouble(t,tabt1,tabr1,cont1)
            R2=Rsun*10.d0**SevalDouble(t,tabt2,tabr2,cont2)    
            R3=Rsun*10.d0**SevalDouble(t,tabt3,tabr3,cont3)
            y(19)=a1
            y(20)=a2
            m1=mp0
            m2=ms0
            m3=SevalDouble(t,tabt3,tabm3,cont3)
            do h=1,3
               y(h)=y(h)*e1/sqrt(dotp(ev,ev))
               y(h+3)=y(h+3)*sqrt(1.d0-e1**2)/j1  
               y(12+h)=y(12+h)*om10*ut*(K1*M1*R1**2)/Sp1
               y(15+h)=y(15+h)*om20*ut*(K2*M2*R2**2)/Sp2
            end do        
            DT=DT0
 66         continue
         end if

         do h=1,3
            ev(h)=y(h)           
            ev2(h)=y(h+6)     
         end do
         a1=y(19)
         a2=y(20)
         e1=sqrt(dotp(ev,ev))
         e2=sqrt(dotp(ev2,ev2))

c         print*,t*ut/1.d6,a1,e1

         if(a1.lt.1d-3.and.js1.eq.1.and.js2.eq.1)then
            print*,a1,'Merger in LK'
            col=3
            goto 841
         end if 

         if(isnan(y(1)))then
              print*,"Code breaks",t*ut/1.d6,m1,m2,m3,R1,R2,y
              if(dtp.ge.1e-5)then
               dtp=dtp/10.d0
               scm1 = 0.d0
               scm2 = 0.d0
               scm3 = 0.d0
               print*,'dtp',dtp
               goto 88
              endif
      open(61,file=trim(directory)//"/breaks.dat",
     &action='write',access='append')
            j=0
            do while(scm1(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw1=scm1(j,2)
            j=0
            do while(scm2(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw2=scm2(j,2)
            j=0
            do while(scm3(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw3=scm3(j,2)
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
               close(61)
               goto 999
         end if

         MA=2.8d0*((1.d0+M3/(M1+M2))*(1.d0+e2)/sqrt(1.d0-e2))
     &**(2.d0/5.d0) !Mardling+Aarseth
         if(y(20)*(1.d0-e2).lt.MA*y(19))then
            print*,'System dynamically unstable at (Myr)',t*ut/1.d6
            open(61,file=trim(directory)//"/unstable.dat",
     &action='write',access='append')
            j=0
            do while(scm1(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw1=scm1(j,2)
            j=0
            do while(scm2(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw2=scm2(j,2)
            j=0
            do while(scm3(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw3=scm3(j,2)
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
            close(61)
            goto 999
         end if

         R3=Rsun*10.d0**SevalDouble(t,tabt3,tabr3,cont3)
         if(R3.gt.Rmax3)R3=Rmax3   !avoid overpredicting R 
c         if(y(20)*(1.d0-sqrt(dotp(ev2,ev2)))-R3.lt.
c     &     max(R1,R2)+y(19)*(1+sqrt(dotp(ev,ev)))+Rcr)then  
         if(R3.gt.y(20)*(1.d0-e2)*Roche(m3,m1+m2))then
            print*,'RLO in the outer binary at (Myr)',t*ut/1.d6      
            open(61,file=trim(directory)//"/3-RLO.dat",
     &action='write',access='append')
            j=0
            do while(scm1(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw1=scm1(j,2)
            j=0
            do while(scm2(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw2=scm2(j,2)
            j=0
            do while(scm3(j,1).le.t*ut/1.d6)
               j=j+1     
            end do
            kw3=scm3(j,2)
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
            close(61)
            goto 999
         end if             
         
!     evolve orbit
c         if((t-time_counter)*ut/1.d6.ge.0.01d0)then
c            time_counter=t
c        do h=1,3
c            ev(h)=y(h)           
c            ev2(h)=y(h+6)     
c         end do
c         a1=y(19)
c         a2=y(20)
c         e1=sqrt(dotp(ev,ev))
c         e2=sqrt(dotp(ev2,ev2))
c         open(44,file="evolution.dat",action='write',access='append')
c         write(44,*)t*ut/1.d6,m1,m2,m3,e1,e2,a1,a2,R1,R2,R3
c         close(44)
         call RK78(IR1,t,dt,y,nq,tol,evolution)
         if(t.ge.tmax/(ut/1.d6).or.
     &(js1.eq.1.and.js2.eq.1.and.js3.eq.1.and.
     &   t.gt.tco1*1.01d0.and.t.gt.tco2*1.01d0.and.
     &   t.gt.tco3*1.01d0).or.
     &(js1.eq.1.and.js2.eq.1.and.js3.eq.0.and.
     &   kwfin3.ne.13.and.kwfin3.ne.14.and.
     $   t.gt.tco1*1.01d0.and.t.gt.tco2*1.01d0.and.
     &   t.gt.tRmax3*1.01d0))then
            print*,'Integration time reached (Myr)',t*ut/1.d6
            goto 841
         elseif(((t+dt)*ut/1.d6).gt.tmax)then
            print*,'hit1',t,dt
            dt=tmax/(ut/1.d6)-t ! make sure that last-step doesn't exceed tmax
            print*,'hit2',t,dt
         endif
      end do

*      print*,"Integration time reached",tmax
 841  continue
      call cpu_time(finish)
      print*,'Time =',finish-start,' seconds.'
      print*,'RLO',RLO
      print*,'CE',CE
      print*,'KCE1',KCE1
      print*,'KCE2',KCE2
      j=0
      do while(scm1(j,1).le.t*ut/1.d6)
         j=j+1     
      end do
      kw1=scm1(j,2)
      j=0
      do while(scm2(j,1).le.t*ut/1.d6)
         j=j+1     
      end do
      kw2=scm2(j,2)
      j=0
      do while(scm3(j,1).le.t*ut/1.d6)
         j=j+1     
      end do
      kw3=scm3(j,2)
      if(col.eq.3)then
       print*,'CO-CO Merger at (Myr)',t*ut/1.d6
       print*,m1,m2,y(19),sqrt(y(1)**2+y(2)**2+y(3)**2)
       open(61,file=trim(directory)//"/merger.dat",
     &action='write',access='append')
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
       close(61)
      elseif(col.eq.4)then
       print*,'Mass-less remnant at (Myr)',tfin
       open(61,file=trim(directory)//
     &"/Massless-remnant.dat",action='write',access='append')
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
       close(61)
      elseif(col.eq.1)then       
         print*,'stellar collision at (Myr)',tfin
         open(61,file=trim(directory)//"/collision.dat",
     &action='write',access='append')
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
         close(61)
      elseif(col.eq.2)then       
         print*,'Disruption at (Myr)',tfin
         open(61,file=trim(directory)//"/disrupt.dat",
     &action='write',access='append')
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
         close(61)
      elseif(a2.eq.-1.d0.and.js1.eq.1.and.js2.eq.1)then
       print*,'Stable CO-binary, without tertiary companion',tfin
       open(61,file=trim(directory)//
     &"/CO-binary-alone.dat",action='write',access='append')
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
       close(61)
      elseif(a2.eq.-1.d0)then
       print*,'No CO-binary in the centre, without tertiary companion',
     &         tfin
       open(61,file=trim(directory)//
     &"/no-CO-binary-alone.dat",action='write',access='append')
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
       close(61)
      elseif(js1.eq.1.and.js2.eq.1.and.js3.eq.1)then
       print*,'stable CO-triple'
       open(61,file=trim(directory)//"/CO-triple.dat",
     &action='write',access='append')
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
       close(61)
      elseif(js1.eq.1.and.js2.eq.1.and.js3.ne.1)then
       print*,'stable CO-binary plus star'
       open(61,file=trim(directory)//"/CO-binary.dat",
     &action='write',access='append')
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
       close(61)
      else
       print*,'No CO-binary in the centre'
       open(61,file=trim(directory)//"/no-CO-binary.dat",
     &action='write',access='append')
       write(61,*)t*ut/1.d6,
     &           y0,y,
     &           m10,m20,m30,m1,m2,m3,
     &           R10,R20,R1,R2,
     &           kwfin1,kwfin2,kwfin3,
     &           tco1,tco2,tco3,
     &           js1,js2,js3,
     &           col,CE,KCE1,KCE2,
     &           vstore1,vstore2,vstore3,
     &           fb1,fb2,fb3,
     &           RLO,iterator,
     &           kw1,kw2,kw3
       close(61)
      end if
 999  continue
      end program spinTriple
