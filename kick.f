***
      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,vs)
      implicit none
*
      integer kw,k
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      INTEGER jkick
      integer bhflag
      real*8 m1,m2,m1n,ecc,sep,jorb,ecc2,harv
      real*8 pi,twopi,gmrkm,yearsc,rsunkm,mco
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mm,em,dif,der,del,r,ps(4),mcmax
      real*8 u1,u2,vk,v(4),s,theta,phi
      real*8 sphi,cphi,stheta,ctheta,salpha,calpha
      real*8 vr,vr2,vk2,vn2,hn2,fb
      real*8 mu,cmu,vs(3),v1,v2,mx1,mx2,angle
      real*8 sigma,sigmaS
      real*8 sepBlaauw,eccBlaauw
      real*8 vk_old,sep_old
      COMMON /VALUE4/ sigma,bhflag
      COMMON /VALUEJ/ angle
      common /fall/fb
      COMMON/mco/mcmax
      real ran3,xx
      external ran3
*
c      print*,'Before kick',sep,ecc
      sep_old=sep
      mco=mcmax

      do k = 1,3
         vs(k) = 0.d0
      enddo
*     if(kw.eq.14.and.bhflag.eq.0) goto 95
*
c      if(m1n.gt.3.)then        
c         sigmaS=sigma*1.5/m1n
c      else
c         sigmaS=sigma
c      end if
      
      
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
      gmrkm = 1.906125d+05
*
* Find the initial separation by randomly choosing a mean anomaly.
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
         call RANDOM_NUMBER(harv)
         xx =harv !RAN3(idum)!* RAN3(idum)
         mm = xx*twopi
         em = mm
 2       dif = em - ecc*SIN(em) - mm
         if(ABS(dif/mm).le.1.0d-04) goto 3
         der = 1.d0 - ecc*COS(em)
         del = dif/der
         em = em - del
         goto 2
 3       continue
         r = sep*(1.d0 - ecc*COS(em))
         do 22 k = 1,2
         u1 = RAN3(idum)
         u2 = RAN3(idum)
* Generate two velocities from polar coordinates S & THETA.
         s = r
         theta = twopi*u2
         ps(2*k-1) = s*COS(theta)
         ps(2*k) = s*SIN(theta)
 22   continue
         
*
* Find the initial relative velocity vector.
         salpha = SQRT((sep*sep*(1.d0-ecc*ecc))/(r*(2.d0*sep-r)))
         calpha = (-1.d0*ecc*SIN(em))/SQRT(1.d0-ecc*ecc*COS(em)*COS(em))
         vr2 = gmrkm*(m1+m2)*(2.d0/r - 1.d0/sep)
         vr = SQRT(vr2)
      else
         vr = 0.d0
         vr2 = 0.d0
         salpha = 0.d0
         calpha = 0.d0
      endif
*
* Generate Kick Velocity using Maxwellian Distribution (Phinney 1992).
* Use Henon's method for pairwise components (Douglas Heggie 22/5/97).
c      do 20 k = 1,2
c         u1 = RAN3(idum)
c         u2 = RAN3(idum)
* Generate two velocities from polar coordinates S & THETA.
c         s = sigmaS*SQRT(-2.d0*LOG(1.d0 - u1))
c         theta = twopi*u2
c         v(2*k-1) = s*COS(theta)
c         v(2*k) = s*SIN(theta)
c 20   continue
c      vk2 = v(1)**2 + v(2)**2 + v(3)**2
c      vk = SQRT(vk2)
     
!     kick velcoity
      if(bhflag.eq.1)then       !momoentum conserving kicks            
         if(m1n.gt.3.d0)then        
            sigmaS=sigma*1.5d0/m1n !convert unit below
         else
            sigmaS=sigma
         end if
         do 21 k = 1,2
            call RANDOM_NUMBER(harv)
            u1 = harv
            call RANDOM_NUMBER(harv)
            u2 = harv
*     Generate two velocities from polar coordinates S & THETA.
            s = sigmaS*DSQRT(-2.d0*DLOG(1.d0 - u1))
            theta = twopi*u2
            v(2*k-1) = s*DCOS(theta)
            v(2*k) = s*DSIN(theta)
 21      continue
         vk2 = v(1)**2 + v(2)**2 + v(3)**2
         vk = DSQRT(vk2)
      end if

      if(bhflag.eq.2)then       ! kicks
         do 420 k = 1,2
            call RANDOM_NUMBER(harv)
            u1 = harv
            call RANDOM_NUMBER(harv)
            u2 = harv
*     Generate two velocities from polar coordinates S & THETA.
            s = sigma*SQRT(-2.d0*LOG(1.d0 - u1))
            theta = twopi*u2
            v(2*k-1) = s*COS(theta)
            v(2*k) = s*SIN(theta)
 420     continue
         vk2 = v(1)**2 + v(2)**2 + v(3)**2
         vk = SQRT(vk2)
         vk=(1.d0-fb)*vk
         vk2=vk**2
      end if
      
      if(bhflag.eq.0)then
         vk2 = 0.d0
         vk = 0.d0
      endif

      if((kw.eq.14.and.bhflag.eq.0).or.kw.lt.0)then
         vk2 = 0.d0
         vk = 0.d0
         if(kw.lt.0) kw = 13         
      endif
      
      if(kw.eq.13)then       ! full NS kicks
         print*,'NS kick'
         do 520 k = 1,2
            call random_number(harv)
            u1 = harv
            call random_number(harv)
            u2 = harv
*     Generate two velocities from polar coordinates S & THETA.
            s = sigma*SQRT(-2.d0*LOG(1.d0 - u1))
            theta = twopi*u2
            v(2*k-1) = s*COS(theta)
            v(2*k) = s*SIN(theta)
 520     continue
         vk2 = v(1)**2 + v(2)**2 + v(3)**2
         vk = SQRT(vk2)
      end if 
          
      sphi = -1.d0 + 2.d0*u1
      phi = ASIN(sphi)
      cphi = COS(phi)
      stheta = SIN(theta)
      ctheta = COS(theta)
      
*     WRITE(66,*)' KICK VK PHI THETA ',vk,phi,theta
      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90
***********************************************************************
*
* Determine the magnitude of the new relative velocity.
      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha-stheta*cphi*calpha)

* Calculate the new semi-major axis.
      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
      
      sep = 1.d0/sep
c      print*,sep*0.0045
*     if(sep.le.0.d0)then
*        ecc = 1.1d0
*        goto 90
*     endif
* Determine the magnitude of the cross product of the separation vector
* and the new relative velocity.
      v1 = vk2*sphi*sphi
      v2 = (vk*ctheta*cphi-vr*salpha)**2
      hn2 = r*r*(v1 + v2)
* Calculate the new eccentricity.
      ecc2 = 1.d0 - hn2/(gmrkm*sep*(m1n+m2))
      ecc2 = MAX(ecc2,0.d0)
      ecc = SQRT(ecc2)
* Calculate the new orbital angular momentum taking care to convert
* hn to units of Rsun^2/yr.
      jorb = (m1n*m2/(m1n+m2))*SQRT(hn2)*(yearsc/rsunkm)
* Determine the angle between the new and old orbital angular
* momentum vectors.
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      mu = ACOS(cmu)
c      print*,'ang',mu*180.d0/pi
* Calculate the components of the velocity of the new centre-of-mass.
 90   continue
      if(ecc.le.1.0)then
* Calculate the components of the velocity of the new centre-of-mass.
         mx1 = vk*m1n/(m1n+m2)
         mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
         vs(1) = mx1*ctheta*cphi + mx2*salpha
         vs(2) = mx1*stheta*cphi + mx2*calpha
         vs(3) = mx1*sphi
      else
* Calculate the relative hyperbolic velocity at infinity (simple method).
c         sep = r/(ecc-1.d0)
*        cmu = SQRT(ecc-1.d0)
*        mu = ATAN(cmu)
c         mu = ACOS(1.d0/ecc)
c         vr2 = gmrkm*(m1n+m2)/sep
c         vr = SQRT(vr2)
c         vs(1) = vr*SIN(mu)
c         vs(2) = vr*COS(mu)
c         vs(3) = 0.d0
c         ecc = MIN(ecc,99.99d0)
      endif
      RETURN
      END
***
