***
      SUBROUTINE kick(kw,m1,m1n,m2,ecc,sep,jorb,vs)

      implicit none
*
      integer k
      INTEGER idum
      COMMON /VALUE3/ idum
      INTEGER idum2,iy,ir(32)
      COMMON /RAND3/ idum2,iy,ir
      integer bhflag
      integer, intent(inout) :: kw
      real*8, intent(inout) :: m1,m1n,m2,jorb
      real*8, intent(inout) :: ecc,sep,vs(3)
      real*8 ecc2,mrem,mean_ej
      real*8 pi,twopi,gmrkm,yearsc,rsunkm
      parameter(yearsc=3.1557d+07,rsunkm=6.96d+05)
      real*8 mm,em,dif,der,del,r
      real*8 u1,u2,vk,v(4),s,theta,phi,ss
      real*8 sphi,cphi,stheta,ctheta,salpha,calpha
      real*8 vr,vr2,vk2,vn2,hn2
      real*8 mu,cmu,v1,v2,mx1,mx2
      real*8 maxwellian
      real*8 sigma1,sigma2

      COMMON /VALUE4/ sigma1,sigma2,bhflag
      real ran3,xx
      external ran3
*
      INTEGER ceflag,tflag,ifflag,nsflag,wdflag,piflag
      COMMON /FLAGS/ ceflag,tflag,ifflag,nsflag,wdflag,piflag
*
      INTEGER directcollapse,ECS
      real*8 ffb,mfin
      real*8 ev(3),jv(3)
      COMMON /KICKSN/ ffb,directcollapse,ECS,mfin,ev,jv
*
* J. S.
      real*8 dotp
      real*8 ux(3),uy(3),uz(3),rv(3),cross(3),vn(3)
*
*     Orbital reference frame
      uz = jv/DSQRT(jv(1)**2+jv(2)**2+jv(3)**2)
      uy(1) = 1.d0
      uy(2) = 1.d0
      uy(3) = -(uy(1)*uz(1)+uy(2)*uz(2))/uz(3)
      uy = uy/DSQRT(uy(1)**2+uy(2)**2+uy(3)**2)
      cross(1) = uy(2) * uz(3) - uy(3) * uz(2)
      cross(2) = uy(3) * uz(1) - uy(1) * uz(3)
      cross(3) = uy(1) * uz(2) - uy(2) * uz(1)
      ux = cross
*
      do k = 1,3
         vs(k) = 0.d0
      enddo
*     if(kw.eq.14.and.bhflag.eq.0) goto 95
*
* To take into account the ECS
      if(ECS.eq.1)then
* Check if the supernova is driven by electron-cupture... 
          maxwellian = sigma2   
      else
* ... or iron core collapse. 
          maxwellian = sigma1
      endif
*
      pi = ACOS(-1.d0)
      twopi = 2.d0*pi
* Conversion factor to ensure velocities are in km/s using mass and
* radius in solar units.
      gmrkm = 1.906125d+05
*
* Find the initial separation by randomly choosing a mean anomaly.
      if(sep.gt.0.d0.and.ecc.ge.0.d0)then
         xx = RAN3(idum)
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
*
* J. S.
c         f = DATAN(DSQRT((1.d0+ecc)/(1.d0-ecc))*DTAN(em/2.d0))*2.d0       
c         jorb = DSQRT(jv(1)**2+jv(2)**2+jv(3)**2)
c         cross(1) = ux(2) * uy(3) - ux(3) * uy(2)
c         cross(2) = ux(3) * uy(1) - ux(1) * uy(3)
c         cross(3) = ux(1) * uy(2) - ux(2) * uy(1)
c         p2 = cross
c         rv = (DCOS(f)*ux+DSIN(f)*uy)*r
         rv = r*uy
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
* Use Henon''s method for pairwise components (Douglas Heggie 22/5/97).
      do 20 k = 1,2
         u1 = RAN3(idum)
         u2 = RAN3(idum)
* Generate two velocities from polar coordinates S & THETA.
* way to avoid negative values in LOG
         ss = 1.d0 - u1
         if(ss.le.0.d0)then
            ss = 1.d-6
         endif
         s = maxwellian*SQRT(-2.d0*LOG(ss))
         theta = twopi*u2
         v(2*k-1) = s*COS(theta)
         v(2*k) = s*SIN(theta)
 20   continue
      vk2 = v(1)**2 + v(2)**2 + v(3)**2
      vk = SQRT(vk2)
*
***
* Module kick velocity for Black Holes...
* Old case:
      if(kw.eq.14.and.bhflag.eq.1)then
         mrem = m1n
         vk = vk*(1.d0 - ffb)
         vk2 = vk*vk
*     ... if they are generated via direct collapse vk = 0
***
* module vick velocity for Black Holes and Neutron Stars
* Fallback case:
      elseif((kw.eq.14.or.kw.eq.13).and.(bhflag.eq.2))then
         mrem = m1n
         vk = vk*(1.d0 - ffb)
         vk2 = vk*vk
***
* Generate Kick Velocity for Black Holes and Neutron Stars
* Giacobbo&Mapelli 2020 case:
      elseif((kw.eq.14.or.kw.eq.13).and.(bhflag.eq.3))then
         mrem = m1n
*
         if(kw.eq.14 .and. ffb.eq.1.d0)then
            mfin = mrem
         endif  
* mean mass ejected <mej> in Msun (in a population of 10^6 single stars)
         if(nsflag.eq.2) mean_ej = 9.d0
         if(nsflag.eq.3) mean_ej = 7.5d0
         vk = vk*((mfin-mrem)/mean_ej)*(1.2d0/mrem)
         vk2 = vk*vk
***   
* NO-modulation of kick kelocity for Black Holes and Neutron Stars
* Full kick case:
      elseif((kw.eq.14.or.kw.eq.13).and.(bhflag.eq.4))then
         mrem = m1n
         vk = vk
         vk2 = vk*vk  
***
      elseif((kw.eq.14.and.bhflag.eq.0).or.kw.lt.0)then
         vk2 = 0.d0
         vk = 0.d0
         if(kw.lt.0) kw = 13
***     
      endif
***
* Correction (Noticed by Peter & Mirek)      
      u1 = RAN3(idum)
      u2 = RAN3(idum)
      theta = twopi*u2
***
      sphi = -1.d0 + 2.d0*u1
      phi = ASIN(sphi)
      cphi = COS(phi)
      stheta = SIN(theta)
      ctheta = COS(theta)
***
      if(sep.le.0.d0.or.ecc.lt.0.d0) goto 90
*
* Determine the magnitude of the new relative velocity.
* (NG 12/17: there was a wrong sign in the brackets)
      vn2 = vk2+vr2-2.d0*vk*vr*(ctheta*cphi*salpha+stheta*cphi*calpha)
*
* J. S.
      vn=(vk*ctheta*cphi-vr*salpha)*ux
     &      +(vk*stheta*cphi-vr*calpha)*uy
     &      +(vk*sphi)*uz
      PRINT*,DSQRT(vn2)
* Calculate the new semi-major axis.
      sep = 2.d0/r - vn2/(gmrkm*(m1n+m2))
      sep = 1.d0/sep
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
*
* J. S.
      cross(1) = rv(2) * vn(3) - rv(3) * vn(2)
      cross(2) = rv(3) * vn(1) - rv(1) * vn(3)
      cross(3) = rv(1) * vn(2) - rv(2) * vn(1)
      jv = cross
      cross(1) = vn(2) * jv(3) - vn(3) * jv(2)
      cross(2) = vn(3) * jv(1) - vn(1) * jv(3)
      cross(3) = vn(1) * jv(2) - vn(2) * jv(1)
      cross = cross
      jv = jv/DSQRT(jv(1)**2+jv(2)**2+jv(3)**2)*DSQRT(1.d0-ecc2)
      ev = cross/(m1n+m2)/gmrkm-rv/r
* Determine the angle between the new and old orbital angular
* momentum vectors.
      cmu = (vr*salpha-vk*ctheta*cphi)/SQRT(v1 + v2)
      mu = ACOS(cmu)
* Calculate the components of the velocity of the new centre-of-mass.
 90   continue
      if(ecc.le.1.0)then
* Calculate the components of the velocity of the new centre-of-mass.
         mx1 = vk*m1n/(m1n+m2)
         mx2 = vr*(m1-m1n)*m2/((m1n+m2)*(m1+m2))
         vs(1) = mx1*ctheta*cphi + mx2*salpha
         vs(2) = mx1*stheta*cphi + mx2*calpha
         vs(3) = mx1*sphi
*
* J.S.: Want to have vs in the frame of y
         vs = vs(1)*ux + vs(2)*uy + vs(3)*uz
      else
* Calculate the relative hyperbolic velocity at infinity (simple method).
         sep = r/(ecc-1.d0)
*        cmu = SQRT(ecc-1.d0)
*        mu = ATAN(cmu)
         mu = ACOS(1.d0/ecc)
         vr2 = gmrkm*(m1n+m2)/sep
         vr = SQRT(vr2)
         vs(1) = vr*SIN(mu)
         vs(2) = vr*COS(mu)
         vs(3) = 0.d0
         ecc = MIN(ecc,99.99d0)
      endif
*
 95   continue
* Save some information about the colculation of the natal kick velocity
c      if(kw.eq.13)then
c         WRITE(33,*)vk,phi,theta,ecc,sep,mfin,mrem,m1n
c      elseif(kw.eq.14)then
c         WRITE(44,*)vk,phi,theta,ecc,sep,mfin,mrem,m1n
c      endif
*
      RETURN
      END
***
      FUNCTION dotp(a, b)
      real*8 dotp
      real*8, DIMENSION(3), INTENT(IN) :: a, b
      dotp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      END FUNCTION dotp