      MODULE commonV
      real*8 c,pi,ut,Rsun
      real*8 Rmax1,Rmax2,Rmax3
      real*8 m12,m123
      real*8 K1,K2,Kq1,Kq2,KL1,KL2
      real*8 m10,m20,m30,R10,R20
      real*8 y0(20)
      real*8, parameter:: LK = 1.d0
      real*8, parameter:: GR = 1.d0
      real*8, parameter:: GW = 1.d0
      real*8, parameter:: Sitter = 1.d0
      real*8, parameter:: Thirring = 1.d0
      real*8, parameter:: Rot = 1.d0
      real*8, parameter:: Tides = 1.d0
      real*8, parameter:: Mass = 1.d0
      real*8 mco1,mco2,mco3,tco1,tco2,tco3
      integer kwfin1,kwfin2,kwfin3,js1,js2,js3
      CHARACTER(len=32) :: arg1,arg2,arg3,arg4,arg5,arg6,arg7,directory
      SAVE c,pi,ut,Rsun
      SAVE K1,K2,Kq1,Kq2,KL1,KL2
      SAVE Rmax1,Rmax2,Rmax3
      SAVE m10,m20,m30,R10,R20
      SAVE kwfin1,kwfin2,kwfin3
      SAVE y0
      SAVE arg1,arg2,arg3,arg4,arg5,arg6,arg7,directory
      END MODULE commonV
     
      MODULE tables
      integer, parameter:: Nmax=1000000
      integer cont1,cont2,cont3
      real*8,dimension(Nmax):: tabt1,tabm1,tabr1,tabdm1
      real*8,dimension(Nmax):: tabt2,tabm2,tabr2,tabdm2
      real*8,dimension(Nmax):: tabt3,tabm3,tabr3,tabdm3      
      real*8,dimension(Nmax):: p1,p2,p3,p4,p5,p6,p7,p8,p9
      real*8,dimension(Nmax):: d1,d2,d3,d4,d5,d6,d7,d8,d9
      real*8,dimension(Nmax):: t1,t2,t3,t4,t5,t6,t7,t8,t9
      real*8,dimension(Nmax):: tabmenv1,tabrenv1,tablumin1,kwv1
      real*8,dimension(Nmax):: tabmenv2,tabrenv2,tablumin2,kwv2
      real*8,dimension(Nmax):: tabmenv3,tabrenv3,tablumin3,kwv3
      real*8,dimension(Nmax):: p10,p11,p12,p13,p14,p15,p16,p17,p18
      real*8,dimension(Nmax):: d10,d11,d12,d13,d14,d15,d16,d17,d18
      real*8,dimension(Nmax):: t10,t11,t12,t13,t14,t15,t16,t17,t18
      SAVE tabt1,tabm1,tabr1,tabdm1,p1,p2,p3,p4,p5,p6,p7,p8,p9
      SAVE tabt2,tabm2,tabr2,tabdm2,d1,d2,d3,d4,d5,d6,d7,d8,d9
      SAVE tabt3,tabm3,tabr3,tabdm3,t1,t2,t3,t4,t5,t6,t7,t8,t9
      SAVE tabmenv1,tabrenv1,tablumin1,kwv1
      SAVE tabmenv2,tabrenv2,tablumin2,kwv2
      SAVE tabmenv3,tabrenv3,tablumin3,kwv3
      SAVE p10,p11,p12,p13,p14,p15,p16,p17,p18
      SAVE d10,d11,d12,d13,d14,d15,d16,d17,d18
      SAVE t10,t11,t12,t13,t14,t15,t16,t17,t18
      SAVE cont1,cont2,cont3
      END MODULE tables

