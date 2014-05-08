      SUBROUTINE FNFW(XI,XIDOT,FM,FD)
*
*
*       NFW halo force.
*       ---------------
*
      IMPLICIT REAL*8  (A-H,O-Z)
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
      REAL*8  XI(3),XIDOT(3),FM(3),FD(3)
*
*
*       Obtain force and first derivative for NFW potential.
      RTEMP = SQRT(XI(1)**2 + XI(2)**2 + (XI(3)/ZDUM(4))**2)
      RDTEMP = (XI(1)*XIDOT(1) + XI(2)*XIDOT(2) +
     &           XI(3)/ZDUM(4)*XIDOT(3))/RTEMP
      H2 = ZDUM(2)/(RTEMP**2)
      BLOG = LOG(1.0D0+RTEMP/ZDUM(3))
      CTEMP = RTEMP*ZDUM(3)+RTEMP**2
      CNFW = BLOG/RTEMP - 1.0/(ZDUM(3)+RTEMP)
      BNFW = (ZDUM(3)*BLOG + RTEMP*(BLOG-1.0))/CTEMP
      ANFW = RDTEMP*(-3.0*BLOG*(ZDUM(3)**2) + (4.0-3.0*BLOG)
     &   *(RTEMP**2) + RTEMP*ZDUM(3)*(3.0-6.0*BLOG))/(CTEMP**2)
*
      DO 10 K = 1,2
          FM(K) = -H2*CNFW*XI(K)
          FD(K) = -H2*(XIDOT(K)*BNFW + XI(K)*ANFW)
   10 CONTINUE

*       Add flattening along z-axis
      FM(3) = -H2*CNFW*XI(3)/(ZDUM(4)**2)
      FD(3) = -H2*(XIDOT(3)*BNFW + XI(3)*ANFW)/(ZDUM(4)**2)
*
      RETURN
*
      END
