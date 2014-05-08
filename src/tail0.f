      SUBROUTINE TAIL0(I)
*
*
*      Initialization of tidal tail member.
*      ------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/GALAXY/ GMG,RG(3),VG(3),FG(3),FGD(3),TG,
     &               OMEGA,DISK,A,B,V02,RL2,GMB,AR,GAM,ZDUM(7)
      REAL*8  XI(3),XIDOT(3),FIRR(3),FD(3)
      SAVE IWARN
      DATA IWARN /0/
*
*
*       Transform to galactic coordinates.
      DO 5 K = 1,3
          XI(K) = X(K,I) + RG(K)
          XIDOT(K) = XDOT(K,I) + VG(K)
          FIRR(K) = 0.0
          FD(K) = 0.0
    5 CONTINUE
*
*       Obtain relevant force components.
      CALL XTRNLT(XI,XIDOT,FIRR,FD)
*
*       Update membership (first member in NTTOT = ITAIL0).
      IF (NTAIL.EQ.0) THEN
          NSTAIL = 0
*       Allow extra space for ten KS solutions.
          ITAIL0 = NZERO + MIN(KMAX,NBIN0+10)
          NTTOT = ITAIL0 - 1
      END IF
*
*       Include safety test and 10 warnings for maximum membership.
      IF (NTTOT.GE.NMAX) THEN
          IWARN = IWARN + 1
          IF (IWARN.LE.10) THEN
              WRITE (6,10)  NTTOT
   10         FORMAT (' WARNING!    MAXIMUM TIDAL TAIL ',I6)
          END IF
          RETURN
      END IF
*
*       Increase the memberships (note initial value ITAIL0 - 1).
      NTTOT = NTTOT + 1
      NTAIL = NTAIL + 1
      J = NTTOT
*
*       Copy escaper data and initialize integration variables.
      NAME(J) = NAME(I)
      BODY(J) = BODY(I)
      FF = 0.0
      FFD = 0.0
      DO 20 K = 1,3
          X(K,J) = XI(K)
          X0(K,J) = XI(K)
          XDOT(K,J) = XIDOT(K)
          X0DOT(K,J) = XIDOT(K)
          F(K,J) = 0.5*FIRR(K)
          FDOT(K,J) = ONE6*FD(K)
          FI(K,J) = FIRR(K)
          FIDOT(K,J) = FD(K)
          FF = FF + FIRR(K)**2
          FFD = FFD + FD(K)**2
   20 CONTINUE
*
*       Form quantized time-step (saves prediction at nice output times).
      DT = 0.5*ETAI*SQRT(FF/FFD)
      CALL STEPK(DT,DTN)
      STEP(J) = DTN
      T0(J) = TIME
      TNEW(J) = T0(J) + STEP(J)
*
      RETURN
*
      END
