      SUBROUTINE BRAKE4(I1,I2,KCASE,DT)
*
*
*       GR analytical orbit shrinkage.
*       ------------------------------
*
      INCLUDE 'common6.h'
      COMMON/POSTN/  CVEL,TAUGR,RZ1,GAMMAZ,TKOZ,EMAX,TSP,KZ24,IGR,IPN
      COMMON/KSPAR/  ISTAT(KMAX)
      REAL*8 M1,M2,UI(4),UIDOT(4)
      SAVE ITER
      DATA ITER /0/
*
*
*       Check relativistic conditions (at least one >= NS).
      IF (MAX(KSTAR(I1),KSTAR(I2)).LT.13) GO TO 100
*
*       Specify the basic elements from BH or KS treatments.
      M1 = BODY(I1)
      M2 = BODY(I2)
      IPAIR = KVEC(I1)
      I = N + IPAIR
      SEMI = -0.5*BODY(I)/H(IPAIR)
      E2 = (1.0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(E2)
*
*       Replace the BH value of CVEL by CLIGHT for standard binaries.
      IF (CVEL.EQ.0.0D0) CVEL = CLIGHT
*
*       Form da/dt & de/dt according to Peters 1964.
      ADOT =  64.0/5.0*M1*M2*(M1+M2)/(CVEL**5*SEMI**3*(1.0 - E2)**3.5)
      ADOT = ADOT*(1.0 + 73.0/24.0*E2 + 37.0/96.0*E2**2)
      EDOT = 304.0/15.0*ECC*M1*M2*(M1 + M2)/(CVEL**5*SEMI**4)
      EDOT = EDOT/(1.0 - E2)**2.5*(1.0 + 121.0/304.0*E2)
      TGR = SEMI/ADOT
*
*       Skip for long time-scale (low probability of unperturbed orbit).
      IF (TGR.GT.500.0) GO TO 100
*
*       Set local KS vectors.
      DO 2 K = 1,4
          UI(K) = U0(K,IPAIR)
          UIDOT(K) = UDOT(K,IPAIR)
    2 CONTINUE
*
*       Evaluate the Einstein shift per orbit and check time-step.
      DW = 3.0*TWOPI*(BODY(I1) + BODY(I2))/(SEMI*CVEL**2*(1.0-E2))
      TK = TWOPI*SEMI*SQRT(SEMI/(BODY(I1) + BODY(I2)))
*       Adopt time-scale of 2 % relative change (subject to c.m. step).
      DT = MIN(0.02*SEMI/ADOT,STEP(I))
      THETA = DW*DT/TK
*       Impose limit of time-step if THETA > TWOPI.
      IF (THETA.GT.TWOPI) THEN
          DT = TWOPI*TK/DW
      END IF
      STEP(I1) = DT
*     THETA = DMOD(THETA,TWOPI)   ! original condition.
*
*       Rotate KS orbit by THETA/2 (period halving).
      CALL KSROT(UI,UIDOT,THETA)
*
*       Copy back to KS common variables.
      DO 15 K = 1,4
          U(K,IPAIR) = UI(K)
          U0(K,IPAIR) = UI(K)
          UDOT(K,IPAIR) = UIDOT(K)
   15 CONTINUE
*
*       Define RZ and obtain the incremental change in SEMI & ECC.
      RZ = 8.0*(M1 + M2)/CVEL**2
      SEMI1 = MAX(SEMI - ADOT*DT,RZ)
      ECC1 = MAX(ECC - EDOT*DT,0.0D0)
      IF (SEMI1.LT.0.5*SEMI) THEN
          WRITE (6,20)  NAME(I1), SEMI, SEMI1, ADOT*DT
   20     FORMAT (' PN WARNING    NM A A1 ADOT*DT ',I6,1P,3E10.2)
      END IF
*
*       Update binding energy and collision energy.
      HI = H(IPAIR)
      H(IPAIR) = -0.5*BODY(I)/SEMI1
      ZMU = BODY(I1)*BODY(I2)/BODY(I) 
      ECOLL = ECOLL + ZMU*(HI - H(IPAIR))
*
      IF (ITER.LT.10000) THEN
          WRITE (93,30)  NAME(I1), TIME+TOFF, ECC1, SEMI1
   30     FORMAT (' ',I7,F10.3,F9.5,1P,E12.4)
          CALL FLUSH(93)
      END IF
*
*       Change KS variables at original ecc and modify ECC at H = const.
      CALL EXPAND2(IPAIR,SEMI)
      CALL KSPERI(IPAIR)
*       Note simplified versions of standard routines (new STEP(I1) given).
      CALL DEFORM2(IPAIR,ECC,ECC1)
*
      ITER = ITER + 1
      IF (ITER.LT.1000.OR.MOD(ITER,1000).EQ.0) THEN
          WRITE (94,40)  TIME+TOFF, ECC, THETA, DT, TGR, SEMI
   40     FORMAT (' GR SHRINK    T E TH DT TGR A ',
     &                           F11.4,F9.5,1P,3E9.1,E12.4)
          CALL FLUSH(94)
      END IF
*
*       Check KS termination with added perturber to activate PN.
      IF (TGR.LT.0.5) THEN
*       Note thst first order Peters formulation is not valid for strong GR.
          JP = LIST(2,I)
          LIST(1,I1) = 1
          LIST(2,I1) = JP
*       Set PN indicator for ARCHAIN (TGR limit means small TZ).
          IPN = 3
          WRITE (6,44)  JP, NAME(JP), STEP(I1), STEP(I), SEMI, TGR
   44     FORMAT (' ENFORCED PERTURB    JP NM S1 SI A TZ',2I6,1P,5E10.2)
          GO TO 100
      END IF
*
*       Activate coalescence condition using local index.
      JPHASE = 0
      IF (SEMI1.LT.100.0*RZ.AND.TGR.LT.0.1) THEN
          WRITE (6,45)  KSTAR(I1), KSTAR(I2), RADIUS(I1), RADIUS(I2),
     &                  SEMI1
   45     FORMAT (' PN COAL    K* R1 R2 A  ',2I4,1P,3E10.2)
          CALL FLUSH(6)
          IQCOLL = -2
          JPHASE = -1
*         KSPAIR = IPAIR
*         CALL CMBODY(SEMI1,2)
*
*       Include optional kick velocity of 3*VRMS km/s for coalescence recoil.
          IF (KZ(43).GT.0) THEN
*       Initialize basic variables at start of new step.
              VI20 = 0.0
              DO 48 K = 1,3
                  X0(K,I) = X(K,I)
                  X0DOT(K,I) = XDOT(K,I)
                  VI20 = VI20 + XDOT(K,I)**2
   48         CONTINUE
*
              VF = 3.0*(VRMS/VSTAR)/SQRT(VI20)
              DO 50 K = 1,3
                  XDOT(K,I) = VF*XDOT(K,I)
                  X0DOT(K,I) = XDOT(K,I)
   50         CONTINUE
              ECD0 = ECDOT
              ECDOT = ECDOT + 0.5*BODY(I)*VI20*(1.0 - VF**2)
              VESC = 3.0*VRMS
              WRITE (6,60)  VF, ECD0-ECDOT, VESC
   60         FORMAT (' COALESCENCE KICK    VF ECDOT VESC ',
     &                                      F7.3,F10.6,F6.1)
*       Form neighbour list and new polynomials.
              RS0 = RS(I)
              CALL NBLIST(I,RS0)
              CALL FPOLY1(I,I,0)
              CALL FPOLY2(I,I,0)
          END IF
      END IF
*
      IF (ISTAT(KCASE).EQ.0) ISTAT(KCASE) = JPHASE
*
  100 RETURN
*
      END
