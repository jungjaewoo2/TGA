       SUBROUTINE CORECT(DF,FPAR,FX,IERROR,IHOLD,IPAR,IWORK,
     1 NVAR,RWORK,STEPX,WK,XR,LRW,LIW,ICRIT,SLNAME)
C
C***********************************************************************
C
C  CORECT PERFORMS THE CORRECTION OF AN APPROXIMATE
C  SOLUTION OF THE EQUATION F(XR)=0.0. FOR THE CORRECTION  EITHER
C  NEWTON'S METHOD OR THE CHORD NEWTON METHOD IS USED. IN THE LATTER
C  CASE THE JACOBIAN IS  EVALUATED ONLY AT THE STARTING POINT.
C  IF B IS THE VALUE OF XR(IHOLD) FOR THE INPUT STARTING POINT,
C  THEN THE AUGMENTING EQUATION IS XR(IHOLD)=B, THAT IS, THE IHOLD-TH
C  COMPONENT OF XR IS TO BE HELD FIXED.  THE AUGMENTED SYSTEM TO BE
C  SOLVED IS THEN DFA(XR,IHOLD)*(-DELX)=FA(XR)
C
C  ACCEPTANCE AND REJECTION CRITERIA-
C
C  LET NCOR DENOTE THE CURRENT ITERATION INDEX; THAT IS,  NCOR=0 FOR
C  THE PREDICTED POINT WHICH SERVES AS STARTING POINT.
C  After each iteration, the current point will be accepted
C  if any of the following conditions hold:
C
C  STRONG ACCEPTANCE CRITERION
C
C  1.  FNRM.LE.ABSERR AND STEPX.LE.(ABSERR+RELERR*XNRM))
C
C  WEAK ACCEPTANCE CRITERIA
C
C  2.  (NCOR.EQ.0) AND
C      FNRM.LE.(.5*ABSERR)
C  3.  FNRM.LE.EPSATE OR STEPX.LE.EPSATE
C  4.  (NCOR.GE.2) AND
C      (FNRM+FNRML).LE.ABSERR AND
C      STEPX.LE.8*(ABSERR+RELERR*XNRM)
C  5.  (NCOR.GE.2) AND
C      FNRM.LE.8.0*ABSERR AND
C      (STEPX+STEPXL).LE.(ABSERR+RELERR*XNRM)
C
C  AFTER EACH STEP OF THE ITERATION, THE PROCESS IS TO BE ABORTED IF
C  ANY OF THE FOLLOWING CONDITIONS HOLDS
C
C  1.  FNRM.GT.(FMP*FNRML+ABSERR)
C  2.  (NCOR.GE.2) AND (ICRIT.EQ.0) AND
C      (STEPX.GT.(FMP*STEPXL+ABSERR))
C  3.  NCOR.GE.MAXCOR (NUMBER OF ITERATIONS EXCEEDED).
C
C  HERE FMP=2.0 FOR NCOR.EQ.1, FMP=1.05 FOR NCOR.GT.1
C
C  ERROR CONDITIONS
C         IERROR=1  DATA OR STORAGE ERROR
C         IERROR=2  ERROR IN FUNCTION OR DERIVATIVE CALL
C         IERROR=3  SOLVER FAILED
C         IERROR=4  UNSUCCESSFUL ITERATION
C         IERROR=5  TOO MANY CORRECTOR STEPS
C
      DOUBLE PRECISION EIGHT
      DOUBLE PRECISION HALF
      DOUBLE PRECISION ONE
      DOUBLE PRECISION ONEFIV
      DOUBLE PRECISION TWO
      DOUBLE PRECISION ZERO
C
      PARAMETER (EIGHT=8.0)
      PARAMETER (HALF=0.5)
      PARAMETER (ONE=1.0)
      PARAMETER (ONEFIV=1.05)
      PARAMETER (TWO=2.0)
      PARAMETER (ZERO=0.0)
C
      EXTERNAL  DF
      EXTERNAL  FX
      EXTERNAL  IDAMAX
      EXTERNAL  DAXPY
      EXTERNAL  SLNAME
      EXTERNAL  DNRM2
C
      INTRINSIC ABS
C
      INTEGER   LIW
      INTEGER   LRW
      INTEGER   NVAR
C
      DOUBLE PRECISION ABSERR
      DOUBLE PRECISION DETS
      DOUBLE PRECISION EPSATE
      DOUBLE PRECISION FMP
      DOUBLE PRECISION FNRM
      DOUBLE PRECISION FNRML
      DOUBLE PRECISION FPAR(*)
      INTEGER   I
      INTEGER   ICRIT
      INTEGER   IERROR
      INTEGER   IFMAX
      INTEGER   IHOLD
      INTEGER   IPAR(*)
      INTEGER   IDAMAX
      INTEGER   IWORK(LIW)
      INTEGER   IWRITE
      INTEGER   IXMAX
      INTEGER   JOB
      INTEGER   KSMAX
      INTEGER   LOUNIT
      INTEGER   MAXCOR
      INTEGER   MAXNEW
      INTEGER   MODNEW
      DOUBLE PRECISION RELERR
      DOUBLE PRECISION RWORK(LRW)
      DOUBLE PRECISION SKALE
      DOUBLE PRECISION DNRM2
      DOUBLE PRECISION STEPX
      DOUBLE PRECISION STEPXL
      DOUBLE PRECISION TLSTEP
      DOUBLE PRECISION WK(NVAR)
      DOUBLE PRECISION XNRM
      DOUBLE PRECISION XR(NVAR)
      DOUBLE PRECISION XVALUE
C
C  Initialize.
C
      ABSERR=RWORK(1)
      RELERR=RWORK(2)
      EPSATE=EIGHT*RWORK(8)
      MODNEW=IWORK(4)
      IWRITE=IWORK(7)
      LOUNIT=IWORK(8)
      MAXCOR=IWORK(17)
      IERROR=0
      IWORK(28)=0
      MAXNEW=MAXCOR
      IF(MODNEW.NE.0)MAXNEW=2*MAXCOR
      FMP=TWO
      STEPX=ZERO
      XVALUE=XR(IHOLD)
C
C  STORE INITIAL FUNCTION VALUE IN THE WORK VECTOR WK
C
      CALL FX(NVAR,FPAR,IPAR,XR,WK,IERROR)
      IWORK(22)=IWORK(22)+1
      IF(IERROR.NE.0)THEN
        WRITE(LOUNIT,*)
     *  'CORECT - Error flag from user function routine.'
        RETURN
        ENDIF
      WK(NVAR)=XR(IHOLD)-XVALUE
      IFMAX=IDAMAX(NVAR,WK,1)
      FNRM=DNRM2(NVAR-1,WK,1)
      IXMAX=IDAMAX(NVAR,XR,1)
      XNRM=DNRM2(NVAR,XR,1)
      IF(IWRITE.GE.2)WRITE(LOUNIT,70)IWORK(28),FNRM,IFMAX
      IF(FNRM.LE.HALF*ABSERR)RETURN
C
C  NEWTON ITERATION LOOP BEGINS
C
      DO 100 I=1,MAXNEW
        IWORK(28)=I
C
C  SOLVE SYSTEM FPRIME*(-DELX)=FX FOR NEWTON DIRECTION DELX
C
        JOB=0
        IF(I.NE.1.AND.I.NE.MAXCOR.AND.MODNEW.EQ.1)JOB=1
        IF(MODNEW.EQ.2)JOB=1
        CALL SLNAME(DETS,FX,DF,FPAR,IERROR,IHOLD,IPAR,IWORK,LIW,
     1  JOB,NVAR,RWORK,LRW,XR,WK)
        IF(IERROR.NE.0)THEN
          WRITE(LOUNIT,12)IERROR
          RETURN
          ENDIF
C
C  UPDATE X, COMPUTE SIZE OF STEP AND OF X
C
        SKALE=-ONE
        CALL DAXPY(NVAR,SKALE,WK,1,XR,1)
        STEPXL=STEPX
        KSMAX=IDAMAX(NVAR,WK,1)
        STEPX=ABS(WK(KSMAX))
        IXMAX=IDAMAX(NVAR,XR,1)
        XNRM=DNRM2(NVAR,XR,1)
C
C  Compute F(X) and take its norm.
C
        CALL FX(NVAR,FPAR,IPAR,XR,WK,IERROR)
        IWORK(22)=IWORK(22)+1
        IF(IERROR.NE.0) THEN
          WRITE(LOUNIT,*)
     *    'CORECT - Error flag from user function routine.'
          RETURN
          ENDIF
        WK(NVAR)=XR(IHOLD)-XVALUE
        FNRML=FNRM
        IFMAX=IDAMAX(NVAR,WK,1)
        FNRM=DNRM2(NVAR-1,WK,1)
        IF(IWRITE.GE.2) THEN
          WRITE(LOUNIT,80)IWORK(28),STEPX,KSMAX
          WRITE(LOUNIT,70)IWORK(28),FNRM,IFMAX
          ENDIF
C
C  Check for strong acceptance of function and stepsize.
C
        TLSTEP=ABSERR+RELERR*XNRM
        IF(FNRM.LE.ABSERR.AND.STEPX.LE.TLSTEP)RETURN
C
C  Check for weak acceptance of function and stepsize.
C
        IF(FNRM.LE.EPSATE.OR.STEPX.LE.EPSATE)RETURN
        IF(IWORK(28).GT.1)THEN
          IF((FNRM+FNRML).LE.ABSERR.AND.STEPX.LE.EIGHT*TLSTEP)RETURN
          IF(FNRM.LE.EIGHT*ABSERR.AND.(STEPX+STEPXL).LE.TLSTEP)RETURN
          ENDIF
C
C  Decide if iteration should be aborted
C
        IF(IWORK(28).GT.1)THEN
          IF(ICRIT.LT.1.AND.STEPX.GT.(FMP*STEPXL+ABSERR)) THEN
            IERROR=4
            IF(IWRITE.GE.2)WRITE(LOUNIT,*)
     *      'CORECT - Size of correction not decreasing.'
            RETURN
            ENDIF
          ENDIF
        IF(ICRIT.LT.2.AND.FNRM.GT.(FMP*FNRML+ABSERR)) THEN
           IERROR=4
           IF(IWRITE.GE.2)WRITE(LOUNIT,*)
     *     'CORECT - Residual is not decreasing.'
           RETURN
           ENDIF
        FMP=ONEFIV
  100   CONTINUE
C
C  Reached maximum number of steps without acceptance or rejection.
C
      IERROR=5
      IF(IWRITE.GE.2)WRITE(LOUNIT,*)'CORECT - Convergence too slow.'
      RETURN
70    FORMAT(' CORECT - Residual ',I6,'=',G14.6,' I=',I6)
80    FORMAT(' CORECT - Step     ',I6,15X,G14.6,' I=',I6)
12    FORMAT(' CORECT - Error flag=',I6,' from solver.')
      END
