      SUBROUTINE START(DF,FPAR,FX,IERROR,IPAR,IPC,IWORK,IWRITE,LIW,
     *LOUNIT,LRW,NVAR,RWORK,TC,WORK1,XC,XF,XR,SLNAME)
C
C***********************************************************************
C
C  On first call, make sure that the input point satisfies the
C  nonlinear equations.  If not, call CORECT and attempt to force this.
C  
      DOUBLE PRECISION ONE
      DOUBLE PRECISION ZERO
C
      PARAMETER (ONE=1.0)
      PARAMETER (ZERO=0.0)
C
      EXTERNAL  COQUAL
      EXTERNAL  CORECT
      EXTERNAL  DF
      EXTERNAL  FX
      EXTERNAL  IDAMAX
      EXTERNAL  DAXPY
      EXTERNAL  DCOPY
      EXTERNAL  SLNAME
C
      INTRINSIC ABS
C
      INTEGER   LIW
      INTEGER   LRW
      INTEGER   NVAR
C
      DOUBLE PRECISION DETS
      DOUBLE PRECISION FPAR(*)
      INTEGER   I
      INTEGER   ICRIT
      INTEGER   IERROR
      INTEGER   IMAX
      INTEGER   IPAR(*)
      INTEGER   IPC
      INTEGER   IDAMAX
      INTEGER   IWORK(LIW)
      INTEGER   IWRITE
      INTEGER   JOB
      INTEGER   LOUNIT
      INTEGER   MODSAV
      DOUBLE PRECISION RWORK(LRW)
      DOUBLE PRECISION SKALE
      DOUBLE PRECISION STEPX
      DOUBLE PRECISION TC(NVAR)
      DOUBLE PRECISION WORK1(NVAR)
      DOUBLE PRECISION XC(NVAR)
      DOUBLE PRECISION XF(NVAR)
      DOUBLE PRECISION XR(NVAR)
C
C  If user is requesting that Jacobian be used as long as possible,
C  then go ahead, generate and factor the first one now.
C
      IF(IWORK(4).EQ.2)THEN
        JOB=2
        CALL SLNAME(DETS,FX,DF,FPAR,IERROR,IPC,IPAR,IWORK,LIW,
     1  JOB,NVAR,RWORK,LRW,XR,WORK1)
        RWORK(17)=DETS
        IF(IERROR.NE.0)THEN
          WRITE(LOUNIT,*)
     *    'START  - Could not factor initial jacobian.'
          RETURN
          ENDIF
        ENDIF
      IF(IWRITE.GE.2)WRITE(LOUNIT,1040)IPC
1040  FORMAT(' START  - Correct initial point, fixing index ',I5)
C
C  Set up a pseudo-tangent vector and other data for first time,
C  check that starting point satisfies the equations.
C
      DO 50 I=1,NVAR
        TC(I)=ZERO
50      CONTINUE
      CALL DCOPY(NVAR,XR,1,XC,1)
      TC(IPC)=ONE
      MODSAV=IWORK(4)
      ICRIT=1
C
60    CONTINUE
      CALL DCOPY(NVAR,XC,1,XR,1)
      CALL CORECT(DF,FPAR,FX,IERROR,IPC,IPAR,IWORK,NVAR,RWORK,STEPX,
     *WORK1,XR,LRW,LIW,ICRIT,SLNAME)
      IWORK(25)=IWORK(25)+IWORK(28)
C
C  If possible, retry with ICRIT=2.
C
      IF(IERROR.NE.0.AND.ICRIT.EQ.1)THEN
        IF(IWRITE.GE.1)WRITE(LOUNIT,*)
     *  'START -  Retry starting point correction'
        ICRIT=2
        GO TO 60
      ELSEIF(IERROR.NE.0)THEN
        ICRIT=1
        ENDIF
C
      IF(IERROR.NE.0.AND.IWORK(4).GT.0)THEN
        IWORK(4)=IWORK(4)-1
        IERROR=0
1110  FORMAT(' START  - Retrying starting point with IWORK(4)=',I6)
        IF(IWRITE.GE.1)WRITE(LOUNIT,1110)IWORK(4)
        GO TO 60
        ENDIF
      IWORK(4)=MODSAV
      IF(IERROR.NE.0)THEN
        WRITE(LOUNIT,*)
     *  'START  - Starting point correction failed.'
        RETURN
        ENDIF
      SKALE=-ONE
      CALL DAXPY(NVAR,SKALE,XR,1,XC,1)
      IMAX=IDAMAX(NVAR,XC,1)
      RWORK(15)=ABS(XC(IMAX))
      CALL DCOPY(NVAR,XR,1,XC,1)
      CALL DCOPY(NVAR,XR,1,XF,1)
      CALL COQUAL(STEPX,IWORK,LIW,RWORK,LRW)
      RWORK(14)=RWORK(13)
      IWORK(27)=IWORK(27)+1
      IWORK(10)=1
      IWORK(1)=1
      RETURN
      END 
