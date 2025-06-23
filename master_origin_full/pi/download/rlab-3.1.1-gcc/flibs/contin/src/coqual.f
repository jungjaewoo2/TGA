      SUBROUTINE COQUAL(STEPX,IWORK,LIW,RWORK,LRW)
C
C***********************************************************************
C
C  COQUAL computes the factor QUAL which is based on the 'quality'
C  (number of steps, last step/total step) of the corrector iteration.  
C  This factor is used as part of the stepsize algorithm.  
C
C  See the paper "On Steplength Algorithms for a Class of Continuation
C  Methods", listed in the documentation.
C
C  STEPX  = Input, size of the last correction.
C  CORDIS = Input, distance between first and last point.
C  MAXCOR = Input, maximum number of corrections allowed.
C  MODNEW = Input, defines type of corrector process used.
C  NCOR   = Input, actual number of corrections used.
C  QUAL   = Output, the quality factor, between 1/8 and 8,
C           with 1/8 poor, 1 average, 8 superior.  Used to
C           modify the next stepsize.
C
      DOUBLE PRECISION EIGHT
      DOUBLE PRECISION ONE
      DOUBLE PRECISION THGIE
C
      PARAMETER (EIGHT=8.0)
      PARAMETER (ONE=1.0)
      PARAMETER (THGIE=0.125)
C
      INTEGER   LIW
      INTEGER   LRW
C
      DOUBLE PRECISION BASE
      DOUBLE PRECISION BOT
      DOUBLE PRECISION CORDIS
      DOUBLE PRECISION EPSATE
      DOUBLE PRECISION ESAB
      DOUBLE PRECISION EXPO
      INTEGER   IBOT
      INTEGER   ITOP
      INTEGER   IWORK(LIW)
      INTEGER   LOUNIT
      INTEGER   MAXCOR
      INTEGER   MODNEW
      INTEGER   NAVE
      INTEGER   NMAX
      DOUBLE PRECISION QUAL
      DOUBLE PRECISION RWORK(LRW)
      DOUBLE PRECISION STEPX
      DOUBLE PRECISION TERM
      DOUBLE PRECISION TEST
      DOUBLE PRECISION TOP
C
      CORDIS=RWORK(15)
      EPSATE=EIGHT*RWORK(8)
      MODNEW=IWORK(4)
      MAXCOR=IWORK(17)
C
C  SEE IF MINIMAL NUMBER OF STEPS TAKEN
C
      IF((IWORK(28).LE.1).OR.(CORDIS.LE.EPSATE)) THEN
        QUAL=EIGHT
        GO TO 88
        ENDIF
C
C  SEE IF AVERAGE NUMBER OF STEPS TAKEN
C
      IF(MODNEW.EQ.0)THEN
        NAVE=(MAXCOR-1)/2
      ELSE
        NAVE=MAXCOR
        ENDIF
      IF(IWORK(28).EQ.NAVE)THEN
        QUAL=ONE
        GO TO 88
        ENDIF
C
C  SEE IF MAXIMUM NUMBER OF STEPS TAKEN
C
      IF(MODNEW.EQ.0)THEN
        NMAX=MAXCOR
      ELSE
        NMAX=2*MAXCOR
        ENDIF
      IF(IWORK(28).GE.NMAX)THEN
        QUAL=THGIE
        GO TO 88
        ENDIF
C
C  COMPUTE QUAL
C
C  FOR MODNEW=0, LET W=(STEPX/CORDIS),
C                    IEXP=1/(2**(NCOR-1)-1)
C                    U=W**IEXP
C                    JEXP=2**(NCOR-NAVE)
C              THEN  QUAL=(U+1+(1/U)) / (U**JEXP+1+(1/U**JEXP))
C
C  OTHERWISE,    EXP=(NCOR-NAVE)/(NCOR-1)
C                W=(STEPX/CORDIS)
C                QUAL=W**EXP
C
      IF(MODNEW.EQ.0)THEN
        ITOP=2**(IWORK(28)-NAVE)
        IBOT=2**(IWORK(28)-1)-1
        EXPO=ONE/FLOAT(IBOT)
        BASE=(STEPX/CORDIS)**EXPO
        TOP=BASE+ONE+(ONE/BASE)
        TERM=BASE**ITOP
        BOT=TERM+ONE+(ONE/TERM)
        QUAL=TOP/BOT
        IF(QUAL.GT.EIGHT)QUAL=EIGHT
        IF(QUAL.LT.THGIE)QUAL=THGIE
        GO TO 88
        ENDIF
C
C  TRY TO ANTICIPATE UNDERFLOW AND OVERFLOW
C
      ITOP=IWORK(28)-NAVE
      IBOT=IWORK(28)-1
      EXPO=FLOAT(IBOT)/FLOAT(ITOP)
      TEST=EIGHT**EXPO
      BASE=STEPX/CORDIS
      ESAB=CORDIS/STEPX
      IF((IWORK(28).LT.NAVE.AND.TEST.GT.BASE)
     1  .OR.(IWORK(28).GT.NAVE.AND.TEST.LT.BASE))THEN
        QUAL=EIGHT
        GO TO 88
        ENDIF
      IF((IWORK(28).LT.NAVE.AND.TEST.GT.ESAB)
     1  .OR.(IWORK(28).GT.NAVE.AND.TEST.LT.ESAB))THEN
        QUAL=THGIE
        GO TO 88
        ENDIF
C
C  COMPUTE QUAL
C
      EXPO=FLOAT(ITOP)/FLOAT(IBOT)
      QUAL=BASE**EXPO
      IF(QUAL.GT.EIGHT)QUAL=EIGHT
      IF(QUAL.LT.THGIE)QUAL=THGIE
88    CONTINUE
      RWORK(23)=QUAL
      IF(IWORK(7).GE.3)THEN
        LOUNIT=IWORK(8)
        WRITE(LOUNIT,1000)QUAL
        ENDIF
      RETURN
1000  FORMAT(' COQUAL - Corrector convergence quality factor=',G14.6)
      END
