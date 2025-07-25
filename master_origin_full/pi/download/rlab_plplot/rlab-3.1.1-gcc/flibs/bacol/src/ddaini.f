      SUBROUTINE DDAINI (X, Y, YPRIME, NEQ, RES, JAC, H, WT, IDID, RPAR,
     +   IPAR, PHI, DELTA, E, WM, IWM, HMIN, UROUND, NONNEG, NTEMP)
C***BEGIN PROLOGUE  DDAINI
C***SUBSIDIARY
C***PURPOSE  Initialization routine for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDAINI-S, DDAINI-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
c-----------------------------------------------------------------------
c     This routine has been modified for the purpose of scaling the
c     interation matrix.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, September 4, 2001.
c
c-----------------------------------------------------------------------
C***DESCRIPTION
C-----------------------------------------------------------------
C     DDAINI TAKES ONE STEP OF SIZE H OR SMALLER
C     WITH THE BACKWARD EULER METHOD, TO
C     FIND YPRIME.  X AND Y ARE UPDATED TO BE CONSISTENT WITH THE
C     NEW STEP.  A MODIFIED DAMPED NEWTON ITERATION IS USED TO
C     SOLVE THE CORRECTOR ITERATION.
C
C     THE INITIAL GUESS FOR YPRIME IS USED IN THE
C     PREDICTION, AND IN FORMING THE ITERATION
C     MATRIX, BUT IS NOT INVOLVED IN THE
C     ERROR TEST. THIS MAY HAVE TROUBLE
C     CONVERGING IF THE INITIAL GUESS IS NO
C     GOOD, OR IF G(X,Y,YPRIME) DEPENDS
C     NONLINEARLY ON YPRIME.
C
C     THE PARAMETERS REPRESENT:
C     X --         INDEPENDENT VARIABLE
C     Y --         SOLUTION VECTOR AT X
C     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
C     NEQ --       NUMBER OF EQUATIONS
C     H --         STEPSIZE. IMDER MAY USE A STEPSIZE
C                  SMALLER THAN H.
C     WT --        VECTOR OF WEIGHTS FOR ERROR
C                  CRITERION
C     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS
C                  IDID= 1 -- YPRIME WAS FOUND SUCCESSFULLY
C                  IDID=-12 -- DDAINI FAILED TO FIND YPRIME
C     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS
C                  THAT ARE NOT ALTERED BY DDAINI
C     PHI --       WORK SPACE FOR DDAINI
C     DELTA,E --   WORK SPACE FOR DDAINI
C     WM,IWM --    REAL AND INTEGER ARRAYS STORING
C                  MATRIX INFORMATION
C
C-----------------------------------------------------------------
C***ROUTINES CALLED  DDAJAC, DDANRM, DDASLV
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C   901030  Minor corrections to declarations.  (FNF)
C***END PROLOGUE  DDAINI
C
      INTEGER  NEQ, IDID, IPAR(*), IWM(*), NONNEG, NTEMP
      DOUBLE PRECISION
     *   X, Y(*), YPRIME(*), H, WT(*), RPAR(*), PHI(NEQ,*), DELTA(*),
     *   E(*), WM(*), HMIN, UROUND
      EXTERNAL  RES, JAC
C
      EXTERNAL  DDAJAC, DDANRM, DDASLV
      DOUBLE PRECISION  DDANRM
C
      INTEGER  I, IER, IRES, JCALC, LNJE, LNRE, M, MAXIT, MJAC, NCF,
     *   NEF, NSF
      DOUBLE PRECISION
     *   CJ, DAMP, DELNRM, ERR, OLDNRM, R, RATE, S, XOLD, YNORM
      LOGICAL  CONVGD
C
      PARAMETER (LNRE=12)
      PARAMETER (LNJE=13)
C
      DATA MAXIT/10/,MJAC/5/
      DATA DAMP/0.75D0/
C
C
C---------------------------------------------------
C     BLOCK 1.
C     INITIALIZATIONS.
C---------------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT  DDAINI
      IDID=1
      NEF=0
      NCF=0
      NSF=0
      XOLD=X
      YNORM=DDANRM(NEQ,Y,WT,RPAR,IPAR)
C
C     SAVE Y AND YPRIME IN PHI
      DO 100 I=1,NEQ
         PHI(I,1)=Y(I)
100      PHI(I,2)=YPRIME(I)
C
C
C----------------------------------------------------
C     BLOCK 2.
C     DO ONE BACKWARD EULER STEP.
C----------------------------------------------------
C
C     SET UP FOR START OF CORRECTOR ITERATION
200   CJ=1.0D0/H
      X=X+H
C
C     PREDICT SOLUTION AND DERIVATIVE
      DO 250 I=1,NEQ
250     Y(I)=Y(I)+H*YPRIME(I)
C
      JCALC=-1
      M=0
      CONVGD=.TRUE.
C
C
C     CORRECTOR LOOP.
300   IWM(LNRE)=IWM(LNRE)+1
      IRES=0
C
      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES.LT.0) GO TO 430
C
C
C     EVALUATE THE ITERATION MATRIX
      IF (JCALC.NE.-1) GO TO 310
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      CALL DDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,
     *   IER,WT,E,WM,IWM,RES,IRES,
     *   UROUND,JAC,RPAR,IPAR,NTEMP)
C
      S=1000000.D0
      IF (IRES.LT.0) GO TO 430
      IF (IER.NE.0) GO TO 430
      NSF=0
C
C
C
C     MULTIPLY RESIDUAL BY DAMPING FACTOR
310   CONTINUE
      DO 320 I=1,NEQ
320      DELTA(I)=DELTA(I)*DAMP
C
C     COMPUTE A NEW ITERATE (BACK SUBSTITUTION)
C     STORE THE CORRECTION IN DELTA
C
c-----------------------------------------------------------------------
      CALL DDASLV(NEQ,DELTA,WM,IWM,CJ)
c-----------------------------------------------------------------------
C
C     UPDATE Y AND YPRIME
      DO 330 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
330      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C
C     TEST FOR CONVERGENCE OF THE ITERATION.
C
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.LE.100.D0*UROUND*YNORM)
     *   GO TO 400
C
      IF (M.GT.0) GO TO 340
         OLDNRM=DELNRM
         GO TO 350
C
340   RATE=(DELNRM/OLDNRM)**(1.0D0/M)
      IF (RATE.GT.0.90D0) GO TO 430
      S=RATE/(1.0D0-RATE)
C
350   IF (S*DELNRM .LE. 0.33D0) GO TO 400
C
C
C     THE CORRECTOR HAS NOT YET CONVERGED. UPDATE
C     M AND AND TEST WHETHER THE MAXIMUM
C     NUMBER OF ITERATIONS HAVE BEEN TRIED.
C     EVERY MJAC ITERATIONS, GET A NEW
C     ITERATION MATRIX.
C
      M=M+1
      IF (M.GE.MAXIT) GO TO 430
C
      IF ((M/MJAC)*MJAC.EQ.M) JCALC=-1
      GO TO 300
C
C
C     THE ITERATION HAS CONVERGED.
C     CHECK NONNEGATIVITY CONSTRAINTS
400   IF (NONNEG.EQ.0) GO TO 450
      DO 410 I=1,NEQ
410      DELTA(I)=MIN(Y(I),0.0D0)
C
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.GT.0.33D0) GO TO 430
C
      DO 420 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
420      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
      GO TO 450
C
C
C     EXITS FROM CORRECTOR LOOP.
430   CONVGD=.FALSE.
450   IF (.NOT.CONVGD) GO TO 600
C
C
C
C-----------------------------------------------------
C     BLOCK 3.
C     THE CORRECTOR ITERATION CONVERGED.
C     DO ERROR TEST.
C-----------------------------------------------------
C
      DO 510 I=1,NEQ
510      E(I)=Y(I)-PHI(I,1)
      ERR=DDANRM(NEQ,E,WT,RPAR,IPAR)
C
      IF (ERR.LE.1.0D0) RETURN
C
C
C
C--------------------------------------------------------
C     BLOCK 4.
C     THE BACKWARD EULER STEP FAILED. RESTORE X, Y
C     AND YPRIME TO THEIR ORIGINAL VALUES.
C     REDUCE STEPSIZE AND TRY AGAIN, IF
C     POSSIBLE.
C---------------------------------------------------------
C
600   CONTINUE
      X = XOLD
      DO 610 I=1,NEQ
         Y(I)=PHI(I,1)
610      YPRIME(I)=PHI(I,2)
C
      IF (CONVGD) GO TO 640
      IF (IER.EQ.0) GO TO 620
         NSF=NSF+1
         H=H*0.25D0
         IF (NSF.LT.3.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
620   IF (IRES.GT.-2) GO TO 630
         IDID=-12
         RETURN
630   NCF=NCF+1
      H=H*0.25D0
      IF (NCF.LT.10.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
C
640   NEF=NEF+1
      R=0.90D0/(2.0D0*ERR+0.0001D0)
      R=MAX(0.1D0,MIN(0.5D0,R))
      H=H*R
      IF (ABS(H).GE.HMIN.AND.NEF.LT.10) GO TO 690
         IDID=-12
         RETURN
690      GO TO 200
C
C-------------END OF SUBROUTINE DDAINI----------------------
      END
