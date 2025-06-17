      SUBROUTINE BANJAC(EPS,FCOL,FPAR,FPRIME,FROW,FX,IERROR,IPAR,IPC,
     * IWORK,JAC,LIW,LOUNIT,NBAND,NEQN,NVAR,X,XTEMP,WORK1,WORK2)
C
C**********************************************************************
C
C  BANJAC estimates the jacobian matrix FPRIME of the function FX, 
C  using forward or central finite differences.  BANJAC is called by 
C  BANSLV when the user has specified the jacobian option as 1 or 2.
C
C  EPS    Input, a tolerance used for shifting the X values.  A value
C         of the square root of the machine precision is usually
C         appropriate.
C
C  FCOL   Output, REAL FCOL(NEQN), the last column of the approximate 
C         jacobian, which is allowed to be "full".  This comprises
C         matrix entries FPRIME(1,NVAR) through FPRIME(NEQN,NVAR).
C
C  FPAR   Input/output, REAL FPAR(*), user parameter vector, not
C         touched by this routine, but passed on to user routines.
C
C  FPRIME Output, REAL FPRIME(NBAND,NEQN), is the array into which the
C         the banded portion of the computed jacobian will be stored.
C         The LINPACK general band format is used, assigning entry (I,J)
C         to FPRIME(I-J+ML+MU+1,J), where ML and MU are the lower and
C         upper half bandwidths respectively.
C
C  FROW   Output, REAL FROW(NVAR), storage for the last (augmenting) row
C         of the jacobian, which will be all zero except for a 1 in
C         location IPC.
C
C  FX     Input, EXTERNAL FX, the name of the routine which evaluates the 
C         function.
C
C         FX computes the value of the nonlinear function.  This name must be 
C         declared EXTERNAL in the calling program.  FX should evaluate the 
C         NVAR-1 function components at the input point X, and store the result 
C         in the vector FVEC.  An augmenting equation will be stored in entry 
C         NVAR of FVEC by the PITCON program.
C
C         FX should have the following form:
C
C           SUBROUTINE FX(NVAR,FPAR,IPAR,X,FVEC,IERROR)
C
C           NVAR   Input, INTEGER NVAR, number of variables.
C
C           FPAR   Input/output, REAL FPAR(*), array of user parameters.
C
C           IPAR   Input/output, INTEGER IPAR(*), array of user parameters.
C
C           X      Input, REAL X(NVAR), the point at which function evaluation
C                  is required.
C
C           FVEC   Output, REAL FVEC(NVAR), the value of the function at point 
C                  X.  Only the first NVAR-1 entries of FVEC are to be set by 
C                  the routine.  PITCON sets the final value itself.
C 
C           IERROR Output, INTEGER IERROR, error flag.  FX should set this to 0 
C                  if there are no problems, or to 2 if there is a problem.
C
C  IERROR Output, INTEGER IERROR, error flag.  A nonzero value means that
C         there was an error in the user routine FX, or in BANJAC itself.
C         In either case, the jacobian has not been computed.
C
C  IPAR   Input, INTEGER IPAR(*), a user parameter vector passed to FX.  
C         However, because this is a problem with a banded jacobian, entries
C         IPAR(1) and IPAR(2) are read by this routine.  IPAR(1) contains
C         ML, the lower half bandwidth of the jacobian, and IPAR(2) contains
C         MU, the upper half bandwidth of the jacobian.
C
C  IPC    Input, INTEGER IPC, the index of the current continuation parameter,
C         which is needed to determine the form of FROW.
C
C  IWORK  Input, INTEGER IWORK(LIW), work and statistics vector.  Only
C         required here so that we can count the number of function
C         evaluations.
C
C  JAC    Input, INTEGER JAC, the user requested jacobian option.  For
C         our purposes, the only two values of interest are:
C
C           1 = estimate jacobian with forward differences,
C           2 = estimate jacobian with central differences (twice the work)
C
C  LIW    Input, INTEGER LIW, the dimension of IWORK.
C
C  LOUNIT Input, INTEGER LOUNIT, the FORTRAN output unit for messages.
C
C  NBAND  Input, INTEGER NBAND, the first dimension of the jacobian matrix
C         FPRIME, NBAND=ML+MU+1.
C
C  NEQN   Input, INTEGER NEQN, the number of equations, equal to NVAR-1.
C
C  NVAR   Input, INTEGER NVAR, the number of variables.
C
C  X      Input, REAL X(NVAR), the point at which the jacobian is desired.
C
C  XTEMP,
C  WORK1,
C  WORK2  Work arrays, REAL XTEMP(NVAR), WORK1(NVAR), WORK2(NVAR).
C
      DOUBLE PRECISION ONE
      DOUBLE PRECISION TWO
C
      PARAMETER (ONE=1.0)
      PARAMETER (TWO=2.0)
C
      EXTERNAL  FX
      EXTERNAL  DAXPY
      EXTERNAL  DCOPY
      EXTERNAL  DSCAL
C
      INTRINSIC ABS
      INTRINSIC MAX
      INTRINSIC MIN
C
      INTEGER   LIW
      INTEGER   NBAND
      INTEGER   NEQN
      INTEGER   NVAR
C
      DOUBLE PRECISION EPS
      DOUBLE PRECISION FCOL(NEQN)
      DOUBLE PRECISION FPAR(*)
      DOUBLE PRECISION FPRIME(NBAND,NEQN)
      DOUBLE PRECISION FROW(NVAR)
      INTEGER   IBAND
      INTEGER   IERROR
      INTEGER   IHI
      INTEGER   ILO
      INTEGER   IPAR(2)
      INTEGER   IPC
      INTEGER   IROW
      INTEGER   IWORK(LIW)
      INTEGER   J
      INTEGER   JAC
      INTEGER   KCALL
      INTEGER   LOUNIT
      INTEGER   MBAND
      INTEGER   ML
      INTEGER   MU
      DOUBLE PRECISION SKALE
      DOUBLE PRECISION X(NVAR)
      DOUBLE PRECISION XJAC
      DOUBLE PRECISION XTEMP(NVAR)
      DOUBLE PRECISION WORK1(NVAR)
      DOUBLE PRECISION WORK2(NVAR)
C
      ML=IPAR(1)
      MU=IPAR(2)
      MBAND=ML+MU+1
      IF(JAC.EQ.1)THEN
        CALL FX(NVAR,FPAR,IPAR,X,WORK2,IERROR)
        IWORK(22)=IWORK(22)+1
        IF(IERROR.NE.0)RETURN
        ENDIF
      XJAC=ONE
      IF(JAC.EQ.2)XJAC=TWO
      DO 40 KCALL=1,MBAND
        CALL DCOPY(NVAR,X,1,XTEMP,1)
        DO 10 J=KCALL,NEQN,MBAND
          XTEMP(J)=X(J)+EPS*(ONE+ABS(X(J)))
10        CONTINUE
        CALL FX(NVAR,FPAR,IPAR,XTEMP,WORK1,IERROR)
        IWORK(22)=IWORK(22)+1
        IF(IERROR.NE.0)RETURN
        IF(JAC.EQ.2)THEN
          CALL DCOPY(NVAR,X,1,XTEMP,1)
          DO 20 J=KCALL,NEQN,MBAND
            XTEMP(J)=X(J)-EPS*(ONE+ABS(X(J)))
20          CONTINUE
          CALL FX(NVAR,FPAR,IPAR,XTEMP,WORK2,IERROR)
          IWORK(22)=IWORK(22)+1
          IF(IERROR.NE.0)RETURN
          ENDIF
        DO 30 J=KCALL,NEQN,MBAND
          ILO=MAX(1,J-MU)
          IHI=MIN(NEQN,J+ML)
          IROW=ILO-J+ML+MU+1
          IBAND=IHI-ILO+1
          SKALE=-ONE
          CALL DAXPY(IBAND,SKALE,WORK2(ILO),1,WORK1(ILO),1)
          SKALE=ONE/(XJAC*EPS*(ONE+ABS(X(J))))
          CALL DSCAL(IBAND,SKALE,WORK1(ILO),1)
          SKALE=ONE
          CALL DAXPY(IBAND,SKALE,WORK1(ILO),1,FPRIME(IROW,J),1)
30        CONTINUE
40      CONTINUE
C
C  Compute last column of jacobian, rows 1 to NEQN
C
      CALL DCOPY(NVAR,X,1,XTEMP,1)
      XTEMP(NVAR)=X(NVAR)+EPS*(ONE+ABS(X(NVAR)))
      CALL FX(NVAR,FPAR,IPAR,XTEMP,WORK1,IERROR)
      IWORK(22)=IWORK(22)+1
      IF(IERROR.NE.0)RETURN
      IF(JAC.EQ.2)THEN
        XTEMP(NVAR)=X(NVAR)-EPS*(ONE+ABS(X(NVAR)))
        CALL FX(NVAR,FPAR,IPAR,XTEMP,WORK2,IERROR)
        IWORK(22)=IWORK(22)+1
        IF(IERROR.NE.0)RETURN  
        ENDIF
      SKALE=-ONE
      CALL DAXPY(NEQN,SKALE,WORK2,1,WORK1,1)
      SKALE=ONE/(XJAC*EPS*(ONE+ABS(X(NVAR))))
      CALL DSCAL(NEQN,SKALE,WORK1,1)
      SKALE=ONE
      CALL DAXPY(NEQN,SKALE,WORK1,1,FCOL,1)
C
C  Do last row, J=1,NVAR
C
      FROW(IPC)=FROW(IPC)+ONE
      RETURN
      END
      FUNCTION REPS()
C
C********************************************************************
C
C  Compute the machine relative precision.
C
C  If you would prefer a simpler computation, you can replace this routine 
C  by a line that reads "REPS=1.0E-7" for example, or "REPS=R1MACH(4)"
C  if you have the PORT/SLATEC machine constant routines.  
C
      DOUBLE PRECISION HALF
      DOUBLE PRECISION ONE
      DOUBLE PRECISION TWO
      DOUBLE PRECISION ZERO
C
      PARAMETER (HALF=0.5)
      PARAMETER (ONE=1.0)
      PARAMETER (TWO=2.0)
      PARAMETER (ZERO=0.0)
C
      DOUBLE PRECISION EPS
      INTEGER   I
      DOUBLE PRECISION RADD
      DOUBLE PRECISION RADDOLD
      DOUBLE PRECISION REPS
      DOUBLE PRECISION RMANT
      DOUBLE PRECISION RMANTOLD
      DOUBLE PRECISION TEMP1
      DOUBLE PRECISION TEMP2
C
      RMANT=ONE
      RADD=ONE
      DO 10 I=1,100
        RADDOLD=RADD
        RADD=HALF*RADD
        RMANTOLD=RMANT
        RMANT=RMANT+RADD
        IF((RMANT.EQ.RMANTOLD).OR.
     *     (RMANT-RADD.NE.RMANTOLD))GO TO 20
        IF(RMANT.EQ.TWO)THEN
          RADDOLD=RADD
          GO TO 20
          ENDIF
10      CONTINUE
C
C  Loop doesn't terminate after 100 steps!
C  We'll just use that value.
C
20    CONTINUE
C
      RMANT=RMANTOLD
      EPS=RADDOLD
C
C  Additional requirement on the value of EPS is that
C  (1.0+EPS)-EPS.EQ.(1.0-EPS)+EPS.EQ.1.0)
C
30    CONTINUE
      TEMP1=ONE+EPS
      TEMP2=ONE-EPS
      IF(TEMP1-EPS.NE.ONE.OR.TEMP2+EPS.NE.ONE)THEN
        EPS=TWO*EPS
        GO TO 30
        ENDIF
C
      REPS=EPS
      RETURN
      END 
