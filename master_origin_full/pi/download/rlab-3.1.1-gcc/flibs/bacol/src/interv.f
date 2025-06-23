      SUBROUTINE BINTERV ( XT, LXT, X, ILEFT, MFLAG )
C-----------------------------------------------------------------------
C THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE
C EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.
C SEE REFERENCE BELOW.
C
C COMPUTES LARGEST ILEFT IN (1,LXT) SUCH THAT XT(ILEFT) .LE. X.  THE
C PROGRAM STARTS THE SEARCH FOR ILEFT WITH THE VALUE OF ILEFT THAT WAS
C RETURNED AT THE PREVIOUS CALL (AND WAS SAVED IN THE LOCAL VARIABLE
C ILO) TO MINIMIZE THE WORK IN THE COMMON CASE THAT THE VALUE OF X ON
C THIS CALL IS CLOSE TO THE VALUE OF X ON THE PREVIOUS CALL.  SHOULD
C THIS ASSUMPTION NOT BE VALID, THEN THE PROGRAM LOCATES ILO AND IHI
C SUCH THAT XT(ILO) .LE. X .LT. XT(IHI) AND, ONCE THEY ARE FOUND USES
C BISECTION TO FIND THE CORRECT VALUE FOR ILEFT.
C
C LAST MODIFIED BY RONG WANG, JAN 9, 2001.
C
C REFERENCE
C
C    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.
C      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.
C
C PACKAGE ROUTINES CALLED..  NONE
C USER ROUTINES CALLED..     NONE
C CALLED BY..                COLPNT,INITAL,VALUES
C FORTRAN FUNCTIONS USED..   NONE
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS
      INTEGER LXT,ILEFT,MFLAG
      DOUBLE PRECISION XT(LXT),X
C-----------------------------------------------------------------------
C LOCAL VARIABLES
      INTEGER ILO,IHI,ISTEP,MIDDLE

      SAVE


C-----------------------------------------------------------------------
      IF(MFLAG.EQ.-2) ILO = 1
      IHI = ILO + 1
      IF (IHI .LT. LXT) GO TO 20
      IF (X .GE. XT(LXT)) GO TO 110
      IF (LXT .LE. 1) GO TO 90
      ILO = LXT - 1
      GO TO 21
   20 IF (X .GE. XT(IHI)) GO TO 40
   21 IF (X .GE. XT(ILO)) GO TO 100
C-----------------------------------------------------------------------
C NOW X .LT. XT(IHI).  FIND LOWER BOUND.
C-----------------------------------------------------------------------
   30 ISTEP = 1
   31 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO .LE. 1) GO TO 35
      IF (X .GE. XT(ILO)) GO TO 50
      ISTEP = ISTEP*2
      GO TO 31
   35 ILO = 1
      IF (X .LT. XT(1)) GO TO 90
      GO TO 50
C-----------------------------------------------------------------------
C NOW X .GE. XT(ILO).  FIND UPPER BOUND.
C-----------------------------------------------------------------------
   40 ISTEP = 1
   41 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI .GE. LXT) GO TO 45
      IF (X .LT. XT(IHI)) GO TO 50
      ISTEP = ISTEP*2
      GO TO 41
   45 IF (X .GE. XT(LXT)) GO TO 110
      IHI = LXT
C-----------------------------------------------------------------------
C NOW XT(ILO) .LE. X .LT. XT(IHI).  NARROW THE INTERVAL.
C-----------------------------------------------------------------------
   50 MIDDLE = (ILO + IHI)/2
      IF (MIDDLE .EQ. ILO) GO TO 100
C-----------------------------------------------------------------------
C NOTE..  IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1.
C-----------------------------------------------------------------------
      IF (X .LT. XT(MIDDLE)) GO TO 53
      ILO = MIDDLE
      GO TO 50
   53 IHI = MIDDLE
      GO TO 50
C-----------------------------------------------------------------------
C SET OUTPUT AND RETURN.
C-----------------------------------------------------------------------

   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END
