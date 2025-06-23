      SUBROUTINE BBSPLVD ( XT, K, X, ILEFT, VNIKX, NDERIV )
C-----------------------------------------------------------------------
C THIS SUBROUTINE IS PART OF THE B-SPLINE PACKAGE FOR THE STABLE
C EVALUATION OF ANY B-SPLINE BASIS FUNCTION OR DERIVATIVE VALUE.
C SEE REFERENCE BELOW.
C
C CALCULATES THE VALUE AND THE FIRST NDERIV-1 DERIVATIVES OF ALL
C B-SPLINES WHICH DO NOT VANISH AT X.  THE ROUTINE FILLS THE TWO-
C DIMENSIONAL ARRAY VNIKX(J,IDERIV), J=IDERIV, ... ,K WITH NONZERO
C VALUES OF B-SPLINES OF ORDER K+1-IDERIV, IDERIV=NDERIV, ... ,1, BY
C REPEATED CALLS TO BBSPLVN.
C
C LAST MODIFIED BY RONG WANG, JAN 8, 2001.
C
C REFERENCE
C
C    DEBOOR, C., PACKAGE FOR CALCULATING WITH B-SPLINES, SIAM J.
C      NUMER. ANAL., VOL. 14, NO. 3, JUNE 1977, PP. 441-472.
C
C PACKAGE ROUTINES CALLED..  BBSPLVN
C USER ROUTINES CALLED..     NONE
C CALLED BY..                COLPNT,INITAL,VALUES
C FORTRAN FUNCTIONS USED..   DBLE,MAX
C-----------------------------------------------------------------------
C SUBROUTINE PARAMETERS
      INTEGER K,NDERIV,ILEFT
      DOUBLE PRECISION X
      DOUBLE PRECISION XT(*),VNIKX(K,NDERIV)
C-----------------------------------------------------------------------
C LOCAL VARIABLES
      INTEGER KO,IDERIV,IDERVM,KMD,JM1,IPKMD,JLOW
      DOUBLE PRECISION A(20,20)
      DOUBLE PRECISION FKMD,DIFF,V
C-----------------------------------------------------------------------
C CONSTANT
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO = 0.D0)
      PARAMETER (ONE  = 1.D0)
C-----------------------------------------------------------------------
C LOOP INDICES
      INTEGER I,J,M,L
C-----------------------------------------------------------------------
      KO = K + 1 - NDERIV
      CALL BBSPLVN(XT,KO,1,X,ILEFT,VNIKX(NDERIV,NDERIV))
      IF (NDERIV .LE. 1) GO TO 120
      IDERIV = NDERIV
      DO 20 I=2,NDERIV
        IDERVM = IDERIV-1
        DO 10 J=IDERIV,K
   10     VNIKX(J-1,IDERVM) = VNIKX(J,IDERIV)
        IDERIV = IDERVM
        CALL BBSPLVN(XT,0,2,X,ILEFT,VNIKX(IDERIV,IDERIV))
   20   CONTINUE
      DO 40 I=1,K
        DO 30 J=1,K
   30     A(I,J) = ZERO
   40   A(I,I) = ONE
      KMD = K
      DO 110 M=2,NDERIV
        KMD = KMD - 1
        FKMD =  DBLE(KMD)
        I = ILEFT
        J = K
   50   JM1 = J-1
        IPKMD = I + KMD
        DIFF = XT(IPKMD) -XT(I)
        IF (JM1 .EQ. 0) GO TO 80
        IF (DIFF .EQ. ZERO) GO TO 70
        DO 60 L=1,J
   60     A(L,J) = (A(L,J) - A(L,J-1))/DIFF*FKMD
   70   J = JM1
        I = I - 1
        GO TO 50
   80   IF (DIFF .EQ. ZERO) GO TO 90
        A(1,1) = A(1,1)/DIFF*FKMD
   90   DO 110 I=1,K
          V = ZERO
          JLOW = MAX(I,M)
          DO 100 J=JLOW,K
  100       V = A(I,J)*VNIKX(J,M) + V
  110     VNIKX(I,M) = V
  120 RETURN
      END
