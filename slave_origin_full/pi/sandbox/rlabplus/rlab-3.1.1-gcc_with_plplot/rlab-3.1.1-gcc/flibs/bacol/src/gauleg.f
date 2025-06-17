
      SUBROUTINE GAULEG(N, NSQ, PTS, WTS, WORK, FLAG)
C $ID: GAULEG.F,V 1.11 1992/06/25 15:09:31 KEAST EXP $
C LAST MODIFIED BY RONG WANG, 2001/03/08
C
C     GAULEG RETURNS THE POINTS AND WEIGHTS FOR GAUSS-LEGENDRE
C     QUADRATURE OR GAUSS-LOBATTO QUADRATURE OVER THE INTERVAL 
C     [-1,1] OR [0,1].
C
C     ON INPUT:
C        
C        N     : THE NUMBER OF GAUSS-LEGENDRE POINTS.
C        NSQ   : EQUAL TO N*N, TO HANDLE STORAGE FOR EIGENVECTORS.
C        PTS   : DOUBLE PRECISION (N).
C        WTS   : DOUBLE PRECISION (N).  IF FLAG (SEE BELOW) IS 1 OR 3,
C                WTS IS USED ONLY FOR TEMPORARY WORKSPACE.
C        WORK  : DOUBLE PRECISION (NSQ), WORK SPACE FOR CALL TO IMTQL2,
C                IF WEIGHTS ARE REQUIRED. IF FLAG IS 1 OR 3, WORK IS
C                NOT REFERENCED, AND MAY BE DECLARED AS SCALAR IN THE 
C                CALLING PROGRAM.
C        FLAG  : SPECIFIES WHETHER WEIGHTS ARE ALSO REQUIRED, AND
C                WHETHER GAUSS-LEGENDRE OR LOBATTO POINTS ARE WANTED.
C                   FLAG  = 1: GAUSS-LEGENDRE POINTS ONLY OVER [-1,1];
C                         = 2: GAUSS-LEGENDRE POINTS ONLY OVER [0,1];
C                         = 3: GAUSS-LEGENDRE POINTS AND WEIGHTS OVER
C                              [-1,1];
C                         = 4: GAUSS-LEGENDRE POINTS AND WEIGHTS OVER
C                              [0,1];
C                         = 5: LOBATTO POINTS ONLY OVER [-1,1];
C                         = 6: LOBATTO POINTS ONLY OVER [0,1];
C                         = 7: LOBATTO POINTS AND WEIGHTS OVER [-1,1];
C                         = 8: LOBATTO POINTS AND WEIGHTS OVER [0,1];
C                FOR ANY OTHER VALUE, THE DEFAULT IS GAUSS-LEGENDRE
C                POINTS AND WEIGHTS OVER [-1,1].
C
C     ON OUTPUT:
C
C        FOR FLAG <> 1 OR 2:
C          PTS : PTS(I) IS THE ITH GAUSS-LEGENDRE POINT IN [-1,1],
C                PTS(I) < PTS(I+1), I = 1,2,..,N-1.
C          WTS : WTS(I) IS THE ITH GAUSS-LEGENDRE WEIGHT IF FLAG <> 1.
C
C        FOR FLAG = 3 OR 4:
C          PTS : PTS(I) IS THE ITH LOBATTO POINT IN [-1,1],
C                PTS(I) < PTS(I+1), I = 1,2,..,N-1.
C                CLEARLY, PTS(1) = -1.0, PTS(N) = 1.0.
C          WTS : WTS(I) IS THE ITH LOBATTO WEIGHT IF FLAG = 4.
C
C        WORK  : WORK, USED TO STORE EIGENVECTORS, UNREFERENCED IF FLAG
C                IS 1 OR 3.
C
C     SUBROUTINES USED:
C
C        IMQTL1: EISPACK ROUTINE TO COMPUTE THE EIGENVALUES OF A 
C                SYMMETRIC TRIDIAGONAL MATRIX.
C
C        IMQTL2: EISPACK ROUTINE TO COMPUTE THE EIGNEVECTORS AND 
C                EIGENVALUES OF A SYMMETRIC TRIDIAGONAL MATRIX.
C
C     FUNCTIONS USED:
C
C        PYTHAG: EISPACK FUNCTION TO COMPUTE EUCLIDEAN NORM.
C
C     INTRINSIC FUNCTIONS USED:
C
C        SQRT, DBLE.
C
C     VERSION: JUNE 22 1992, PAT KEAST.
C
C     DECLARATIONS:
C
C        PARAMETERS:
C
      INTEGER          N, NSQ, FLAG
      DOUBLE PRECISION PTS(N), WTS(N), WORK(NSQ)
C
C        LOCAL VARIABLES:

      INTEGER          IFAIL, J, NM2
      DOUBLE PRECISION FOUR, THREE, TWO, ONE, ZERO
*     .. EXTERNAL FUNCTIONS ..
      EXTERNAL         IMTQL2 
      PARAMETER ( FOUR = 4.0D0, THREE = 3.0D0, TWO = 2.0D0, 
     *            ONE = 1.0D0, ZERO = 0.0D0 )
  
      DO 10 J = 1,N
         PTS(J) = ZERO
   10 CONTINUE
 
      DO 20 J = 1,NSQ
         WORK(J) = ZERO
   20 CONTINUE

      DO 30 J = 1,N
         WORK((J-1)*N+J) = ONE
   30 CONTINUE
 
      IF ( FLAG .LE.4 .OR. FLAG .GT. 8 ) THEN
C        GAUS-LEGENDRE POINTS AND WEIGHTS.
         DO 40 J = 1,N-1
            WTS(J+1) = DBLE(J)/SQRT(DBLE(4*J*J)-ONE)
   40    CONTINUE
 
         IF ( FLAG .EQ. 1 .OR. FLAG .EQ. 2 ) THEN
C           COMPUTE ONLY THE GAUSS-LEGENDRE POINTS OVER [-1,1].
            CALL IMTQL1(N, PTS, WTS, IFAIL )
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 2 ) THEN
               DO 45 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
   45          CONTINUE
            ENDIF
         ELSE
C           COMPUTE BOTH POINTS AND WEIGHTS OVER [-1,1].
            IFAIL = 1
C
            CALL IMTQL2(N, N, PTS, WTS, WORK, IFAIL)
    
            DO 50 J = 1,N
               WTS(J) = TWO*WORK((J-1)*N+1)*WORK((J-1)*N+1)
   50       CONTINUE
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 4 ) THEN
               DO 55 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
                  WTS(J) = WTS(J)/TWO
   55          CONTINUE
            ENDIF
         ENDIF
C
      ELSE
C        THE LOBATTO POINTS AND WEIGHTS.

C        FIRST, COMPUTE THE ORDER N-2 JACOBI POINTS AND/OR WEIGHTS.
         NM2 = N-2
         DO 60 J = 1,NM2-1
            WTS(J+1) = SQRT(DBLE(J*(J+2))/DBLE((2*J+1)*(2*J+3)))
   60    CONTINUE

         IF ( FLAG .EQ. 5 .OR. FLAG .EQ. 6) THEN
C           COMPUTE ONLY THE GAUSS-LOBATTO POINTS OVER [-1,1].
            CALL IMTQL1(NM2, PTS(2), WTS, IFAIL )
            PTS(1) = -ONE
            PTS(N) = ONE
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 6 ) THEN
               DO 65 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
   65          CONTINUE
            ENDIF
         ELSE
C           COMPUTE BOTH POINTS AND WEIGHTS.
            IFAIL = 1
C
            CALL IMTQL2(NM2, NM2, PTS(2), WTS, WORK, IFAIL)
            PTS(1) = -ONE
            PTS(N) = ONE
    
            DO 70 J = 2,N-1
               WTS(J) = (FOUR/THREE)*WORK((J-2)*NM2+1)*WORK((J-2)*NM2+1)
     *                  /(ONE - PTS(J)*PTS(J))
   70       CONTINUE
            WTS(1) = ZERO
            DO 80 J = 2,N-1
               WTS(1) = WTS(1) - WTS(J)
   80       CONTINUE
            WTS(1) = ONE + WTS(1)/TWO
            WTS(N) = WTS(1)
C           SCALE THE VALUES TO THE INTERVAL [0,1].
            IF ( FLAG .EQ. 8 ) THEN
               DO 85 J = 1,N
                  PTS(J) = (PTS(J)+ONE)/TWO
                  WTS(J) = WTS(J)/TWO
   85          CONTINUE
            ENDIF
         ENDIF
C
         RETURN
      ENDIF
*
      RETURN
      END
