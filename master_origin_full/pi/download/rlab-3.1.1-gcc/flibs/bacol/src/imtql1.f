      SUBROUTINE IMTQL1(N,D,E,IERR) 
C 
      INTEGER I,J,L,M,N,II,MML,IERR 
      DOUBLE PRECISION D(N),E(N) 
      DOUBLE PRECISION B,C,F,G,P,R,S,TST1,TST2,PYTHAG 
C 
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE IMTQL1, 
C     NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON, 
C     AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE. 
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971). 
C 
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC 
C     TRIDIAGONAL MATRIX BY THE IMPLICIT QL METHOD. 
C 
C     ON INPUT 
C 
C        N IS THE ORDER OF THE MATRIX. 
C 
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX. 
C 
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX 
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY. 
C 
C      ON OUTPUT 
C 
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN 
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND 
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE 
C          THE SMALLEST EIGENVALUES. 
C 
C        E HAS BEEN DESTROYED. 
C 
C        IERR IS SET TO 
C          ZERO       FOR NORMAL RETURN, 
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN 
C                     DETERMINED AFTER 30 ITERATIONS. 
C 
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) . 
C 
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW, 
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
C 
C     THIS VERSION DATED APRIL 1983. 
C 
C     ------------------------------------------------------------------ 
C 
      IERR = 0 
      IF (N .EQ. 1) GO TO 1001 
C 
      DO 100 I = 2, N 
  100 E(I-1) = E(I) 
C 
      E(N) = 0.0D0 
C 
      DO 290 L = 1, N 
         J = 0 
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT .......... 
  105    DO 110 M = L, N 
            IF (M .EQ. N) GO TO 120 
            TST1 = DABS(D(M)) + DABS(D(M+1)) 
            TST2 = TST1 + DABS(E(M)) 
            IF (TST2 .EQ. TST1) GO TO 120 
  110    CONTINUE 
C 
  120    P = D(L) 
         IF (M .EQ. L) GO TO 215 
         IF (J .EQ. 30) GO TO 1000 
         J = J + 1 
C     .......... FORM SHIFT .......... 
         G = (D(L+1) - P) / (2.0D0 * E(L)) 
         R = PYTHAG(G,1.0D0) 
         G = D(M) - P + E(L) / (G + DSIGN(R,G)) 
         S = 1.0D0 
         C = 1.0D0 
         P = 0.0D0 
         MML = M - L 
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- .......... 
         DO 200 II = 1, MML 
            I = M - II 
            F = S * E(I) 
            B = C * E(I) 
            R = PYTHAG(F,G) 
            E(I+1) = R 
            S = F / R 
            C = G / R 
            G = D(I+1) - P 
            R = (D(I) - G) * S + 2.0D0 * C * B 
            P = S * R 
            D(I+1) = G + P 
            G = C * R - B 
  200    CONTINUE 
C 
         D(L) = D(L) - P 
         E(L) = G 
         E(M) = 0.0D0 
         GO TO 105 
C     .......... ORDER EIGENVALUES .......... 
  215    IF (L .EQ. 1) GO TO 250 
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- .......... 
         DO 230 II = 2, L 
            I = L + 2 - II 
            IF (P .GE. D(I-1)) GO TO 270 
            D(I) = D(I-1) 
  230    CONTINUE 
C 
  250    I = 1 
  270    D(I) = P 
  290 CONTINUE 
C 
      GO TO 1001 
C     .......... SET ERROR -- NO CONVERGENCE TO AN 
C                EIGENVALUE AFTER 30 ITERATIONS .......... 
 1000 IERR = L 
 1001 RETURN 
      END 
