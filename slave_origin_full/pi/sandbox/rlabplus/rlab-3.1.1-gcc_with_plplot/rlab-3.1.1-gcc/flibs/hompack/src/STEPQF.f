        SUBROUTINE STEPQF(N,NFE,IFLAG,START,CRASH,HOLD,H,
     &            WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Q,R,
     &            F0,F1,Z0,DZ,W,T,SSPAR,PAR,IPAR)
C
C SUBROUTINE  STEPQF  TAKES ONE STEP ALONG THE ZERO CURVE OF THE 
C HOMOTOPY MAP  RHO(LAMBDA,X)  USING A PREDICTOR-CORRECTOR ALGORITHM.
C THE PREDICTOR USES A HERMITE CUBIC INTERPOLANT, AND THE CORRECTOR 
C RETURNS TO THE ZERO CURVE USING A QUASI-NEWTON ALGORITHM, REMAINING
C IN A HYPERPLANE PERPENDICULAR TO THE MOST RECENT TANGENT VECTOR.
C  STEPQF  ALSO ESTIMATES A STEP SIZE  H  FOR THE NEXT STEP ALONG THE 
C ZERO CURVE.
C
C
C ON INPUT:
C 
C N = DIMENSION OF  X. 
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C START = .TRUE. ON FIRST CALL TO  STEPQF, .FALSE. OTHERWISE.
C         SHOULD NOT BE MODIFIED BY THE USER AFTER THE FIRST CALL.
C
C HOLD = ||Y - YOLD|| ; SHOULD NOT BE MODIFIED BY THE USER.
C
C H = UPPER LIMIT ON LENGTH OF STEP THAT WILL BE ATTEMPTED.  H  MUST
C    BE SET TO A POSITIVE NUMBER ON THE FIRST CALL TO  STEPQF.
C    THEREAFTER,  STEPQF  CALCULATES AN OPTIMAL VALUE FOR  H, AND  H
C    SHOULD NOT BE MODIFIED BY THE USER.
C
C WK = APPROXIMATE CURVATURE FOR THE LAST STEP (COMPUTED BY PREVIOUS
C    CALL TO  STEPQF).  UNDEFINED ON FIRST CALL.  SHOULD NOT BE
C    MODIFIED BY THE USER.
C  
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION
C    IS CONSIDERED TO HAVE CONVERGED WHEN A POINT  Z=(LAMBDA,X)  IS 
C    FOUND SUCH THAT
C       ||DZ|| .LE. RELERR*||Z|| + ABSERR,
C    WHERE  DZ  IS THE LAST QUASI-NEWTON STEP.
C
C S  = (APPROXIMATE) ARC LENGTH ALONG THE HOMOTOPY ZERO CURVE UP TO
C    Y(S) = (LAMBDA(S), X(S)).
C
C Y(1:N+1) = PREVIOUS POINT (LAMBDA(S),X(S)) FOUND ON THE ZERO CURVE
C    OF THE HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  Y.  INPUT IN THIS VECTOR IS NOT USED ON THE FIRST CALL
C    TO  STEPQF.
C
C YOLD(1:N+1) = A POINT BEFORE  Y  ON THE ZERO CURVE OF THE HOMOTOPY
C    MAP.  INPUT IN THIS VECTOR IS NOT USED ON THE FIRST CALL TO 
C    STEPQF.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE 
C    HOMOTOPY MAP AT  YOLD.
C
C A(1:N) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C Q(1:N+1,1:N+1) =  Q  OF THE QR FACTORIZATION OF
C    THE AUGMENTED JACOBIAN MATRIX AT  Y.
C
C R((N+1)*(N+2)/2) = THE UPPER TRIANGLE  R  OF THE QR 
C    FACTORIZATION, STORED BY COLUMNS.
C
C F0(1:N+1), F1(1:N+1), Z0(1:N+1), DZ(1:N+1), W(1:N+1), T(1:N+1) ARE
C    WORK ARRAYS.  
C 
C SSPAR(1:4) = PARAMETERS USED FOR COMPUTATION OF THE OPTIMAL STEP SIZE.  
C    SSPAR(1) = HMIN, SSPAR(2) = HMAX, SSPAR(3) = BMIN, SSPAR(4) = BMAX.  
C    THE OPTIMAL STEP  H  IS RESTRICTED SUCH THAT 
C       HMIN .LE. H .LE. HMAX, AND  BMIN*HOLD .LE. H .LE. BMAX*HOLD.
C
C PAR(1:*)  AND  IPAR(1:*)  ARE ARRAYS FOR (OPTIONAL) USER PARAMETERS,
C    WHICH ARE SIMPLY PASSED THROUGH TO THE USER WRITTEN SUBROUTINES
C    RHO, RHOJAC.
C
C
C ON OUTPUT:
C
C NFE HAS BEEN UPDATED.
C
C IFLAG
C
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF A JACOBIAN MATRIX WITH RANK <  N  HAS OCCURRED.  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE ITERATION FAILED TO CONVERGE. 
C
C START = .FALSE. ON A NORMAL RETURN.
C
C CRASH 
C
C    = .FALSE. ON A NORMAL RETURN.
C
C    = .TRUE. IF THE STEP SIZE  H  WAS TOO SMALL.  H  HAS BEEN
C      INCREASED TO AN ACCEPTABLE VALUE, WITH WHICH  STEPQF  MAY BE
C      CALLED AGAIN.
C
C    = .TRUE. IF  RELERR  AND/OR  ABSERR  WERE TOO SMALL.  THEY HAVE
C      BEEN INCREASED TO ACCEPTABLE VALUES, WITH WHICH  STEPQF  MAY
C      BE CALLED AGAIN.
C
C HOLD = ||Y-YOLD||.
C
C H = OPTIMAL VALUE FOR NEXT STEP TO BE ATTEMPTED.  NORMALLY  H  SHOULD
C     NOT BE MODIFIED BY THE USER.
C
C WK = APPROXIMATE CURVATURE FOR THE STEP TAKEN BY  STEPQF.
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE ZERO CURVE OF THE HOMOTOPY 
C     MAP UP TO THE LATEST POINT FOUND, WHICH IS RETURNED IN  Y.
C
C RELERR, ABSERR  ARE UNCHANGED ON A NORMAL RETURN.  THEY ARE POSSIBLY
C     CHANGED IF  CRASH  = .TRUE. (SEE DESCRIPTION OF  CRASH  ABOVE).
C
C Y, YP, YOLD, YPOLD  CONTAIN THE TWO MOST RECENT POINTS AND TANGENT
C     VECTORS FOUND ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C Q, R  STORE THE QR FACTORIZATION OF THE AUGMENTED JACOBIAN MATRIX 
C     EVALUATED AT  Y.
C
C
C CALLS  D1MACH, DAXPY, DCOPY, DDOT, DGEMV, DGEQRF, DNRM2, DORGQR, 
C     DSCAL, DTPSV, F (OR RHO), FJAC (OR RHOJAC), TANGQF, UPQRQF.
C
C ***** DECLARATIONS *****
C
C     FUNCTION DECLARATIONS  
C
        DOUBLE PRECISION D1MACH, DDOT, DNRM2, QOFS
C
C     LOCAL VARIABLES
C
        DOUBLE PRECISION ALPHA, DD001, DD0011, DD01, DD011, DELS, ETA, 
     &    FOURU, GAMMA, HFAIL, HTEMP, IDLERR, ONE, P0, P1, PP0, PP1, 
     &    TEMP, TWOU, WKOLD, ZERO         
        INTEGER I, ITCNT, LITFH, J, JP1, NP1
        LOGICAL FAILED
C
C     SCALAR ARGUMENTS 
C
        INTEGER N, NFE, IFLAG
        LOGICAL START, CRASH
        DOUBLE PRECISION HOLD, H, WK, RELERR, ABSERR, S
C
C     ARRAY DECLARATIONS
C
        DOUBLE PRECISION Y(N+1), YP(N+1), YOLD(N+1), YPOLD(N+1), 
     &    A(N), Q(N+1,N+1), R((N+1)*(N+2)/2), F0(N+1), F1(N+1),
     &    Z0(N+1), DZ(N+1), W(N+1), T(N+1), SSPAR(4), PAR(1)
        INTEGER IPAR(1)
C
        SAVE
C
C ***** END OF DECLARATIONS *****
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
        DD01(P0,P1,DELS) = (P1-P0)/DELS
        DD001(P0,PP0,P1,DELS) = (DD01(P0,P1,DELS)-PP0)/DELS
        DD011(P0,P1,PP1,DELS) = (PP1-DD01(P0,P1,DELS))/DELS
        DD0011(P0,PP0,P1,PP1,DELS) = (DD011(P0,P1,PP1,DELS) -
     &                                DD001(P0,PP0,P1,DELS))/DELS
        QOFS(P0,PP0,P1,PP1,DELS,S) = ((DD0011(P0,PP0,P1,PP1,DELS)*
     &    (S-DELS) + DD001(P0,PP0,P1,DELS))*S + PP0)*S + P0
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
C
C ***** INITIALIZATION *****
C
C ETA = PARAMETER FOR BROYDEN'S UPDATE.
C LITFH = MAXIMUM NUMBER OF QUASI-NEWTON ITERATIONS ALLOWED.
C
        ONE = 1.0
        ZERO = 0.0
        TWOU = 2.0*D1MACH(4)
        FOURU = TWOU + TWOU
        NP1 = N+1
        FAILED = .FALSE.
        CRASH = .TRUE.
        ETA = 50.0*TWOU
        LITFH = 2*(INT(ABS(LOG10(ABSERR+RELERR)))+1)
C
C CHECK THAT ALL INPUT PARAMETERS ARE CORRECT.
C
C     THE ARCLENGTH  S MUST BE NONNEGATIVE.
C
        IF (S .LT. 0.0) RETURN
C
C     IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE.
C   
        IF (H .LT. FOURU*(1.0+S)) THEN
          H=FOURU*(1.0 + S)
          RETURN
        END IF
C
C     IF ERROR TOLERANCES ARE TOO SMALL, INCREASE THEM TO ACCEPTABLE 
C     VALUES.
C
        TEMP=DNRM2(NP1,Y,1) + 1.0
        IF (.5*(RELERR*TEMP+ABSERR) .LT. TWOU*TEMP) THEN
          IF (RELERR .NE. 0.0) THEN
            RELERR = FOURU*(1.0+FOURU)
            TEMP = 0.0
            ABSERR = MAX(ABSERR,TEMP)
          ELSE
            ABSERR=FOURU*TEMP
          END IF
          RETURN
        END IF
C
C     INPUT PARAMETERS WERE ALL ACCEPTABLE.
C
        CRASH = .FALSE.
C
C COMPUTE  YP  ON FIRST CALL.
C NOTE:  DZ  IS USED SIMPLY AS A WORK ARRAY HERE.
C
        IF (START) THEN
          CALL TANGQF(Y,YP,YPOLD,A,Q,R,W,DZ,T,N,IFLAG,NFE,PAR,IPAR)
          IF (IFLAG .GT. 0) RETURN
        END IF
C
C F0 = (RHO(Y), YP*Y) TRANSPOSE (DIFFERENT FOR EACH PROBLEM TYPE).
C
         IF (IFLAG .EQ. -2) THEN
C
C          CURVE TRACKING PROBLEM.
C
           CALL HOM1RHO(A,Y(1),Y(2),F0,PAR,IPAR)
         ELSE IF (IFLAG .EQ. -1) THEN
C
C          ZERO FINDING PROBLEM.
C
           CALL HOM1F(Y(2),F0)
           DO 5 I=1,N
             F0(I) = Y(1)*F0(I) + (1.0-Y(1))*(Y(I+1)-A(I))
  5        CONTINUE
         ELSE
C
C          FIXED POINT PROBLEM.
C
           CALL HOM1F(Y(2),F0)
           DO 10 I=1,N
             F0(I) = Y(1)*(A(I)-F0(I))+Y(I+1)-A(I)
  10       CONTINUE
         END IF
C
C        DEFINE LAST ROW OF F0  =  YP*Y.
C
           F0(NP1) = DDOT(NP1,YP,1,Y,1)
C
C ***** END OF INITIALIZATION *****
C
C ***** COMPUTE PREDICTOR POINT Z0 *****
C
  20    IF (START) THEN
C           
C         COMPUTE Z0 WITH LINEAR PREDICTOR USING Y, YP --
C         Z0 = Y+H*YP.
C
          CALL DCOPY(NP1,Y,1,Z0,1)
          CALL DAXPY(NP1,H,YP,1,Z0,1)
C         
        ELSE
C
C         COMPUTE Z0 WITH CUBIC PREDICTOR.
C
          DO 30 I=1,NP1
            Z0(I) = QOFS(YOLD(I),YPOLD(I),Y(I),YP(I),HOLD,HOLD+H) 
  30      CONTINUE
C
        END IF
C
C F1 = (RHO(Z0), YP*Z0) TRANSPOSE.
C
        IF (IFLAG .EQ. -2) THEN
          CALL HOM1RHO(A,Z0(1),Z0(2),F1,PAR,IPAR)
        ELSE IF (IFLAG .EQ. -1) THEN
          CALL HOM1F(Z0(2),F1)
          DO 40 I=1,N
            F1(I) = Z0(1)*F1(I) + (1.0-Z0(1))*(Z0(I+1)-A(I))
  40      CONTINUE
        ELSE
          CALL HOM1F(Z0(2),F1)
          DO 50 I=1,N
            F1(I) = Z0(1)*(A(I)-F1(I))+Z0(I+1)-A(I)
  50      CONTINUE
        END IF
        F1(NP1) = DDOT(NP1,YP,1,Z0,1)
C
C ***** END OF PREDICTOR SECTION *****
C
C ***** SET-UP FOR QUASI-NEWTON ITERATION *****
C
        IF (FAILED) THEN
C        
C GENERATE Q = AUGMENTED JACOBIAN MATRIX FOR POINT Z0=(LAMBDA,X).
C        
          IF (IFLAG .EQ. -2) THEN
C
C           CURVE TRACKING PROBLEM:
C           D(RHO) = (D RHO(A,LAMBDA,X)/D LAMBDA, D RHO(A,LAMBDA,X)/DX).
C
            DO 60 J = 1,NP1
              CALL HOM1RHOJ(A,Z0(1),Z0(2),Q(1,J),J,PAR,IPAR)
  60        CONTINUE
          ELSE IF (IFLAG .EQ. -1) THEN
C
C           ZERO FINDING PROBLEM:
C           D(RHO) = (F(X) - X + A, LAMBDA*DF(X) + (1-LAMBDA)*I).
C
            CALL HOM1F(Z0(2),Q(1,1))
            DO 70 I=1,N
              Q(I,1) = A(I) - Z0(I+1) + Q(I,1)
  70        CONTINUE
            DO 80 J= 1,N
              JP1 = J+1
              CALL HOM1FJAC(Z0(2),Q(1,JP1),J)
              CALL DSCAL(N, Z0(1), Q(1,JP1), 1)
              Q(J,JP1) = 1.0 - Z0(1) + Q(J,JP1)
  80        CONTINUE
          ELSE 
C 
C           FIXED POINT PROBLEM:
C           D(RHO) = (A - F(X), I - LAMBDA*DF(X)).
C
            CALL HOM1F(Z0(2),Q(1,1))
            CALL DSCAL(N,-ONE,Q(1,1),1)
            CALL DAXPY(N,ONE,A,1,Q(1,1),1)
            DO 90 J=1,N
              JP1 = J+1
              CALL HOM1FJAC(Z0(2),Q(1,JP1),J)
              CALL DSCAL(N, -Z0(1), Q(1,JP1), 1)
              Q(J,JP1) = 1.0 + Q(J,JP1)
  90        CONTINUE
          END IF
C
C       DEFINE LAST ROW OF Q = YP.
C
          CALL DCOPY(NP1, YP, 1, Q(NP1,1), NP1)
C
C       COUNT JACOBIAN EVALUATION.
C
          NFE = NFE+1
C
C DO FIRST QUASI-NEWTON STEP.
C
C       FACTOR AUG.
C
          CALL DGEQRF(NP1,NP1,Q,NP1,T,W,NP1,I)
C
C       PACK UPPER TRIANGLE INTO ARRAY R.
C
          DO 100  I=1,NP1
            CALL DCOPY(I,Q(1,I),1,R((I*(I-1))/2 + 1),1)
 100      CONTINUE
C
C       CHECK FOR SINGULARITY.
C
          J = 1
          DO 110 I = 1, N
            IF( R(J+I-1) .EQ. ZERO ) THEN
              IFLAG = 4
              RETURN
            END IF
            J = J + I
 110      CONTINUE
C
C       EXPAND HOUSEHOLDER REFLECTIONS INTO FULL MATRIX Q .
C
          CALL DORGQR(NP1, NP1, N, Q, NP1, T, W, NP1, I)
C
C       COMPUTE NEWTON STEP.
C
          CALL DCOPY(N,F1,1,T,1)
          CALL DSCAL(N,-ONE,T,1)
          T(NP1) = 0.0
          CALL DGEMV('T',NP1,NP1,ONE,Q,NP1,T,1,ZERO,DZ,1)
          CALL DTPSV('U', 'N', 'N', NP1, R, DZ, 1)
C
C       TAKE STEP AND SET F0 = F1.
C
          CALL DAXPY(NP1, ONE, DZ, 1, Z0, 1)
          CALL DCOPY(NP1, F1, 1, F0, 1)
C
C       F1 = (RHO(Z0), YP*Z0) TRANSPOSE.
C
          IF (IFLAG .EQ. -2) THEN
            CALL HOM1RHO(A,Z0(1),Z0(2),F1,PAR,IPAR)
          ELSE IF (IFLAG .EQ. -1) THEN
            CALL HOM1F(Z0(2),F1)
            DO 120 I=1,N
              F1(I) = Z0(1)*F1(I) + (1.0-Z0(1))*(Z0(I+1)-A(I))
  120       CONTINUE
          ELSE
            CALL HOM1F(Z0(2),F1)
            DO 130 I=1,N
              F1(I) = Z0(1)*(A(I)-F1(I))+Z0(I+1)-A(I)
  130       CONTINUE
          END IF
          F1(NP1) = DDOT(NP1,YP,1,Z0,1)
C
        ELSE
C
C IF NOT FAILED THEN DEFINE  DZ=Z0-Y  PRIOR TO MAIN LOOP.
C
          CALL DCOPY(NP1,Z0,1,DZ,1)
          CALL DAXPY(NP1,-ONE,Y,1,DZ,1)
        END IF
C
C ***** END OF PREPARATION FOR QUASI-NEWTON ITERATION *****
C
C ***** QUASI-NEWTON ITERATION *****
C
        DO 160 ITCNT = 1,LITFH
C
C PERFORM UPDATE FOR NEWTON STEP JUST TAKEN.
C
          CALL UPQRQF(NP1,ETA,DZ,F0,F1,Q,R,W,T)
C
C COMPUTE NEXT NEWTON STEP.
C
          CALL DCOPY(N,F1,1,T,1)
          CALL DSCAL(N,-ONE,T,1)
          T(NP1) = 0.0
          CALL DGEMV('T',NP1,NP1,ONE,Q,NP1,T,1,ZERO,DZ,1)
          CALL DTPSV('U', 'N', 'N', NP1, R, DZ, 1)
C
C TAKE STEP.
C
          CALL DAXPY(NP1, ONE, DZ, 1, Z0, 1)
C
C CHECK FOR CONVERGENCE.
C
          IF (DNRM2(NP1,DZ,1) .LE. RELERR*DNRM2(NP1,Z0,1)+ABSERR) THEN
             GO TO 180
          END IF
C
C IF NOT CONVERGED, PREPARE FOR NEXT ITERATION.
C
C       F0 = F1.
C
          CALL DCOPY(NP1, F1, 1, F0, 1)
C
C       F1 = (RHO(Z0), YP*Z0) TRANSPOSE.
C
          IF (IFLAG .EQ. -2) THEN
            CALL HOM1RHO(A,Z0(1),Z0(2),F1,PAR,IPAR)
          ELSE IF (IFLAG .EQ. -1) THEN
            CALL HOM1F(Z0(2),F1)
            DO 140 I=1,N
              F1(I) = Z0(1)*F1(I) + (1.0-Z0(1))*(Z0(I+1)-A(I))
  140       CONTINUE
          ELSE
            CALL HOM1F(Z0(2),F1)
            DO 150 I=1,N
              F1(I) = Z0(1)*(A(I)-F1(I))+Z0(I+1)-A(I)
  150       CONTINUE
          END IF
          F1(NP1) = DDOT(NP1,YP,1,Z0,1)
C
  160   CONTINUE
C
C ***** END OF QUASI-NEWTON LOOP *****
C
C ***** DIDN'T CONVERGE OR TANGENT AT NEW POINT DID NOT MAKE
C       AN ACUTE ANGLE WITH YPOLD -- TRY AGAIN WITH A SMALLER H *****
C      
  170   FAILED = .TRUE.
        HFAIL = H
        IF (H .LE. FOURU*(1.0 + S)) THEN
          IFLAG = 6
          RETURN
        ELSE
          H = .5 * H
        END IF
        GO TO 20
C
C ***** END OF CONVERGENCE FAILURE SECTION *****
C
C ***** CONVERGED -- MOP UP AND RETURN *****
C
C COMPUTE TANGENT & AUGMENTED JACOBIAN AT  Z0.
C NOTE:  DZ  AND  F1  ARE USED SIMPLY AS WORK ARRAYS HERE.
C
  180   CALL TANGQF(Z0,T,YP,A,Q,R,W,DZ,F1,N,IFLAG,NFE,PAR,IPAR)
        IF (IFLAG .GT. 0) RETURN
C
C CHECK THAT COMPUTED TANGENT  T  MAKES AN ANGLE NO LARGER THAN
C 60 DEGREES WITH CURRENT TANGENT  YP.  (I.E. COS OF ANGLE < .5)
C IF NOT, STEP SIZE WAS TOO LARGE, SO THROW AWAY Z0, AND TRY
C AGAIN WITH A SMALLER STEP.
C
        ALPHA = DDOT(NP1,T,1,YP,1)
        IF (ALPHA .LT. 0.5) GOTO 170
        ALPHA = ACOS(ALPHA)
C
C SET UP VARIABLES FOR NEXT CALL.
C
        CALL DCOPY(NP1,Y,1,YOLD,1)
        CALL DCOPY(NP1,Z0,1,Y,1)
        CALL DCOPY(NP1,YP,1,YPOLD,1)
        CALL DCOPY(NP1,T,1,YP,1)
C
C UPDATE ARCLENGTH   S = S + ||Y-YOLD||.
C
        HTEMP = HOLD
        CALL DAXPY(NP1,-ONE,YOLD,1,Z0,1)
        HOLD = DNRM2(NP1,Z0,1)
        S = S+HOLD
C
C COMPUTE OPTIMAL STEP SIZE. 
C   IDLERR = DESIRED ERROR FOR NEXT PREDICTOR STEP.
C   WK = APPROXIMATE CURVATURE = 2*SIN(ALPHA/2)/HOLD  WHERE 
C        ALPHA = ARCCOS(YP*YPOLD).
C   GAMMA = EXPECTED CURVATURE FOR NEXT STEP, COMPUTED BY 
C        EXTRAPOLATING FROM CURRENT CURVATURE  WK, AND LAST 
C        CURVATURE  WKOLD.  GAMMA IS FURTHER REQUIRED TO BE 
C        POSITIVE.
C
        WKOLD = WK
        IDLERR = SQRT(SQRT(ABSERR + RELERR*DNRM2(NP1,Y,1)))
C
C     IDLERR SHOULD BE NO BIGGER THAN 1/2 PREVIOUS STEP.
C
        IDLERR = MIN(.5*HOLD,IDLERR)
        WK = 2.0*ABS(SIN(.5*ALPHA))/HOLD
        IF (START) THEN
           GAMMA = WK
        ELSE 
           GAMMA = WK + HOLD/(HOLD+HTEMP)*(WK-WKOLD)
        END IF
        GAMMA = MAX(GAMMA, 0.01*ONE)
        H = SQRT(2.0*IDLERR/GAMMA)
C
C     ENFORCE RESTRICTIONS ON STEP SIZE SO AS TO ENSURE STABILITY.
C        HMIN <= H <= HMAX, BMIN*HOLD <= H <= BMAX*HOLD.
C
        H = MIN(MAX(SSPAR(1),SSPAR(3)*HOLD,H),SSPAR(4)*HOLD,SSPAR(2))
        IF (FAILED) H = MIN(HFAIL,H)
        START = .FALSE.
C
C ***** END OF MOP UP SECTION *****
C
        RETURN
C
C ***** END OF SUBROUTINE STEPQF *****
        END
