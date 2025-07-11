      SUBROUTINE ROOTNF(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
     &   YPOLD,A,QR,ALPHA,TZ,PIVOT,W,WP,PAR,IPAR)
C
C ROOTNF  FINDS THE POINT  YBAR = (1, XBAR)  ON THE ZERO CURVE OF THE
C HOMOTOPY MAP.  IT STARTS WITH TWO POINTS YOLD=(LAMBDAOLD,XOLD) AND
C Y=(LAMBDA,X) SUCH THAT  LAMBDAOLD < 1 <= LAMBDA , AND ALTERNATES
C BETWEEN HERMITE CUBIC INTERPOLATION AND NEWTON ITERATION UNTIL
C CONVERGENCE.
C
C ON INPUT:
C
C N = DIMENSION OF X AND THE HOMOTOPY MAP.
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION IS
C    CONSIDERED TO HAVE CONVERGED WHEN A POINT Y=(LAMBDA,X) IS FOUND 
C    SUCH THAT
C
C    |Y(1) - 1| <= RELERR + ABSERR              AND
C
C    ||Z|| <= RELERR*||X|| + ABSERR  ,          WHERE
C
C    (?,Z) IS THE NEWTON STEP TO Y=(LAMBDA,X).
C
C Y(1:N+1) = POINT (LAMBDA(S), X(S)) ON ZERO CURVE OF HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY MAP
C    AT  Y .
C
C YOLD(1:N+1) = A POINT DIFFERENT FROM  Y  ON THE ZERO CURVE.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  YOLD .
C
C A(1:*) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C QR(1:N,1:N+2), ALPHA(1:3*N+3), TZ(1:N+1), PIVOT(1:N+1), W(1:N+1), 
C    WP(1:N+1)  ARE WORK ARRAYS USED FOR THE QR FACTORIZATION (IN THE
C    NEWTON STEP CALCULATION) AND THE INTERPOLATION.
C
C PAR(1:*) AND IPAR(1:*) ARE ARRAYS FOR (OPTIONAL) USER PARAMETERS,
C    WHICH ARE SIMPLY PASSED THROUGH TO THE USER WRITTEN SUBROUTINES
C    RHO, RHOJAC.
C
C ON OUTPUT:
C
C N , RELERR , ABSERR , A  ARE UNCHANGED.
C
C NFE  HAS BEEN UPDATED.
C
C IFLAG 
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF A JACOBIAN MATRIX WITH RANK < N HAS OCCURRED.  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE ITERATION FAILED TO CONVERGE.  Y  AND  YOLD  CONTAIN
C        THE LAST TWO POINTS FOUND ON THE ZERO CURVE.
C
C Y  IS THE POINT ON THE ZERO CURVE OF THE HOMOTOPY MAP AT  LAMBDA = 1 .
C
C
C CALLS  D1MACH ,DAXPY, DCOPY, DSCAL, DNRM2 , ROOT , TANGNF .
C
      DOUBLE PRECISION ABSERR,AERR,D1MACH, 
     &   DD001,DD0011,DD01,DD011, DELS,DNRM2,F0,F1,FP0,FP1,
     &   QOFS,QSOUT,RELERR,RERR,ONE,S,SA,SB,SOUT,U
      INTEGER IFLAG,JUDY,JW,LCODE,LIMIT,N,NFE,NP1
      LOGICAL BRACK
C
C ***** ARRAY DECLARATIONS. *****
C
      DOUBLE PRECISION Y(N+1),YP(N+1),YOLD(N+1),
     &   YPOLD(N+1),A(N),QR(N,N+2),ALPHA(3*N+3),TZ(N+1),W(N+1),
     &   WP(N+1),PAR(1)
      INTEGER PIVOT(N+1),IPAR(1)
C
C ***** END OF DIMENSIONAL INFORMATION. *****
C
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(F0,F1,DELS)=(F1-F0)/DELS
      DD001(F0,FP0,F1,DELS)=(DD01(F0,F1,DELS)-FP0)/DELS
      DD011(F0,F1,FP1,DELS)=(FP1-DD01(F0,F1,DELS))/DELS
      DD0011(F0,FP0,F1,FP1,DELS)=(DD011(F0,F1,FP1,DELS) - 
     &                            DD001(F0,FP0,F1,DELS))/DELS
      QOFS(F0,FP0,F1,FP1,DELS,S)=((DD0011(F0,FP0,F1,FP1,DELS)*(S-DELS) +
     &   DD001(F0,FP0,F1,DELS))*S + FP0)*S + F0
C
C
      ONE=1.0
      U=D1MACH(4)
      RERR=MAX(RELERR,U)
      AERR=MAX(ABSERR,0.0D0)
      NP1=N+1
C
C THE LIMIT ON THE NUMBER OF ITERATIONS ALLOWED MAY BE CHANGED BY
C CHANGING THE FOLLOWING PARAMETER STATEMENT:
      LIMIT=2*(INT(ABS(LOG10(AERR+RERR)))+1)
C
      CALL DCOPY(NP1,Y,1,TZ,1)
      CALL DAXPY(NP1,-ONE,YOLD,1,TZ,1)
      DELS=DNRM2(NP1,TZ,1)
C
C USING TWO POINTS AND TANGENTS ON THE HOMOTOPY ZERO CURVE, CONSTRUCT 
C THE HERMITE CUBIC INTERPOLANT Q(S).  THEN USE  ROOT  TO FIND THE S
C CORRESPONDING TO  LAMBDA = 1 .  THE TWO POINTS ON THE ZERO CURVE ARE
C ALWAYS CHOSEN TO BRACKET LAMBDA=1, WITH THE BRACKETING INTERVAL
C ALWAYS BEING [0, DELS].
C
      SA=0.0
      SB=DELS
      LCODE=1
130   CALL ROOT(SOUT,QSOUT,SA,SB,RERR,AERR,LCODE)
      IF (LCODE .GT. 0) GO TO 140
      QSOUT=QOFS(YOLD(1),YPOLD(1),Y(1),YP(1),DELS,SOUT) - 1.0
      GO TO 130
C IF LAMBDA = 1 WERE BRACKETED,  ROOT  CANNOT FAIL.
140   IF (LCODE .GT. 2) THEN
        IFLAG=6
        RETURN
      ENDIF
C
C CALCULATE Q(SA) AS THE INITIAL POINT FOR A NEWTON ITERATION.
      DO 150 JW=1,NP1
        W(JW)=QOFS(YOLD(JW),YPOLD(JW),Y(JW),YP(JW),DELS,SA)
150   CONTINUE
C
C ***** END OF CALCULATION OF CUBIC INTERPOLANT. *****
C
C TANGENT INFORMATION  YP  IS NO LONGER NEEDED.  HEREAFTER,  YP
C REPRESENTS THE MOST RECENT POINT WHICH IS ON THE OPPOSITE SIDE OF
C THE HYPERPLANE  LAMBDA = 1 FROM  Y.
C
C    PREPARE FOR MAIN LOOP.
C
      CALL DCOPY(NP1,YOLD,1,YP,1)
C
C INITIALIZE  BRACK  TO INDICATE THAT THE POINTS  Y  AND  YOLD  BRACKET
C LAMBDA = 1,  THUS  YOLD = YP .
C
      BRACK = .TRUE.
C
C ***** MAIN LOOP. *****
      DO 300 JUDY = 1,LIMIT
C CALCULATE NEWTON STEP AT CURRENT ESTIMATE  W .
      CALL TANGNF(SA,W,WP,YPOLD,A,QR,ALPHA,TZ,PIVOT,NFE,N,IFLAG,
     &            PAR,IPAR)
      IF (IFLAG .GT. 0) RETURN
C
C NEXT POINT = CURRENT POINT + NEWTON STEP.
C
      CALL DAXPY(NP1,ONE,TZ,1,W,1)
C
C CHECK FOR CONVERGENCE.
C
      IF ((ABS(W(1)-1.0) .LE. RERR+AERR) .AND.
     &    (DNRM2(NP1,TZ,1) .LE. RERR*DNRM2(N,W(2),1)+AERR)) THEN
        CALL DCOPY(NP1,W,1,Y,1)
        RETURN
      ENDIF
C
C PREPARE FOR NEXT ITERATION.
C
      IF (ABS(W(1)-1.0) .LE. RERR+AERR) THEN
         CALL DCOPY(NP1,WP,1,YPOLD,1)
         GOTO 300
      ENDIF
C
C    UPDATE  Y  AND  YOLD .
C
      CALL DCOPY(NP1,Y,1,YOLD,1)
      CALL DCOPY(NP1,W,1,Y,1)
C
C    UPDATE  YP  SUCH THAT  YP  IS THE MOST RECENT POINT
C    OPPOSITE OF  LAMBDA = 1  FROM  Y .  SET  BRACK = .TRUE.  IFF
C    Y  AND  YOLD  BRACKET  LAMBDA = 1  SO THAT  YP = YOLD .
C
          IF ((Y(1)-1.0)*(YOLD(1)-1.0) .GT. 0) THEN
            BRACK = .FALSE.
          ELSE
            BRACK = .TRUE.
            CALL DCOPY(NP1,YOLD,1,YP,1)
          END IF
C
C    COMPUTE  DELS = ||Y-YP||.
C
          CALL DCOPY(NP1,Y,1,TZ,1)
          CALL DAXPY(NP1,-ONE,YP,1,TZ,1)
          DELS=DNRM2(NP1,TZ,1)
C
C       COMPUTE  TZ  FOR THE LINEAR PREDICTOR   W = Y + TZ,
C           WHERE  TZ = SA*(YOLD-Y).
C
          SA = (1.0-Y(1))/(YOLD(1)-Y(1))
          CALL DCOPY(NP1,YOLD,1,TZ,1)
          CALL DAXPY(NP1,-ONE,Y,1,TZ,1)
          CALL DSCAL(NP1,SA,TZ,1)
C
C       TO INSURE STABILITY, THE LINEAR PREDICTION MUST BE NO FARTHER
C       FROM  Y  THAN  YP  IS.  THIS IS GUARANTEED IF  BRACK = .TRUE.
C       IF LINEAR PREDICTION IS TOO FAR AWAY, USE BRACKETING POINTS
C       TO COMPUTE LINEAR PREDICTION.
C
          IF (.NOT. BRACK) THEN
            IF (DNRM2(NP1,TZ,1) .GT. DELS) THEN
C
C             COMPUTE  TZ = SA*(YP-Y).
C
              SA = (1.0-Y(1))/(YP(1)-Y(1))
              CALL DCOPY(NP1,YP,1,TZ,1)
              CALL DAXPY(NP1,-ONE,Y,1,TZ,1)
              CALL DSCAL(NP1,SA,TZ,1)
            END IF
          END IF
C
C       COMPUTE ESTIMATE  W = Y + TZ  AND SAVE OLD TANGENT VECTOR.
C
           CALL DAXPY(NP1,ONE,TZ,1,W,1)
           CALL DCOPY(NP1,WP,1,YPOLD,1)
 300   CONTINUE
C
C ***** END OF MAIN LOOP. *****
C
C THE ALTERNATING OSCULATORY CUBIC INTERPOLATION AND NEWTON ITERATION
C HAS NOT CONVERGED IN  LIMIT  STEPS.  ERROR RETURN.
      IFLAG=6
      RETURN
      END
