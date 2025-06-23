C     This compilation of MJD Powell's optimization-without-derivatives
C     codes was created by M.Kostrun (mkostrun@gmail.com). It comprises
C     of two files:
C       libmjdpowell.f (this file, collection of solvers NEWUOA,BOBYQA,
C                       slightly modified, see below)
C     and
C       libmjdpowell.h (header file for linking to C and C++ codes)
C     MODIFICATIONS: (1) The original Powell's routines were modified so that
C     user can control messages created by the codes (no messages, and max
C     verbosity). This is done for explicit purpose of integrating the 
C     solvers in scripting environments (rlabplus was primary motivation)
C     where the user can choose amount of messaging, and where the messages
C     go (standard output, another terminal, or file).
C     If you link the code using really old gcc's then be aware that they
C     may cause segmentation fault if the same output messaging stream
C     is accessed from C/C++ and FORTRAN.
C     (2) Names of some routines were changed for integration with rlabplus,
C     because some other earlier integrated solvers with rlab had subroutines
C     with the same name.
      SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,LSTDOUT,MAXFUN,W,
     1  CALFUN,OUTFILE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*)
C     Definition of stream control variables 256 appears in calling codes
      EXTERNAL      CALFUN
      INTEGER       LSTDOUT
      CHARACTER*256 OUTFILE
      IF (LSTDOUT .gt. 1) THEN
          IPRINT = 12
          OPEN (IPRINT, FILE=OUTFILE, STATUS='OLD')
      ELSE
          IPRINT = 0
      ENDIF
C     End of it. KMK

C
C     This subroutine seeks the least value of a function of many variables,
C     by a trust region method that forms quadratic models by interpolation.
C     There can be some freedom in the interpolation conditions, which is
C     taken up by minimizing the Frobenius norm of the change to the second
C     derivative of the quadratic model, beginning with a zero matrix. The
C     arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in the
C       interval [N+2,(N+1)(N+2)/2].
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
C       RHOBEG should be about one tenth of the greatest expected change to a
C       variable, and RHOEND should indicate the accuracy that is required in
C       the final values of the variables.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
C
C     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
C     the value of the objective function for the variables X(1),X(2),...,X(N).
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
      NP=N+1
      NPTM=NPT-NP
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
     1      ' the required interval')
          GO TO 20
      END IF
      NDIM=NPT+N
      IXB=1
      IXO=IXB+N
      IXN=IXO+N
      IXP=IXN+N
      IFV=IXP+N*NPT
      IGQ=IFV+NPT
      IHQ=IGQ+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ID=IZMAT+NPT*NPTM
      IVL=ID+N
      IW=IVL+NDIM
C
C     The above settings provide a partition of W for subroutine NEWUOB.
C     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
C     W plus the space that is needed by the last array of NEWUOB.
C
      CALL NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ),W(IPQ),W(IBMAT),
     2  W(IZMAT),NDIM,W(ID),W(IVL),W(IW),CALFUN)
   20 RETURN
      IF (LSTDOUT .gt. 1)
     * CLOSE (IPRINT)
      END
C
C
C
C
C
      SUBROUTINE NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,XBASE,
     1  XOPT,XNEW,XPT,FVAL,GQ,HQ,PQ,BMAT,ZMAT,NDIM,D,VLAG,W,
     2  CALFUN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XBASE(*),XOPT(*),XNEW(*),XPT(NPT,*),FVAL(*),
     1  GQ(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),VLAG(*),W(*)
      EXTERNAL CALFUN
C
C     The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
C       to the corresponding arguments in SUBROUTINE NEWUOA.
C     XBASE will hold a shift of origin that should reduce the contributions
C       from rounding errors to values of the model and Lagrange functions.
C     XOPT will be set to the displacement from XBASE of the vector of
C       variables that provides the least calculated F so far.
C     XNEW will be set to the displacement from XBASE of the vector of
C       variables for the current calculation of F.
C     XPT will contain the interpolation point coordinates relative to XBASE.
C     FVAL will hold the values of F at the interpolation points.
C     GQ will hold the gradient of the quadratic model at XBASE.
C     HQ will hold the explicit second derivatives of the quadratic model.
C     PQ will contain the parameters of the implicit second derivatives of
C       the quadratic model.
C     BMAT will hold the last N columns of H.
C     ZMAT will hold the factorization of the leading NPT by NPT submatrix of
C       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
C       the elements of DZ are plus or minus one, as specified by IDZ.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     D is reserved for trial steps from XOPT.
C     VLAG will contain the values of the Lagrange functions at a new point X.
C       They are part of a product that requires VLAG to be of length NDIM.
C     The array W will be used for working space. Its length must be at least
C       10*NDIM = 10*(NPT+N).
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TENTH=0.1D0
      ZERO=0.0D0
      NP=N+1
      NH=(N*NP)/2
      NPTM=NPT-NP
      NFTEST=MAX0(MAXFUN,1)
C
C     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
C
      DO 20 J=1,N
      XBASE(J)=X(J)
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 IH=1,NH
   30 HQ(IH)=ZERO
      DO 40 K=1,NPT
      PQ(K)=ZERO
      DO 40 J=1,NPTM
   40 ZMAT(K,J)=ZERO
C
C     Begin the initialization procedure. NF becomes one more than the number
C     of function values so far. The coordinates of the displacement of the
C     next initial interpolation point from XBASE are set in XPT(NF,.).
C
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      RECIQ=DSQRT(HALF)/RHOSQ
      NF=0
   50 NFM=NF
      NFMM=NF-N
      NF=NF+1
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              XPT(NF,NFM)=RHOBEG
          ELSE IF (NFM .GT. N) THEN
              XPT(NF,NFMM)=-RHOBEG
          END IF
      ELSE
          ITEMP=(NFMM-1)/N
          JPT=NFM-ITEMP*N-N
          IPT=JPT+ITEMP
          IF (IPT .GT. N) THEN
              ITEMP=JPT
              JPT=IPT-N
              IPT=ITEMP
          END IF
          XIPT=RHOBEG
          IF (FVAL(IPT+NP) .LT. FVAL(IPT+1)) XIPT=-XIPT
          XJPT=RHOBEG
          IF (FVAL(JPT+NP) .LT. FVAL(JPT+1)) XJPT=-XJPT
          XPT(NF,IPT)=XIPT
          XPT(NF,JPT)=XJPT
      END IF
C
C     Calculate the next value of F, label 70 being reached immediately
C     after this calculation. The least function value so far and its index
C     are required.
C
      DO 60 J=1,N
   60 X(J)=XPT(NF,J)+XBASE(J)
      GOTO 310
   70 FVAL(NF)=F
      IF (NF .EQ. 1) THEN
          FBEG=F
          FOPT=F
          KOPT=1
      ELSE IF (F .LT. FOPT) THEN
          FOPT=F
          KOPT=NF
      END IF
C
C     Set the nonzero initial elements of BMAT and the quadratic model in
C     the cases when NF is at most 2*N+1.
C
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              GQ(NFM)=(F-FBEG)/RHOBEG
              IF (NPT .LT. NF+N) THEN
                  BMAT(1,NFM)=-ONE/RHOBEG
                  BMAT(NF,NFM)=ONE/RHOBEG
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              END IF
          ELSE IF (NFM .GT. N) THEN
              BMAT(NF-N,NFMM)=HALF/RHOBEG
              BMAT(NF,NFMM)=-HALF/RHOBEG
              ZMAT(1,NFMM)=-RECIQ-RECIQ
              ZMAT(NF-N,NFMM)=RECIQ
              ZMAT(NF,NFMM)=RECIQ
              IH=(NFMM*(NFMM+1))/2
              TEMP=(FBEG-F)/RHOBEG
              HQ(IH)=(GQ(NFMM)-TEMP)/RHOBEG
              GQ(NFMM)=HALF*(GQ(NFMM)+TEMP)
          END IF
C
C     Set the off-diagonal second derivatives of the Lagrange functions and
C     the initial quadratic model.
C
      ELSE
          IH=(IPT*(IPT-1))/2+JPT
          IF (XIPT .LT. ZERO) IPT=IPT+N
          IF (XJPT .LT. ZERO) JPT=JPT+N
          ZMAT(1,NFMM)=RECIP
          ZMAT(NF,NFMM)=RECIP
          ZMAT(IPT+1,NFMM)=-RECIP
          ZMAT(JPT+1,NFMM)=-RECIP
          HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/(XIPT*XJPT)
      END IF
      IF (NF .LT. NPT) GOTO 50
C
C     Begin the iterative procedure, because the initial model is complete.
C
      RHO=RHOBEG
      DELTA=RHO
      IDZ=1
      DIFFA=ZERO
      DIFFB=ZERO
      ITEST=0
      XOPTSQ=ZERO
      DO 80 I=1,N
      XOPT(I)=XPT(KOPT,I)
   80 XOPTSQ=XOPTSQ+XOPT(I)**2
   90 NFSAV=NF
C
C     Generate the next trust region step and test its length. Set KNEW
C     to -1 if the purpose of the next F will be to improve the model.
C
  100 KNEW=0
      CALL TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,D,W,W(NP),
     1  W(NP+N),W(NP+2*N),CRVMIN)
      DSQ=ZERO
      DO 110 I=1,N
  110 DSQ=DSQ+D(I)**2
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      IF (DNORM .LT. HALF*RHO) THEN
          KNEW=-1
          DELTA=TENTH*DELTA
          RATIO=-1.0D0
          IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
          IF (NF .LE. NFSAV+2) GOTO 460
          TEMP=0.125D0*CRVMIN*RHO*RHO
          IF (TEMP .LE. DMAX1(DIFFA,DIFFB,DIFFC)) GOTO 460
          GOTO 490
      END IF
C
C     Shift XBASE if XOPT may be too far from XBASE. First make the changes
C     to BMAT that do not depend on ZMAT.
C
  120 IF (DSQ .LE. 1.0D-3*XOPTSQ) THEN
          TEMPQ=0.25D0*XOPTSQ
          DO 140 K=1,NPT
          SUM=ZERO
          DO 130 I=1,N
  130     SUM=SUM+XPT(K,I)*XOPT(I)
          TEMP=PQ(K)*SUM
          SUM=SUM-HALF*XOPTSQ
          W(NPT+K)=SUM
          DO 140 I=1,N
          GQ(I)=GQ(I)+TEMP*XPT(K,I)
          XPT(K,I)=XPT(K,I)-HALF*XOPT(I)
          VLAG(I)=BMAT(K,I)
          W(I)=SUM*XPT(K,I)+TEMPQ*XOPT(I)
          IP=NPT+I
          DO 140 J=1,I
  140     BMAT(IP,J)=BMAT(IP,J)+VLAG(I)*W(J)+W(I)*VLAG(J)
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO 180 K=1,NPTM
          SUMZ=ZERO
          DO 150 I=1,NPT
          SUMZ=SUMZ+ZMAT(I,K)
  150     W(I)=W(NPT+I)*ZMAT(I,K)
          DO 170 J=1,N
          SUM=TEMPQ*SUMZ*XOPT(J)
          DO 160 I=1,NPT
  160     SUM=SUM+W(I)*XPT(I,J)
          VLAG(J)=SUM
          IF (K .LT. IDZ) SUM=-SUM
          DO 170 I=1,NPT
  170     BMAT(I,J)=BMAT(I,J)+SUM*ZMAT(I,K)
          DO 180 I=1,N
          IP=I+NPT
          TEMP=VLAG(I)
          IF (K .LT. IDZ) TEMP=-TEMP
          DO 180 J=1,I
  180     BMAT(IP,J)=BMAT(IP,J)+TEMP*VLAG(J)
C
C     The following instructions complete the shift of XBASE, including
C     the changes to the parameters of the quadratic model.
C
          IH=0
          DO 200 J=1,N
          W(J)=ZERO
          DO 190 K=1,NPT
          W(J)=W(J)+PQ(K)*XPT(K,J)
  190     XPT(K,J)=XPT(K,J)-HALF*XOPT(J)
          DO 200 I=1,J
          IH=IH+1
          IF (I .LT. J) GQ(J)=GQ(J)+HQ(IH)*XOPT(I)
          GQ(I)=GQ(I)+HQ(IH)*XOPT(J)
          HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
  200     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 210 J=1,N
          XBASE(J)=XBASE(J)+XOPT(J)
  210     XOPT(J)=ZERO
          XOPTSQ=ZERO
      END IF
C
C     Pick the model step if KNEW is positive. A different choice of D
C     may be made later, if the choice of D by BIGLAG causes substantial
C     cancellation in DENOM.
C
      IF (KNEW .GT. 0) THEN
          CALL BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW,DSTEP,
     1      D,ALPHA,VLAG,VLAG(NPT+1),W,W(NP),W(NP+N))
      END IF
C
C     Calculate VLAG and BETA for the current choice of D. The first NPT
C     components of W_check will be held in W.
C
      DO 230 K=1,NPT
      SUMA=ZERO
      SUMB=ZERO
      SUM=ZERO
      DO 220 J=1,N
      SUMA=SUMA+XPT(K,J)*D(J)
      SUMB=SUMB+XPT(K,J)*XOPT(J)
  220 SUM=SUM+BMAT(K,J)*D(J)
      W(K)=SUMA*(HALF*SUMA+SUMB)
  230 VLAG(K)=SUM
      BETA=ZERO
      DO 250 K=1,NPTM
      SUM=ZERO
      DO 240 I=1,NPT
  240 SUM=SUM+ZMAT(I,K)*W(I)
      IF (K .LT. IDZ) THEN
          BETA=BETA+SUM*SUM
          SUM=-SUM
      ELSE
          BETA=BETA-SUM*SUM
      END IF
      DO 250 I=1,NPT
  250 VLAG(I)=VLAG(I)+SUM*ZMAT(I,K)
      BSUM=ZERO
      DX=ZERO
      DO 280 J=1,N
      SUM=ZERO
      DO 260 I=1,NPT
  260 SUM=SUM+W(I)*BMAT(I,J)
      BSUM=BSUM+SUM*D(J)
      JP=NPT+J
      DO 270 K=1,N
  270 SUM=SUM+BMAT(JP,K)*D(K)
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D(J)
  280 DX=DX+D(J)*XOPT(J)
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     If KNEW is positive and if the cancellation in DENOM is unacceptable,
C     then BIGDEN calculates an alternative model step, XNEW being used for
C     working space.
C
      IF (KNEW .GT. 0) THEN
          TEMP=ONE+ALPHA*BETA/VLAG(KNEW)**2
          IF (DABS(TEMP) .LE. 0.8D0) THEN
              CALL BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1          KNEW,D,W,VLAG,BETA,XNEW,W(NDIM+1),W(6*NDIM+1))
          END IF
      END IF
C
C     Calculate the next value of the objective function.
C
  290 DO 300 I=1,N
      XNEW(I)=XOPT(I)+D(I)
  300 X(I)=XBASE(I)+XNEW(I)
      NF=NF+1
  310 IF (NF .GT. NFTEST) THEN
          NF=NF-1
          IF (IPRINT .GT. 0) THEN
            WRITE(IPRINT,320)
  320       FORMAT (/4X,'Return from NEWUOA because CALFUN has been',
     1        ' called MAXFUN times.')
          ENDIF
          GOTO 530
      END IF
      CALL CALFUN (N,X,F)
      IF (IPRINT .GE. 3) THEN
         WRITE(IPRINT,330) NF,F,(X(I),I=1,N)
  330      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NF .LE. NPT) GOTO 70
      IF (KNEW .EQ. -1) GOTO 530
C
C     Use the quadratic model to predict the change in F due to the step D,
C     and set DIFF to the error of this prediction.
C
      VQUAD=ZERO
      IH=0
      DO 340 J=1,N
      VQUAD=VQUAD+D(J)*GQ(J)
      DO 340 I=1,J
      IH=IH+1
      TEMP=D(I)*XNEW(J)+D(J)*XOPT(I)
      IF (I .EQ. J) TEMP=HALF*TEMP
  340 VQUAD=VQUAD+TEMP*HQ(IH)
      DO 350 K=1,NPT
  350 VQUAD=VQUAD+PQ(K)*W(K)
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM .GT. RHO) NFSAV=NF
C
C     Update FOPT and XOPT if the new F is the least value of the objective
C     function so far. The branch when KNEW is positive occurs if D is not
C     a trust region step.
C
      FSAVE=FOPT
      IF (F .LT. FOPT) THEN
          FOPT=F
          XOPTSQ=ZERO
          DO 360 I=1,N
          XOPT(I)=XNEW(I)
  360     XOPTSQ=XOPTSQ+XOPT(I)**2
      END IF
      KSAVE=KNEW
      IF (KNEW .GT. 0) GOTO 410
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (VQUAD .GE. ZERO) THEN
          IF (IPRINT .GT. 0) THEN
            WRITE(IPRINT,370)
  370       FORMAT (/4X,'Return from NEWUOA because a trust',
     1        ' region step has failed to reduce Q.')
          ENDIF
          GOTO 530
      END IF
      RATIO=(F-FSAVE)/VQUAD
      IF (RATIO .LE. TENTH) THEN
          DELTA=HALF*DNORM
      ELSE IF (RATIO. LE. 0.7D0) THEN
          DELTA=DMAX1(HALF*DELTA,DNORM)
      ELSE
          DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
      END IF
      IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
C
C     Set KNEW to the index of the next interpolation point to be deleted.
C
      RHOSQ=DMAX1(TENTH*DELTA,RHO)**2
      KTEMP=0
      DETRAT=ZERO
      IF (F .GE. FSAVE) THEN
          KTEMP=KOPT
          DETRAT=ONE
      END IF
      DO 400 K=1,NPT
      HDIAG=ZERO
      DO 380 J=1,NPTM
      TEMP=ONE
      IF (J .LT. IDZ) TEMP=-ONE
  380 HDIAG=HDIAG+TEMP*ZMAT(K,J)**2
      TEMP=DABS(BETA*HDIAG+VLAG(K)**2)
      DISTSQ=ZERO
      DO 390 J=1,N
  390 DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
      IF (DISTSQ .GT. RHOSQ) TEMP=TEMP*(DISTSQ/RHOSQ)**3
      IF (TEMP .GT. DETRAT .AND. K .NE. KTEMP) THEN
          DETRAT=TEMP
          KNEW=K
      END IF
  400 CONTINUE
      IF (KNEW .EQ. 0) GOTO 460
C
C     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
C     can be moved. Begin the updating of the quadratic model, starting
C     with the explicit second derivative term.
C
  410 CALL NEWUOUPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
      FVAL(KNEW)=F
      IH=0
      DO 420 I=1,N
      TEMP=PQ(KNEW)*XPT(KNEW,I)
      DO 420 J=1,I
      IH=IH+1
  420 HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
      PQ(KNEW)=ZERO
C
C     Update the other second derivative parameters, and then the gradient
C     vector of the model. Also include the new interpolation point.
C
      DO 440 J=1,NPTM
      TEMP=DIFF*ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 440 K=1,NPT
  440 PQ(K)=PQ(K)+TEMP*ZMAT(K,J)
      GQSQ=ZERO
      DO 450 I=1,N
      GQ(I)=GQ(I)+DIFF*BMAT(KNEW,I)
      GQSQ=GQSQ+GQ(I)**2
  450 XPT(KNEW,I)=XNEW(I)
C
C     If a trust region step makes a small change to the objective function,
C     then calculate the gradient of the least Frobenius norm interpolant at
C     XBASE, and store it in W, using VLAG for a vector of right hand sides.
C
      IF (KSAVE .EQ. 0 .AND. DELTA .EQ. RHO) THEN
          IF (DABS(RATIO) .GT. 1.0D-2) THEN
              ITEST=0
          ELSE
              DO 700 K=1,NPT
  700         VLAG(K)=FVAL(K)-FVAL(KOPT)
              GISQ=ZERO
              DO 720 I=1,N
              SUM=ZERO
              DO 710 K=1,NPT
  710         SUM=SUM+BMAT(K,I)*VLAG(K)
              GISQ=GISQ+SUM*SUM
  720         W(I)=SUM
C
C     Test whether to replace the new quadratic model by the least Frobenius
C     norm interpolant, making the replacement if the test is satisfied.
C
              ITEST=ITEST+1
              IF (GQSQ .LT. 1.0D2*GISQ) ITEST=0
              IF (ITEST .GE. 3) THEN
                  DO 730 I=1,N
  730             GQ(I)=W(I)
                  DO 740 IH=1,NH
  740             HQ(IH)=ZERO
                  DO 760 J=1,NPTM
                  W(J)=ZERO
                  DO 750 K=1,NPT
  750             W(J)=W(J)+VLAG(K)*ZMAT(K,J)
  760             IF (J .LT. IDZ) W(J)=-W(J)
                  DO 770 K=1,NPT
                  PQ(K)=ZERO
                  DO 770 J=1,NPTM
  770             PQ(K)=PQ(K)+ZMAT(K,J)*W(J)
                  ITEST=0
              END IF
          END IF
      END IF
      IF (F .LT. FSAVE) KOPT=KNEW
C
C     If a trust region step has provided a sufficient decrease in F, then
C     branch for another trust region calculation. The case KSAVE>0 occurs
C     when the new function value was calculated by a model step.
C
      IF (F .LE. FSAVE+TENTH*VQUAD) GOTO 100
      IF (KSAVE .GT. 0) GOTO 100
C
C     Alternatively, find out if the interpolation points are close enough
C     to the best point so far.
C
      KNEW=0
  460 DISTSQ=4.0D0*DELTA*DELTA
      DO 480 K=1,NPT
      SUM=ZERO
      DO 470 J=1,N
  470 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM .GT. DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  480 CONTINUE
C
C     If KNEW is positive, then set DSTEP, and branch back for the next
C     iteration, which will generate a "model step".
C
      IF (KNEW .GT. 0) THEN
          DSTEP=DMAX1(DMIN1(TENTH*DSQRT(DISTSQ),HALF*DELTA),RHO)
          DSQ=DSTEP*DSTEP
          GOTO 120
      END IF
      IF (RATIO .GT. ZERO) GOTO 100
      IF (DMAX1(DELTA,DNORM) .GT. RHO) GOTO 100
C
C     The calculations with the current value of RHO are complete. Pick the
C     next values of RHO and DELTA.
C
  490 IF (RHO .GT. RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO .LE. 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO .LE. 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) THEN
                WRITE (IPRINT,500)
  500           FORMAT (5X)
              ENDIF
              WRITE (IPRINT,510) RHO,NF
  510         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              WRITE(IPRINT,520) FOPT,(XBASE(I)+XOPT(I),I=1,N)
  520         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          GOTO 90
      END IF
C
C     Return from the calculation, after another Newton-Raphson step, if
C     it is too short to have been tried before.
C
      IF (KNEW .EQ. -1) GOTO 290
  530 IF (FOPT .LE. F) THEN
          DO 540 I=1,N
  540     X(I)=XBASE(I)+XOPT(I)
          F=FOPT
      END IF
      IF (IPRINT .GE. 1) THEN
          WRITE (IPRINT,550) NF
  550     FORMAT (/4X,'At the return from NEWUOA',5X,
     1      'Number of function values =',I6)
          WRITE (IPRINT,520) F,(X(I),I=1,N)
      END IF
      RETURN
      END
C
C
C
C
C
      SUBROUTINE BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT,
     1  KNEW,D,W,VLAG,BETA,S,WVEC,PROD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),
     1  W(*),VLAG(*),S(*),WVEC(NDIM,*),PROD(NDIM,*)
      DIMENSION DEN(9),DENEX(9),PAR(9)
C
C     N is the number of variables.
C     NPT is the number of interpolation equations.
C     XOPT is the best interpolation point so far.
C     XPT contains the coordinates of the current interpolation points.
C     BMAT provides the last N columns of H.
C     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     KOPT is the index of the optimal interpolation point.
C     KNEW is the index of the interpolation point that is going to be moved.
C     D will be set to the step from XOPT to the new point, and on entry it
C       should be the D that was calculated by the last call of BIGLAG. The
C       length of the initial D provides a trust region bound on the final D.
C     W will be set to Wcheck for the final choice of D.
C     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
C     BETA will be set to the value that will occur in the updating formula
C       when the KNEW-th interpolation point is moved to its new position.
C     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
C       for working space.
C
C     D is calculated in a way that should provide a denominator with a large
C     modulus in the updating formula when the KNEW-th interpolation point is
C     shifted to the new position XOPT+D.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      QUART=0.25D0
      TWO=2.0D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(ONE)
      NPTM=NPT-N-1
C
C     Store the first NPT elements of the KNEW-th column of H in W(N+1)
C     to W(N+NPT).
C
      DO 10 K=1,NPT
   10 W(N+K)=ZERO
      DO 20 J=1,NPTM
      TEMP=ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 20 K=1,NPT
   20 W(N+K)=W(N+K)+TEMP*ZMAT(K,J)
      ALPHA=W(N+KNEW)
C
C     The initial search direction D is taken from the last call of BIGLAG,
C     and the initial S is set below, usually to the direction from X_OPT
C     to X_KNEW, but a different direction to an interpolation point may
C     be chosen, in order to prevent S from being nearly parallel to D.
C
      DD=ZERO
      DS=ZERO
      SS=ZERO
      XOPTSQ=ZERO
      DO 30 I=1,N
      DD=DD+D(I)**2
      S(I)=XPT(KNEW,I)-XOPT(I)
      DS=DS+D(I)*S(I)
      SS=SS+S(I)**2
   30 XOPTSQ=XOPTSQ+XOPT(I)**2
      IF (DS*DS .GT. 0.99D0*DD*SS) THEN
          KSAV=KNEW
          DTEST=DS*DS/SS
          DO 50 K=1,NPT
          IF (K .NE. KOPT) THEN
              DSTEMP=ZERO
              SSTEMP=ZERO
              DO 40 I=1,N
              DIFF=XPT(K,I)-XOPT(I)
              DSTEMP=DSTEMP+D(I)*DIFF
   40         SSTEMP=SSTEMP+DIFF*DIFF
              IF (DSTEMP*DSTEMP/SSTEMP .LT. DTEST) THEN
                  KSAV=K
                  DTEST=DSTEMP*DSTEMP/SSTEMP
                  DS=DSTEMP
                  SS=SSTEMP
              END IF
          END IF
   50     CONTINUE
          DO 60 I=1,N
   60     S(I)=XPT(KSAV,I)-XOPT(I)
      END IF
      SSDEN=DD*SS-DS*DS
      ITERC=0
      DENSAV=ZERO
C
C     Begin the iteration by overwriting S with a vector that has the
C     required length and direction.
C
   70 ITERC=ITERC+1
      TEMP=ONE/DSQRT(SSDEN)
      XOPTD=ZERO
      XOPTS=ZERO
      DO 80 I=1,N
      S(I)=TEMP*(DD*S(I)-DS*D(I))
      XOPTD=XOPTD+XOPT(I)*D(I)
   80 XOPTS=XOPTS+XOPT(I)*S(I)
C
C     Set the coefficients of the first two terms of BETA.
C
      TEMPA=HALF*XOPTD*XOPTD
      TEMPB=HALF*XOPTS*XOPTS
      DEN(1)=DD*(XOPTSQ+HALF*DD)+TEMPA+TEMPB
      DEN(2)=TWO*XOPTD*DD
      DEN(3)=TWO*XOPTS*DD
      DEN(4)=TEMPA-TEMPB
      DEN(5)=XOPTD*XOPTS
      DO 90 I=6,9
   90 DEN(I)=ZERO
C
C     Put the coefficients of Wcheck in WVEC.
C
      DO 110 K=1,NPT
      TEMPA=ZERO
      TEMPB=ZERO
      TEMPC=ZERO
      DO 100 I=1,N
      TEMPA=TEMPA+XPT(K,I)*D(I)
      TEMPB=TEMPB+XPT(K,I)*S(I)
  100 TEMPC=TEMPC+XPT(K,I)*XOPT(I)
      WVEC(K,1)=QUART*(TEMPA*TEMPA+TEMPB*TEMPB)
      WVEC(K,2)=TEMPA*TEMPC
      WVEC(K,3)=TEMPB*TEMPC
      WVEC(K,4)=QUART*(TEMPA*TEMPA-TEMPB*TEMPB)
  110 WVEC(K,5)=HALF*TEMPA*TEMPB
      DO 120 I=1,N
      IP=I+NPT
      WVEC(IP,1)=ZERO
      WVEC(IP,2)=D(I)
      WVEC(IP,3)=S(I)
      WVEC(IP,4)=ZERO
  120 WVEC(IP,5)=ZERO
C
C     Put the coefficents of THETA*Wcheck in PROD.
C
      DO 190 JC=1,5
      NW=NPT
      IF (JC .EQ. 2 .OR. JC .EQ. 3) NW=NDIM
      DO 130 K=1,NPT
  130 PROD(K,JC)=ZERO
      DO 150 J=1,NPTM
      SUM=ZERO
      DO 140 K=1,NPT
  140 SUM=SUM+ZMAT(K,J)*WVEC(K,JC)
      IF (J .LT. IDZ) SUM=-SUM
      DO 150 K=1,NPT
  150 PROD(K,JC)=PROD(K,JC)+SUM*ZMAT(K,J)
      IF (NW .EQ. NDIM) THEN
          DO 170 K=1,NPT
          SUM=ZERO
          DO 160 J=1,N
  160     SUM=SUM+BMAT(K,J)*WVEC(NPT+J,JC)
  170     PROD(K,JC)=PROD(K,JC)+SUM
      END IF
      DO 190 J=1,N
      SUM=ZERO
      DO 180 I=1,NW
  180 SUM=SUM+BMAT(I,J)*WVEC(I,JC)
  190 PROD(NPT+J,JC)=SUM
C
C     Include in DEN the part of BETA that depends on THETA.
C
      DO 210 K=1,NDIM
      SUM=ZERO
      DO 200 I=1,5
      PAR(I)=HALF*PROD(K,I)*WVEC(K,I)
  200 SUM=SUM+PAR(I)
      DEN(1)=DEN(1)-PAR(1)-SUM
      TEMPA=PROD(K,1)*WVEC(K,2)+PROD(K,2)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,4)+PROD(K,4)*WVEC(K,2)
      TEMPC=PROD(K,3)*WVEC(K,5)+PROD(K,5)*WVEC(K,3)
      DEN(2)=DEN(2)-TEMPA-HALF*(TEMPB+TEMPC)
      DEN(6)=DEN(6)-HALF*(TEMPB-TEMPC)
      TEMPA=PROD(K,1)*WVEC(K,3)+PROD(K,3)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,5)+PROD(K,5)*WVEC(K,2)
      TEMPC=PROD(K,3)*WVEC(K,4)+PROD(K,4)*WVEC(K,3)
      DEN(3)=DEN(3)-TEMPA-HALF*(TEMPB-TEMPC)
      DEN(7)=DEN(7)-HALF*(TEMPB+TEMPC)
      TEMPA=PROD(K,1)*WVEC(K,4)+PROD(K,4)*WVEC(K,1)
      DEN(4)=DEN(4)-TEMPA-PAR(2)+PAR(3)
      TEMPA=PROD(K,1)*WVEC(K,5)+PROD(K,5)*WVEC(K,1)
      TEMPB=PROD(K,2)*WVEC(K,3)+PROD(K,3)*WVEC(K,2)
      DEN(5)=DEN(5)-TEMPA-HALF*TEMPB
      DEN(8)=DEN(8)-PAR(4)+PAR(5)
      TEMPA=PROD(K,4)*WVEC(K,5)+PROD(K,5)*WVEC(K,4)
  210 DEN(9)=DEN(9)-HALF*TEMPA
C
C     Extend DEN so that it holds all the coefficients of DENOM.
C
      SUM=ZERO
      DO 220 I=1,5
      PAR(I)=HALF*PROD(KNEW,I)**2
  220 SUM=SUM+PAR(I)
      DENEX(1)=ALPHA*DEN(1)+PAR(1)+SUM
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,2)
      TEMPB=PROD(KNEW,2)*PROD(KNEW,4)
      TEMPC=PROD(KNEW,3)*PROD(KNEW,5)
      DENEX(2)=ALPHA*DEN(2)+TEMPA+TEMPB+TEMPC
      DENEX(6)=ALPHA*DEN(6)+TEMPB-TEMPC
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,3)
      TEMPB=PROD(KNEW,2)*PROD(KNEW,5)
      TEMPC=PROD(KNEW,3)*PROD(KNEW,4)
      DENEX(3)=ALPHA*DEN(3)+TEMPA+TEMPB-TEMPC
      DENEX(7)=ALPHA*DEN(7)+TEMPB+TEMPC
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,4)
      DENEX(4)=ALPHA*DEN(4)+TEMPA+PAR(2)-PAR(3)
      TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,5)
      DENEX(5)=ALPHA*DEN(5)+TEMPA+PROD(KNEW,2)*PROD(KNEW,3)
      DENEX(8)=ALPHA*DEN(8)+PAR(4)-PAR(5)
      DENEX(9)=ALPHA*DEN(9)+PROD(KNEW,4)*PROD(KNEW,5)
C
C     Seek the value of the angle that maximizes the modulus of DENOM.
C
      SUM=DENEX(1)+DENEX(2)+DENEX(4)+DENEX(6)+DENEX(8)
      DENOLD=SUM
      DENMAX=SUM
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      PAR(1)=ONE
      DO 250 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      PAR(2)=DCOS(ANGLE)
      PAR(3)=DSIN(ANGLE)
      DO 230 J=4,8,2
      PAR(J)=PAR(2)*PAR(J-2)-PAR(3)*PAR(J-1)
  230 PAR(J+1)=PAR(2)*PAR(J-1)+PAR(3)*PAR(J-2)
      SUMOLD=SUM
      SUM=ZERO
      DO 240 J=1,9
  240 SUM=SUM+DENEX(J)*PAR(J)
      IF (DABS(SUM) .GT. DABS(DENMAX)) THEN
          DENMAX=SUM
          ISAVE=I
          TEMPA=SUMOLD
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=SUM
      END IF
  250 CONTINUE
      IF (ISAVE .EQ. 0) TEMPA=SUM
      IF (ISAVE .EQ. IU) TEMPB=DENOLD
      STEP=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-DENMAX
          TEMPB=TEMPB-DENMAX
          STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+STEP)
C
C     Calculate the new parameters of the denominator, the new VLAG vector
C     and the new D. Then test for convergence.
C
      PAR(2)=DCOS(ANGLE)
      PAR(3)=DSIN(ANGLE)
      DO 260 J=4,8,2
      PAR(J)=PAR(2)*PAR(J-2)-PAR(3)*PAR(J-1)
  260 PAR(J+1)=PAR(2)*PAR(J-1)+PAR(3)*PAR(J-2)
      BETA=ZERO
      DENMAX=ZERO
      DO 270 J=1,9
      BETA=BETA+DEN(J)*PAR(J)
  270 DENMAX=DENMAX+DENEX(J)*PAR(J)
      DO 280 K=1,NDIM
      VLAG(K)=ZERO
      DO 280 J=1,5
  280 VLAG(K)=VLAG(K)+PROD(K,J)*PAR(J)
      TAU=VLAG(KNEW)
      DD=ZERO
      TEMPA=ZERO
      TEMPB=ZERO
      DO 290 I=1,N
      D(I)=PAR(2)*D(I)+PAR(3)*S(I)
      W(I)=XOPT(I)+D(I)
      DD=DD+D(I)**2
      TEMPA=TEMPA+D(I)*W(I)
  290 TEMPB=TEMPB+W(I)*W(I)
      IF (ITERC .GE. N) GOTO 340
      IF (ITERC .GT. 1) DENSAV=DMAX1(DENSAV,DENOLD)
      IF (DABS(DENMAX) .LE. 1.1D0*DABS(DENSAV)) GOTO 340
      DENSAV=DENMAX
C
C     Set S to half the gradient of the denominator with respect to D.
C     Then branch for the next iteration.
C
      DO 300 I=1,N
      TEMP=TEMPA*XOPT(I)+TEMPB*D(I)-VLAG(NPT+I)
  300 S(I)=TAU*BMAT(KNEW,I)+ALPHA*TEMP
      DO 320 K=1,NPT
      SUM=ZERO
      DO 310 J=1,N
  310 SUM=SUM+XPT(K,J)*W(J)
      TEMP=(TAU*W(N+K)-ALPHA*VLAG(K))*SUM
      DO 320 I=1,N
  320 S(I)=S(I)+TEMP*XPT(K,I)
      SS=ZERO
      DS=ZERO
      DO 330 I=1,N
      SS=SS+S(I)**2
  330 DS=DS+D(I)*S(I)
      SSDEN=DD*SS-DS*DS
      IF (SSDEN .GE. 1.0D-8*DD*SS) GOTO 70
C
C     Set the vector W before the RETURN from the subroutine.
C
  340 DO 350 K=1,NDIM
      W(K)=ZERO
      DO 350 J=1,5
  350 W(K)=W(K)+WVEC(K,J)*PAR(J)
      VLAG(KOPT)=VLAG(KOPT)+ONE
      RETURN
      END
C
C
C
C
C
      SUBROUTINE BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW,
     1  DELTA,D,ALPHA,HCOL,GC,GD,S,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*),
     1  HCOL(*),GC(*),GD(*),S(*),W(*)
C
C     N is the number of variables.
C     NPT is the number of interpolation equations.
C     XOPT is the best interpolation point so far.
C     XPT contains the coordinates of the current interpolation points.
C     BMAT provides the last N columns of H.
C     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     KNEW is the index of the interpolation point that is going to be moved.
C     DELTA is the current trust region bound.
C     D will be set to the step from XOPT to the new point.
C     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
C     HCOL, GC, GD, S and W will be used for working space.
C
C     The step D is calculated in a way that attempts to maximize the modulus
C     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
C     the KNEW-th Lagrange function.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(ONE)
      DELSQ=DELTA*DELTA
      NPTM=NPT-N-1
C
C     Set the first NPT components of HCOL to the leading elements of the
C     KNEW-th column of H.
C
      ITERC=0
      DO 10 K=1,NPT
   10 HCOL(K)=ZERO
      DO 20 J=1,NPTM
      TEMP=ZMAT(KNEW,J)
      IF (J .LT. IDZ) TEMP=-TEMP
      DO 20 K=1,NPT
   20 HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
      ALPHA=HCOL(KNEW)
C
C     Set the unscaled initial direction D. Form the gradient of LFUNC at
C     XOPT, and multiply D by the second derivative matrix of LFUNC.
C
      DD=ZERO
      DO 30 I=1,N
      D(I)=XPT(KNEW,I)-XOPT(I)
      GC(I)=BMAT(KNEW,I)
      GD(I)=ZERO
   30 DD=DD+D(I)**2
      DO 50 K=1,NPT
      TEMP=ZERO
      SUM=ZERO
      DO 40 J=1,N
      TEMP=TEMP+XPT(K,J)*XOPT(J)
   40 SUM=SUM+XPT(K,J)*D(J)
      TEMP=HCOL(K)*TEMP
      SUM=HCOL(K)*SUM
      DO 50 I=1,N
      GC(I)=GC(I)+TEMP*XPT(K,I)
   50 GD(I)=GD(I)+SUM*XPT(K,I)
C
C     Scale D and GD, with a sign change if required. Set S to another
C     vector in the initial two dimensional subspace.
C
      GG=ZERO
      SP=ZERO
      DHD=ZERO
      DO 60 I=1,N
      GG=GG+GC(I)**2
      SP=SP+D(I)*GC(I)
   60 DHD=DHD+D(I)*GD(I)
      SCALE=DELTA/DSQRT(DD)
      IF (SP*DHD .LT. ZERO) SCALE=-SCALE
      TEMP=ZERO
      IF (SP*SP .GT. 0.99D0*DD*GG) TEMP=ONE
      TAU=SCALE*(DABS(SP)+HALF*SCALE*DABS(DHD))
      IF (GG*DELSQ .LT. 0.01D0*TAU*TAU) TEMP=ONE
      DO 70 I=1,N
      D(I)=SCALE*D(I)
      GD(I)=SCALE*GD(I)
   70 S(I)=GC(I)+TEMP*GD(I)
C
C     Begin the iteration by overwriting S with a vector that has the
C     required length and direction, except that termination occurs if
C     the given D and S are nearly parallel.
C
   80 ITERC=ITERC+1
      DD=ZERO
      SP=ZERO
      SS=ZERO
      DO 90 I=1,N
      DD=DD+D(I)**2
      SP=SP+D(I)*S(I)
   90 SS=SS+S(I)**2
      TEMP=DD*SS-SP*SP
      IF (TEMP .LE. 1.0D-8*DD*SS) GOTO 160
      DENOM=DSQRT(TEMP)
      DO 100 I=1,N
      S(I)=(DD*S(I)-SP*D(I))/DENOM
  100 W(I)=ZERO
C
C     Calculate the coefficients of the objective function on the circle,
C     beginning with the multiplication of S by the second derivative matrix.
C
      DO 120 K=1,NPT
      SUM=ZERO
      DO 110 J=1,N
  110 SUM=SUM+XPT(K,J)*S(J)
      SUM=HCOL(K)*SUM
      DO 120 I=1,N
  120 W(I)=W(I)+SUM*XPT(K,I)
      CF1=ZERO
      CF2=ZERO
      CF3=ZERO
      CF4=ZERO
      CF5=ZERO
      DO 130 I=1,N
      CF1=CF1+S(I)*W(I)
      CF2=CF2+D(I)*GC(I)
      CF3=CF3+S(I)*GC(I)
      CF4=CF4+D(I)*GD(I)
  130 CF5=CF5+S(I)*GD(I)
      CF1=HALF*CF1
      CF4=HALF*CF4-CF1
C
C     Seek the value of the angle that maximizes the modulus of TAU.
C
      TAUBEG=CF1+CF2+CF4
      TAUMAX=TAUBEG
      TAUOLD=TAUBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO 140 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
      IF (DABS(TAU) .GT. DABS(TAUMAX)) THEN
          TAUMAX=TAU
          ISAVE=I
          TEMPA=TAUOLD
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=TAU
      END IF
  140 TAUOLD=TAU
      IF (ISAVE .EQ. 0) TEMPA=TAU
      IF (ISAVE .EQ. IU) TEMPB=TAUBEG
      STEP=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-TAUMAX
          TEMPB=TEMPB-TAUMAX
          STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+STEP)
C
C     Calculate the new D and GD. Then test for convergence.
C
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
      DO 150 I=1,N
      D(I)=CTH*D(I)+STH*S(I)
      GD(I)=CTH*GD(I)+STH*W(I)
  150 S(I)=GC(I)+GD(I)
      IF (DABS(TAU) .LE. 1.1D0*DABS(TAUBEG)) GOTO 160
      IF (ITERC .LT. N) GOTO 80
  160 RETURN
      END
C
C
C
C
C
      SUBROUTINE TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,STEP,
     1  D,G,HD,HS,CRVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XOPT(*),XPT(NPT,*),GQ(*),HQ(*),PQ(*),STEP(*),
     1  D(*),G(*),HD(*),HS(*)
C
C     N is the number of variables of a quadratic objective function, Q say.
C     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
C       in order to define the current quadratic model Q.
C     DELTA is the trust region radius, and has to be positive.
C     STEP will be set to the calculated trial step.
C     The arrays D, G, HD and HS will be used for working space.
C     CRVMIN will be set to the least curvature of H along the conjugate
C       directions that occur, except that it is set to zero if STEP goes
C       all the way to the trust region boundary.
C
C     The calculation of STEP begins with the truncated conjugate gradient
C     method. If the boundary of the trust region is reached, then further
C     changes to STEP may be made, each one being in the 2D space spanned
C     by the current STEP and the corresponding gradient of Q. Thus STEP
C     should provide a substantial reduction to Q within the trust region.
C
C     Initialization, which includes setting HD to H times XOPT.
C
      HALF=0.5D0
      ZERO=0.0D0
      TWOPI=8.0D0*DATAN(1.0D0)
      DELSQ=DELTA*DELTA
      ITERC=0
      ITERMAX=N
      ITERSW=ITERMAX
      DO 10 I=1,N
   10 D(I)=XOPT(I)
      GOTO 170
C
C     Prepare for the first line search.
C
   20 QRED=ZERO
      DD=ZERO
      DO 30 I=1,N
      STEP(I)=ZERO
      HS(I)=ZERO
      G(I)=GQ(I)+HD(I)
      D(I)=-G(I)
   30 DD=DD+D(I)**2
      CRVMIN=ZERO
      IF (DD .EQ. ZERO) GOTO 160
      DS=ZERO
      SS=ZERO
      GG=DD
      GGBEG=GG
C
C     Calculate the step to the trust region boundary and the product HD.
C
   40 ITERC=ITERC+1
      TEMP=DELSQ-SS
      BSTEP=TEMP/(DS+DSQRT(DS*DS+DD*TEMP))
      GOTO 170
   50 DHD=ZERO
      DO 60 J=1,N
   60 DHD=DHD+D(J)*HD(J)
C
C     Update CRVMIN and set the step-length ALPHA.
C
      ALPHA=BSTEP
      IF (DHD .GT. ZERO) THEN
          TEMP=DHD/DD
          IF (ITERC .EQ. 1) CRVMIN=TEMP
          CRVMIN=DMIN1(CRVMIN,TEMP)
          ALPHA=DMIN1(ALPHA,GG/DHD)
      END IF
      QADD=ALPHA*(GG-HALF*ALPHA*DHD)
      QRED=QRED+QADD
C
C     Update STEP and HS.
C
      GGSAV=GG
      GG=ZERO
      DO 70 I=1,N
      STEP(I)=STEP(I)+ALPHA*D(I)
      HS(I)=HS(I)+ALPHA*HD(I)
   70 GG=GG+(G(I)+HS(I))**2
C
C     Begin another conjugate direction iteration if required.
C
      IF (ALPHA .LT. BSTEP) THEN
          IF (QADD .LE. 0.01D0*QRED) GOTO 160
          IF (GG .LE. 1.0D-4*GGBEG) GOTO 160
          IF (ITERC .EQ. ITERMAX) GOTO 160
          TEMP=GG/GGSAV
          DD=ZERO
          DS=ZERO
          SS=ZERO
          DO 80 I=1,N
          D(I)=TEMP*D(I)-G(I)-HS(I)
          DD=DD+D(I)**2
          DS=DS+D(I)*STEP(I)
   80     SS=SS+STEP(I)**2
          IF (DS .LE. ZERO) GOTO 160
          IF (SS .LT. DELSQ) GOTO 40
      END IF
      CRVMIN=ZERO
      ITERSW=ITERC
C
C     Test whether an alternative iteration is required.
C
   90 IF (GG .LE. 1.0D-4*GGBEG) GOTO 160
      SG=ZERO
      SHS=ZERO
      DO 100 I=1,N
      SG=SG+STEP(I)*G(I)
  100 SHS=SHS+STEP(I)*HS(I)
      SGK=SG+SHS
      ANGTEST=SGK/DSQRT(GG*DELSQ)
      IF (ANGTEST .LE. -0.99D0) GOTO 160
C
C     Begin the alternative iteration by calculating D and HD and some
C     scalar products.
C
      ITERC=ITERC+1
      TEMP=DSQRT(DELSQ*GG-SGK*SGK)
      TEMPA=DELSQ/TEMP
      TEMPB=SGK/TEMP
      DO 110 I=1,N
  110 D(I)=TEMPA*(G(I)+HS(I))-TEMPB*STEP(I)
      GOTO 170
  120 DG=ZERO
      DHD=ZERO
      DHS=ZERO
      DO 130 I=1,N
      DG=DG+D(I)*G(I)
      DHD=DHD+HD(I)*D(I)
  130 DHS=DHS+HD(I)*STEP(I)
C
C     Seek the value of the angle that minimizes Q.
C
      CF=HALF*(SHS-DHD)
      QBEG=SG+CF
      QSAV=QBEG
      QMIN=QBEG
      ISAVE=0
      IU=49
      TEMP=TWOPI/DFLOAT(IU+1)
      DO 140 I=1,IU
      ANGLE=DFLOAT(I)*TEMP
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      QNEW=(SG+CF*CTH)*CTH+(DG+DHS*CTH)*STH
      IF (QNEW .LT. QMIN) THEN
          QMIN=QNEW
          ISAVE=I
          TEMPA=QSAV
      ELSE IF (I .EQ. ISAVE+1) THEN
          TEMPB=QNEW
      END IF
  140 QSAV=QNEW
      IF (ISAVE .EQ. ZERO) TEMPA=QNEW
      IF (ISAVE .EQ. IU) TEMPB=QBEG
      ANGLE=ZERO
      IF (TEMPA .NE. TEMPB) THEN
          TEMPA=TEMPA-QMIN
          TEMPB=TEMPB-QMIN
          ANGLE=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
      END IF
      ANGLE=TEMP*(DFLOAT(ISAVE)+ANGLE)
C
C     Calculate the new STEP and HS. Then test for convergence.
C
      CTH=DCOS(ANGLE)
      STH=DSIN(ANGLE)
      REDUC=QBEG-(SG+CF*CTH)*CTH-(DG+DHS*CTH)*STH
      GG=ZERO
      DO 150 I=1,N
      STEP(I)=CTH*STEP(I)+STH*D(I)
      HS(I)=CTH*HS(I)+STH*HD(I)
  150 GG=GG+(G(I)+HS(I))**2
      QRED=QRED+REDUC
      RATIO=REDUC/QRED
      IF (ITERC .LT. ITERMAX .AND. RATIO .GT. 0.01D0) GOTO 90
  160 RETURN
C
C     The following instructions act as a subroutine for setting the vector
C     HD to the vector D multiplied by the second derivative matrix of Q.
C     They are called from three different places, which are distinguished
C     by the value of ITERC.
C
  170 DO 180 I=1,N
  180 HD(I)=ZERO
      DO 200 K=1,NPT
      TEMP=ZERO
      DO 190 J=1,N
  190 TEMP=TEMP+XPT(K,J)*D(J)
      TEMP=TEMP*PQ(K)
      DO 200 I=1,N
  200 HD(I)=HD(I)+TEMP*XPT(K,I)
      IH=0
      DO 210 J=1,N
      DO 210 I=1,J
      IH=IH+1
      IF (I .LT. J) HD(J)=HD(J)+HQ(IH)*D(I)
  210 HD(I)=HD(I)+HQ(IH)*D(J)
      IF (ITERC .EQ. 0) GOTO 20
      IF (ITERC .LE. ITERSW) GOTO 50
      GOTO 120
      END
C
C
C
C
C
      SUBROUTINE NEWUOUPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
C
C     The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
C     interpolation point that has index KNEW. On entry, VLAG contains the
C     components of the vector Theta*Wcheck+e_b of the updating formula
C     (6.11), and BETA holds the value of the parameter that has this name.
C     The vector W is used for working space.
C
C     Set some constants.
C
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
C
C     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
C
      JL=1
      DO 20 J=2,NPTM
      IF (J .EQ. IDZ) THEN
          JL=IDZ
      ELSE IF (ZMAT(KNEW,J) .NE. ZERO) THEN
          TEMP=DSQRT(ZMAT(KNEW,JL)**2+ZMAT(KNEW,J)**2)
          TEMPA=ZMAT(KNEW,JL)/TEMP
          TEMPB=ZMAT(KNEW,J)/TEMP
          DO 10 I=1,NPT
          TEMP=TEMPA*ZMAT(I,JL)+TEMPB*ZMAT(I,J)
          ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,JL)
   10     ZMAT(I,JL)=TEMP
          ZMAT(KNEW,J)=ZERO
      END IF
   20 CONTINUE
C
C     Put the first NPT components of the KNEW-th column of HLAG into W,
C     and calculate the parameters of the updating formula.
C
      TEMPA=ZMAT(KNEW,1)
      IF (IDZ .GE. 2) TEMPA=-TEMPA
      IF (JL .GT. 1) TEMPB=ZMAT(KNEW,JL)
      DO 30 I=1,NPT
      W(I)=TEMPA*ZMAT(I,1)
      IF (JL .GT. 1) W(I)=W(I)+TEMPB*ZMAT(I,JL)
   30 CONTINUE
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      TAUSQ=TAU*TAU
      DENOM=ALPHA*BETA+TAUSQ
      VLAG(KNEW)=VLAG(KNEW)-ONE
C
C     Complete the updating of ZMAT when there is only one nonzero element
C     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
C     then the first column of ZMAT will be exchanged with another one later.
C
      IFLAG=0
      IF (JL .EQ. 1) THEN
          TEMP=DSQRT(DABS(DENOM))
          TEMPB=TEMPA/TEMP
          TEMPA=TAU/TEMP
          DO 40 I=1,NPT
   40     ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
          IF (IDZ .EQ. 1 .AND. TEMP .LT. ZERO) IDZ=2
          IF (IDZ .GE. 2 .AND. TEMP .GE. ZERO) IFLAG=1
      ELSE
C
C     Complete the updating of ZMAT in the alternative case.
C
          JA=1
          IF (BETA .GE. ZERO) JA=JL
          JB=JL+1-JA
          TEMP=ZMAT(KNEW,JB)/DENOM
          TEMPA=TEMP*BETA
          TEMPB=TEMP*TAU
          TEMP=ZMAT(KNEW,JA)
          SCALA=ONE/DSQRT(DABS(BETA)*TEMP*TEMP+TAUSQ)
          SCALB=SCALA*DSQRT(DABS(DENOM))
          DO 50 I=1,NPT
          ZMAT(I,JA)=SCALA*(TAU*ZMAT(I,JA)-TEMP*VLAG(I))
   50     ZMAT(I,JB)=SCALB*(ZMAT(I,JB)-TEMPA*W(I)-TEMPB*VLAG(I))
          IF (DENOM .LE. ZERO) THEN
              IF (BETA .LT. ZERO) IDZ=IDZ+1
              IF (BETA .GE. ZERO) IFLAG=1
          END IF
      END IF
C
C     IDZ is reduced in the following case, and usually the first column
C     of ZMAT is exchanged with a later one.
C
      IF (IFLAG .EQ. 1) THEN
          IDZ=IDZ-1
          DO 60 I=1,NPT
          TEMP=ZMAT(I,1)
          ZMAT(I,1)=ZMAT(I,IDZ)
   60     ZMAT(I,IDZ)=TEMP
      END IF
C
C     Finally, update the matrix BMAT.
C
      DO 70 J=1,N
      JP=NPT+J
      W(JP)=BMAT(KNEW,J)
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
      DO 70 I=1,JP
      BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
      IF (I .GT. NPT) BMAT(JP,I-NPT)=BMAT(I,J)
   70 CONTINUE
      RETURN
      END
      SUBROUTINE BOBYQA (N,NPT,X,XL,XU,RHOBEG,RHOEND,LSTDOUT,
     1  MAXFUN,W,CALFUN,OUTFILE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),W(*)
C     Definition of stream control variables 256 appears in gsl_findroot.c
      EXTERNAL CALFUN
      INTEGER       LSTDOUT
      CHARACTER*256 OUTFILE
      IF (LSTDOUT .gt. 1) THEN
          IPRINT = 12
          OPEN (IPRINT, FILE=OUTFILE, STATUS='OLD')
      ELSE
          IPRINT = 0
      ENDIF
C     End of it. KMK

C
C     This subroutine seeks the least value of a function of many variables,
C     by applying a trust region method that forms quadratic models by
C     interpolation. There is usually some freedom in the interpolation
C     conditions, which is taken up by minimizing the Frobenius norm of
C     the change to the second derivative of the model, beginning with the
C     zero matrix. The values of the variables are constrained by upper and
C     lower bounds. The arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT is the number of interpolation conditions. Its value must be in
C       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
C       recommended.
C     Initial values of the variables must be set in X(1),X(2),...,X(N). They
C       will be changed to the values that give the least calculated F.
C     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
C       bounds, respectively, on X(I). The construction of quadratic models
C       requires XL(I) to be strictly less than XU(I) for each I. Further,
C       the contribution to a model from changes to the I-th variable is
C       damaged severely by rounding errors if XU(I)-XL(I) is too small.
C     RHOBEG and RHOEND must be set to the initial and final values of a trust
C       region radius, so both must be positive with RHOEND no greater than
C       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
C       expected change to a variable, while RHOEND should indicate the
C       accuracy that is required in the final values of the variables. An
C       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
C       is less than 2*RHOBEG.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, each new
C       value of RHO is printed, with the best vector of variables so far and
C       the corresponding value of the objective function. Further, each new
C       value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
C     The array W will be used for working space. Its length must be at least
C       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
C
C     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
C     F to the value of the objective function for the current values of the
C     variables X(1),X(2),...,X(N), which are generated automatically in a
C     way that satisfies the bounds given in XL and XU.
C
C     Return if the value of NPT is unacceptable.
C
      NP=N+1
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
          if (IPRINT .gt. 0) WRITE(IPRINT,10)
   10     FORMAT (/4X,'Return from BOBYQA because NPT is not in',
     1      ' the required interval')
          GO TO 40
      END IF
C
C     Partition the working space array, so that different parts of it can
C     be treated separately during the calculation of BOBYQB. The partition
C     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
C     space that is taken by the last array in the argument list of BOBYQB.
C
      NDIM=NPT+N
      IXB=1
      IXP=IXB+N
      IFV=IXP+N*NPT
      IXO=IFV+NPT
      IGO=IXO+N
      IHQ=IGO+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ISL=IZMAT+NPT*(NPT-NP)
      ISU=ISL+N
      IXN=ISU+N
      IXA=IXN+N
      ID=IXA+N
      IVL=ID+N
      IW=IVL+NDIM
C
C     Return if there is insufficient space between the bounds. Modify the
C     initial X if necessary in order to avoid conflicts between the bounds
C     and the construction of the first quadratic model. The lower and upper
C     bounds on moves from the updated X are set now, in the ISL and ISU
C     partitions of W, in order to provide useful and exact information about
C     components of X that become within distance RHOBEG from their bounds.
C
      ZERO=0.0D0
      DO 30 J=1,N
      TEMP=XU(J)-XL(J)
      IF (TEMP .LT. RHOBEG+RHOBEG) THEN
          if (IPRINT .gt. 0) WRITE(IPRINT,20)
   20     FORMAT (/4X,'Return from BOBYQA because one of the',
     1      ' differences XU(I)-XL(I)'/6X,' is less than 2*RHOBEG.')
          GO TO 40
      END IF
      JSL=ISL+J-1
      JSU=JSL+N
      W(JSL)=XL(J)-X(J)
      W(JSU)=XU(J)-X(J)
      IF (W(JSL) .GE. -RHOBEG) THEN
          IF (W(JSL) .GE. ZERO) THEN
              X(J)=XL(J)
              W(JSL)=ZERO
              W(JSU)=TEMP
          ELSE
              X(J)=XL(J)+RHOBEG
              W(JSL)=-RHOBEG
              W(JSU)=DMAX1(XU(J)-X(J),RHOBEG)
          END IF
      ELSE IF (W(JSU) .LE. RHOBEG) THEN
          IF (W(JSU) .LE. ZERO) THEN
              X(J)=XU(J)
              W(JSL)=-TEMP
              W(JSU)=ZERO
          ELSE
              X(J)=XU(J)-RHOBEG
              W(JSL)=DMIN1(XL(J)-X(J),-RHOBEG)
              W(JSU)=RHOBEG
          END IF
      END IF
   30 CONTINUE
C
C     Make the call of BOBYQB.
C
      CALL BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB),
     1  W(IXP),W(IFV),W(IXO),W(IGO),W(IHQ),W(IPQ),W(IBMAT),W(IZMAT),
     2  NDIM,W(ISL),W(ISU),W(IXN),W(IXA),W(ID),W(IVL),W(IW),
     3  CALFUN)
   40 CONTINUE
      IF (LSTDOUT .gt. 1)
     * CLOSE (IPRINT)
      RETURN
      END
      SUBROUTINE ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
     1  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,GLAG,HCOL,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPT(NPT,*),XOPT(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),
     1  SU(*),XNEW(*),XALT(*),GLAG(*),HCOL(*),W(*)
C
C     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
C       the same meanings as the corresponding arguments of BOBYQB.
C     KOPT is the index of the optimal interpolation point.
C     KNEW is the index of the interpolation point that is going to be moved.
C     ADELT is the current trust region bound.
C     XNEW will be set to a suitable new position for the interpolation point
C       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
C       bounds and it should provide a large denominator in the next call of
C       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
C       straight lines through XOPT and another interpolation point.
C     XALT also provides a large value of the modulus of the KNEW-th Lagrange
C       function subject to the constraints that have been mentioned, its main
C       difference from XNEW being that XALT-XOPT is a constrained version of
C       the Cauchy step within the trust region. An exception is that XALT is
C       not calculated if all components of GLAG (see below) are zero.
C     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
C     CAUCHY will be set to the square of the KNEW-th Lagrange function at
C       the step XALT-XOPT from XOPT for the vector XALT that is returned,
C       except that CAUCHY is set to zero if XALT is not calculated.
C     GLAG is a working space vector of length N for the gradient of the
C       KNEW-th Lagrange function at XOPT.
C     HCOL is a working space vector of length NPT for the second derivative
C       coefficients of the KNEW-th Lagrange function.
C     W is a working space vector of length 2N that is going to hold the
C       constrained Cauchy step from XOPT of the Lagrange function, followed
C       by the downhill version of XALT when the uphill step is calculated.
C
C     Set the first NPT components of W to the leading elements of the
C     KNEW-th column of the H matrix.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      CONST=ONE+DSQRT(2.0D0)
      DO 10 K=1,NPT
   10 HCOL(K)=ZERO
      DO 20 J=1,NPT-N-1
      TEMP=ZMAT(KNEW,J)
      DO 20 K=1,NPT
   20 HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
      ALPHA=HCOL(KNEW)
      HA=HALF*ALPHA
C
C     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
C
      DO 30 I=1,N
   30 GLAG(I)=BMAT(KNEW,I)
      DO 50 K=1,NPT
      TEMP=ZERO
      DO 40 J=1,N
   40 TEMP=TEMP+XPT(K,J)*XOPT(J)
      TEMP=HCOL(K)*TEMP
      DO 50 I=1,N
   50 GLAG(I)=GLAG(I)+TEMP*XPT(K,I)
C
C     Search for a large denominator along the straight lines through XOPT
C     and another interpolation point. SLBD and SUBD will be lower and upper
C     bounds on the step along each of these lines in turn. PREDSQ will be
C     set to the square of the predicted denominator for each line. PRESAV
C     will be set to the largest admissible value of PREDSQ that occurs.
C
      PRESAV=ZERO
      DO 80 K=1,NPT
      IF (K .EQ. KOPT) GOTO 80
      DDERIV=ZERO
      DISTSQ=ZERO
      DO 60 I=1,N
      TEMP=XPT(K,I)-XOPT(I)
      DDERIV=DDERIV+GLAG(I)*TEMP
   60 DISTSQ=DISTSQ+TEMP*TEMP
      SUBD=ADELT/DSQRT(DISTSQ)
      SLBD=-SUBD
      ILBD=0
      IUBD=0
      SUMIN=DMIN1(ONE,SUBD)
C
C     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
C
      DO 70 I=1,N
      TEMP=XPT(K,I)-XOPT(I)
      IF (TEMP .GT. ZERO) THEN
          IF (SLBD*TEMP .LT. SL(I)-XOPT(I)) THEN
              SLBD=(SL(I)-XOPT(I))/TEMP
              ILBD=-I
          END IF
          IF (SUBD*TEMP .GT. SU(I)-XOPT(I)) THEN
              SUBD=DMAX1(SUMIN,(SU(I)-XOPT(I))/TEMP)
              IUBD=I
          END IF
      ELSE IF (TEMP .LT. ZERO) THEN
          IF (SLBD*TEMP .GT. SU(I)-XOPT(I)) THEN
              SLBD=(SU(I)-XOPT(I))/TEMP
              ILBD=I
          END IF
          IF (SUBD*TEMP .LT. SL(I)-XOPT(I)) THEN
              SUBD=DMAX1(SUMIN,(SL(I)-XOPT(I))/TEMP)
              IUBD=-I
          END IF
      END IF
   70 CONTINUE
C
C     Seek a large modulus of the KNEW-th Lagrange function when the index
C     of the other interpolation point on the line through XOPT is KNEW.
C
      IF (K .EQ. KNEW) THEN
          DIFF=DDERIV-ONE
          STEP=SLBD
          VLAG=SLBD*(DDERIV-SLBD*DIFF)
          ISBD=ILBD
          TEMP=SUBD*(DDERIV-SUBD*DIFF)
          IF (DABS(TEMP) .GT. DABS(VLAG)) THEN
              STEP=SUBD
              VLAG=TEMP
              ISBD=IUBD
          END IF
          TEMPD=HALF*DDERIV
          TEMPA=TEMPD-DIFF*SLBD
          TEMPB=TEMPD-DIFF*SUBD
          IF (TEMPA*TEMPB .LT. ZERO) THEN
              TEMP=TEMPD*TEMPD/DIFF
              IF (DABS(TEMP) .GT. DABS(VLAG)) THEN
                  STEP=TEMPD/DIFF
                  VLAG=TEMP
                  ISBD=0
              END IF
          END IF
C
C     Search along each of the other lines through XOPT and another point.
C
      ELSE
          STEP=SLBD
          VLAG=SLBD*(ONE-SLBD)
          ISBD=ILBD
          TEMP=SUBD*(ONE-SUBD)
          IF (DABS(TEMP) .GT. DABS(VLAG)) THEN
              STEP=SUBD
              VLAG=TEMP
              ISBD=IUBD
          END IF
          IF (SUBD .GT. HALF) THEN
              IF (DABS(VLAG) .LT. 0.25D0) THEN
                  STEP=HALF
                  VLAG=0.25D0
                  ISBD=0
              END IF
          END IF
          VLAG=VLAG*DDERIV
      END IF
C
C     Calculate PREDSQ for the current line search and maintain PRESAV.
C
      TEMP=STEP*(ONE-STEP)*DISTSQ
      PREDSQ=VLAG*VLAG*(VLAG*VLAG+HA*TEMP*TEMP)
      IF (PREDSQ .GT. PRESAV) THEN
          PRESAV=PREDSQ
          KSAV=K
          STPSAV=STEP
          IBDSAV=ISBD
      END IF
   80 CONTINUE
C
C     Construct XNEW in a way that satisfies the bound constraints exactly.
C
      DO 90 I=1,N
      TEMP=XOPT(I)+STPSAV*(XPT(KSAV,I)-XOPT(I))
   90 XNEW(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
      IF (IBDSAV .LT. 0) XNEW(-IBDSAV)=SL(-IBDSAV)
      IF (IBDSAV .GT. 0) XNEW(IBDSAV)=SU(IBDSAV)
C
C     Prepare for the iterative method that assembles the constrained Cauchy
C     step in W. The sum of squares of the fixed components of W is formed in
C     WFIXSQ, and the free components of W are set to BIGSTP.
C
      BIGSTP=ADELT+ADELT
      IFLAG=0
  100 WFIXSQ=ZERO
      GGFREE=ZERO
      DO 110 I=1,N
      W(I)=ZERO
      TEMPA=DMIN1(XOPT(I)-SL(I),GLAG(I))
      TEMPB=DMAX1(XOPT(I)-SU(I),GLAG(I))
      IF (TEMPA .GT. ZERO .OR. TEMPB .LT. ZERO) THEN
          W(I)=BIGSTP
          GGFREE=GGFREE+GLAG(I)**2
      END IF
  110 CONTINUE
      IF (GGFREE .EQ. ZERO) THEN
          CAUCHY=ZERO
          GOTO 200
      END IF
C
C     Investigate whether more components of W can be fixed.
C
  120 TEMP=ADELT*ADELT-WFIXSQ
      IF (TEMP .GT. ZERO) THEN
          WSQSAV=WFIXSQ
          STEP=DSQRT(TEMP/GGFREE)
          GGFREE=ZERO
          DO 130 I=1,N
          IF (W(I) .EQ. BIGSTP) THEN
              TEMP=XOPT(I)-STEP*GLAG(I)
              IF (TEMP .LE. SL(I)) THEN
                  W(I)=SL(I)-XOPT(I)
                  WFIXSQ=WFIXSQ+W(I)**2
              ELSE IF (TEMP .GE. SU(I)) THEN
                  W(I)=SU(I)-XOPT(I)
                  WFIXSQ=WFIXSQ+W(I)**2
              ELSE
                  GGFREE=GGFREE+GLAG(I)**2
              END IF
          END IF
  130     CONTINUE
          IF (WFIXSQ .GT. WSQSAV .AND. GGFREE .GT. ZERO) GOTO 120
      END IF
C
C     Set the remaining free components of W and all components of XALT,
C     except that W may be scaled later.
C
      GW=ZERO
      DO 140 I=1,N
      IF (W(I) .EQ. BIGSTP) THEN
          W(I)=-STEP*GLAG(I)
          XALT(I)=DMAX1(SL(I),DMIN1(SU(I),XOPT(I)+W(I)))
      ELSE IF (W(I) .EQ. ZERO) THEN
          XALT(I)=XOPT(I)
      ELSE IF (GLAG(I) .GT. ZERO) THEN
          XALT(I)=SL(I)
      ELSE
          XALT(I)=SU(I)
      END IF
  140 GW=GW+GLAG(I)*W(I)
C
C     Set CURV to the curvature of the KNEW-th Lagrange function along W.
C     Scale W by a factor less than one if that can reduce the modulus of
C     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
C     the square of this function.
C
      CURV=ZERO
      DO 160 K=1,NPT
      TEMP=ZERO
      DO 150 J=1,N
  150 TEMP=TEMP+XPT(K,J)*W(J)
  160 CURV=CURV+HCOL(K)*TEMP*TEMP
      IF (IFLAG .EQ. 1) CURV=-CURV
      IF (CURV .GT. -GW .AND. CURV .LT. -CONST*GW) THEN
          SCALE=-GW/CURV
          DO 170 I=1,N
          TEMP=XOPT(I)+SCALE*W(I)
  170     XALT(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
          CAUCHY=(HALF*GW*SCALE)**2
      ELSE
          CAUCHY=(GW+HALF*CURV)**2
      END IF
C
C     If IFLAG is zero, then XALT is calculated as before after reversing
C     the sign of GLAG. Thus two XALT vectors become available. The one that
C     is chosen is the one that gives the larger value of CAUCHY.
C
      IF (IFLAG .EQ. 0) THEN
          DO 180 I=1,N
          GLAG(I)=-GLAG(I)
  180     W(N+I)=XALT(I)
          CSAVE=CAUCHY
          IFLAG=1
          GOTO 100
      END IF
      IF (CSAVE .GT. CAUCHY) THEN
          DO 190 I=1,N
  190     XALT(I)=W(N+I)
          CAUCHY=CSAVE
      END IF
  200 RETURN
      END
      SUBROUTINE BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,
     1  MAXFUN,XBASE,XPT,FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,
     2  SL,SU,XNEW,XALT,D,VLAG,W,CALFUN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),
     1  XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),
     2  SL(*),SU(*),XNEW(*),XALT(*),D(*),VLAG(*),W(*)
      EXTERNAL CALFUN
C
C     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
C       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
C     XBASE holds a shift of origin that should reduce the contributions
C       from rounding errors to values of the model and Lagrange functions.
C     XPT is a two-dimensional array that holds the coordinates of the
C       interpolation points relative to XBASE.
C     FVAL holds the values of F at the interpolation points.
C     XOPT is set to the displacement from XBASE of the trust region centre.
C     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
C     HQ holds the explicit second derivatives of the quadratic model.
C     PQ contains the parameters of the implicit second derivatives of the
C       quadratic model.
C     BMAT holds the last N columns of H.
C     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
C       this factorization being ZMAT times ZMAT^T, which provides both the
C       correct rank and positive semi-definiteness.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
C       All the components of every XOPT are going to satisfy the bounds
C       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
C       XOPT is on a constraint boundary.
C     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
C       vector of variables for the next call of CALFUN. XNEW also satisfies
C       the SL and SU constraints in the way that has just been mentioned.
C     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
C       in order to increase the denominator in the updating of UPDATE.
C     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
C     VLAG contains the values of the Lagrange functions at a new point X.
C       They are part of a product that requires VLAG to be of length NDIM.
C     W is a one-dimensional array that is used for working space. Its length
C       must be at least 3*NDIM = 3*(NPT+N).
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TEN=10.0D0
      TENTH=0.1D0
      TWO=2.0D0
      ZERO=0.0D0
      NP=N+1
      NPTM=NPT-NP
      NH=(N*NP)/2
C
C     The call of BOPRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
C     BMAT and ZMAT for the first iteration, with the corresponding values of
C     of NF and KOPT, which are the number of calls of CALFUN so far and the
C     index of the interpolation point at the trust region centre. Then the
C     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
C     less than NPT. GOPT will be updated if KOPT is different from KBASE.
C
      CALL BOPRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,XPT,
     1  FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT,CALFUN)
      XOPTSQ=ZERO
      DO 10 I=1,N
      XOPT(I)=XPT(KOPT,I)
   10 XOPTSQ=XOPTSQ+XOPT(I)**2
      FSAVE=FVAL(1)
      IF (NF .LT. NPT) THEN
          IF (IPRINT .GT. 0) WRITE(IPRINT,390)
          GOTO 720
      END IF
      KBASE=1
C
C     Complete the settings that are required for the iterative procedure.
C
      RHO=RHOBEG
      DELTA=RHO
      NRESC=NF
      NTRITS=0
      DIFFA=ZERO
      DIFFB=ZERO
      ITEST=0
      NFSAV=NF
C
C     Update GOPT if necessary before the first iteration and after each
C     call of RESCUE that makes a call of CALFUN.
C
   20 IF (KOPT .NE. KBASE) THEN
          IH=0
          DO 30 J=1,N
          DO 30 I=1,J
          IH=IH+1
          IF (I .LT. J) GOPT(J)=GOPT(J)+HQ(IH)*XOPT(I)
   30     GOPT(I)=GOPT(I)+HQ(IH)*XOPT(J)
          IF (NF .GT. NPT) THEN
              DO 50 K=1,NPT
              TEMP=ZERO
              DO 40 J=1,N
   40         TEMP=TEMP+XPT(K,J)*XOPT(J)
              TEMP=PQ(K)*TEMP
              DO 50 I=1,N
   50         GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
          END IF
      END IF
C
C     Generate the next point in the trust region that provides a small value
C     of the quadratic model subject to the constraints on the variables.
C     The integer NTRITS is set to the number "trust region" iterations that
C     have occurred since the last "alternative" iteration. If the length
C     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
C     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
C
   60 CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D,
     1  W,W(NP),W(NP+N),W(NP+2*N),W(NP+3*N),DSQ,CRVMIN)
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      IF (DNORM .LT. HALF*RHO) THEN
          NTRITS=-1
          DISTSQ=(TEN*RHO)**2
          IF (NF .LE. NFSAV+2) GOTO 650
C
C     The following choice between labels 650 and 680 depends on whether or
C     not our work with the current RHO seems to be complete. Either RHO is
C     decreased or termination occurs if the errors in the quadratic model at
C     the last three interpolation points compare favourably with predictions
C     of likely improvements to the model within distance HALF*RHO of XOPT.
C
          ERRBIG=DMAX1(DIFFA,DIFFB,DIFFC)
          FRHOSQ=0.125D0*RHO*RHO
          IF (CRVMIN .GT. ZERO .AND. ERRBIG .GT. FRHOSQ*CRVMIN)
     1       GOTO 650
          BDTOL=ERRBIG/RHO
          DO 80 J=1,N
          BDTEST=BDTOL
          IF (XNEW(J) .EQ. SL(J)) BDTEST=W(J)
          IF (XNEW(J) .EQ. SU(J)) BDTEST=-W(J)
          IF (BDTEST .LT. BDTOL) THEN
              CURV=HQ((J+J*J)/2)
              DO 70 K=1,NPT
   70         CURV=CURV+PQ(K)*XPT(K,J)**2
              BDTEST=BDTEST+HALF*CURV*RHO
              IF (BDTEST .LT. BDTOL) GOTO 650
          END IF
   80     CONTINUE
          GOTO 680
      END IF
      NTRITS=NTRITS+1
C
C     Severe cancellation is likely to occur if XOPT is too far from XBASE.
C     If the following test holds, then XBASE is shifted so that XOPT becomes
C     zero. The appropriate changes are made to BMAT and to the second
C     derivatives of the current model, beginning with the changes to BMAT
C     that do not depend on ZMAT. VLAG is used temporarily for working space.
C
   90 IF (DSQ .LE. 1.0D-3*XOPTSQ) THEN
          FRACSQ=0.25D0*XOPTSQ
          SUMPQ=ZERO
          DO 110 K=1,NPT
          SUMPQ=SUMPQ+PQ(K)
          SUM=-HALF*XOPTSQ
          DO 100 I=1,N
  100     SUM=SUM+XPT(K,I)*XOPT(I)
          W(NPT+K)=SUM
          TEMP=FRACSQ-HALF*SUM
          DO 110 I=1,N
          W(I)=BMAT(K,I)
          VLAG(I)=SUM*XPT(K,I)+TEMP*XOPT(I)
          IP=NPT+I
          DO 110 J=1,I
  110     BMAT(IP,J)=BMAT(IP,J)+W(I)*VLAG(J)+VLAG(I)*W(J)
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO 150 JJ=1,NPTM
          SUMZ=ZERO
          SUMW=ZERO
          DO 120 K=1,NPT
          SUMZ=SUMZ+ZMAT(K,JJ)
          VLAG(K)=W(NPT+K)*ZMAT(K,JJ)
  120     SUMW=SUMW+VLAG(K)
          DO 140 J=1,N
          SUM=(FRACSQ*SUMZ-HALF*SUMW)*XOPT(J)
          DO 130 K=1,NPT
  130     SUM=SUM+VLAG(K)*XPT(K,J)
          W(J)=SUM
          DO 140 K=1,NPT
  140     BMAT(K,J)=BMAT(K,J)+SUM*ZMAT(K,JJ)
          DO 150 I=1,N
          IP=I+NPT
          TEMP=W(I)
          DO 150 J=1,I
  150     BMAT(IP,J)=BMAT(IP,J)+TEMP*W(J)
C
C     The following instructions complete the shift, including the changes
C     to the second derivative parameters of the quadratic model.
C
          IH=0
          DO 170 J=1,N
          W(J)=-HALF*SUMPQ*XOPT(J)
          DO 160 K=1,NPT
          W(J)=W(J)+PQ(K)*XPT(K,J)
  160     XPT(K,J)=XPT(K,J)-XOPT(J)
          DO 170 I=1,J
          IH=IH+1
          HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
  170     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 180 I=1,N
          XBASE(I)=XBASE(I)+XOPT(I)
          XNEW(I)=XNEW(I)-XOPT(I)
          SL(I)=SL(I)-XOPT(I)
          SU(I)=SU(I)-XOPT(I)
  180     XOPT(I)=ZERO
          XOPTSQ=ZERO
      END IF
      IF (NTRITS .EQ. 0) GOTO 210
      GOTO 230
C
C     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
C     more expensive than the previous shift, because new matrices BMAT and
C     ZMAT are generated from scratch, which may include the replacement of
C     interpolation points whose positions seem to be causing near linear
C     dependence in the interpolation conditions. Therefore RESCUE is called
C     only if rounding errors have reduced by at least a factor of two the
C     denominator of the formula for updating the H matrix. It provides a
C     useful safeguard, but is not invoked in most applications of BOBYQA.
C
  190 NFSAV=NF
      KBASE=KOPT
      CALL RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,FVAL,
     1  XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,KOPT,
     2  VLAG,W,W(N+NP),W(NDIM+NP),CALFUN)
C
C     XOPT is updated now in case the branch below to label 720 is taken.
C     Any updating of GOPT occurs after the branch below to label 20, which
C     leads to a trust region iteration as does the branch to label 60.
C
      XOPTSQ=ZERO
      IF (KOPT .NE. KBASE) THEN
          DO 200 I=1,N
          XOPT(I)=XPT(KOPT,I)
  200     XOPTSQ=XOPTSQ+XOPT(I)**2
      END IF
      IF (NF .LT. 0) THEN
          NF=MAXFUN
          IF (IPRINT .GT. 0) WRITE(IPRINT,390)
          GOTO 720
      END IF
      NRESC=NF
      IF (NFSAV .LT. NF) THEN
          NFSAV=NF
          GOTO 20
      END IF
      IF (NTRITS .GT. 0) GOTO 60
C
C     Pick two alternative vectors of variables, relative to XBASE, that
C     are suitable as new positions of the KNEW-th interpolation point.
C     Firstly, XNEW is set to the point on a line through XOPT and another
C     interpolation point that minimizes the predicted value of the next
C     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
C     and SU bounds. Secondly, XALT is set to the best feasible point on
C     a constrained version of the Cauchy step of the KNEW-th Lagrange
C     function, the corresponding value of the square of this function
C     being returned in CAUCHY. The choice between these alternatives is
C     going to be made when the denominator is calculated.
C
  210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT,
     1  KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
      DO 220 I=1,N
  220 D(I)=XNEW(I)-XOPT(I)
C
C     Calculate VLAG and BETA for the current choice of D. The scalar
C     product of D with XPT(K,.) is going to be held in W(NPT+K) for
C     use when VQUAD is calculated.
C
  230 DO 250 K=1,NPT
      SUMA=ZERO
      SUMB=ZERO
      SUM=ZERO
      DO 240 J=1,N
      SUMA=SUMA+XPT(K,J)*D(J)
      SUMB=SUMB+XPT(K,J)*XOPT(J)
  240 SUM=SUM+BMAT(K,J)*D(J)
      W(K)=SUMA*(HALF*SUMA+SUMB)
      VLAG(K)=SUM
  250 W(NPT+K)=SUMA
      BETA=ZERO
      DO 270 JJ=1,NPTM
      SUM=ZERO
      DO 260 K=1,NPT
  260 SUM=SUM+ZMAT(K,JJ)*W(K)
      BETA=BETA-SUM*SUM
      DO 270 K=1,NPT
  270 VLAG(K)=VLAG(K)+SUM*ZMAT(K,JJ)
      DSQ=ZERO
      BSUM=ZERO
      DX=ZERO
      DO 300 J=1,N
      DSQ=DSQ+D(J)**2
      SUM=ZERO
      DO 280 K=1,NPT
  280 SUM=SUM+W(K)*BMAT(K,J)
      BSUM=BSUM+SUM*D(J)
      JP=NPT+J
      DO 290 I=1,N
  290 SUM=SUM+BMAT(JP,I)*D(I)
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D(J)
  300 DX=DX+D(J)*XOPT(J)
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     If NTRITS is zero, the denominator may be increased by replacing
C     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
C     rounding errors have damaged the chosen denominator.
C
      IF (NTRITS .EQ. 0) THEN
          DENOM=VLAG(KNEW)**2+ALPHA*BETA
          IF (DENOM .LT. CAUCHY .AND. CAUCHY .GT. ZERO) THEN
              DO 310 I=1,N
              XNEW(I)=XALT(I)
  310         D(I)=XNEW(I)-XOPT(I)
              CAUCHY=ZERO
              GO TO 230
          END IF
          IF (DENOM .LE. HALF*VLAG(KNEW)**2) THEN
              IF (NF .GT. NRESC) GOTO 190
              IF (IPRINT .GT. 0) WRITE(IPRINT,320)
  320         FORMAT (/5X,'Return from BOBYQA because of much',
     1          ' cancellation in a denominator.')
              GOTO 720
          END IF
C
C     Alternatively, if NTRITS is positive, then set KNEW to the index of
C     the next interpolation point to be deleted to make room for a trust
C     region step. Again RESCUE may be called if rounding errors have damaged
C     the chosen denominator, which is the reason for attempting to select
C     KNEW before calculating the next value of the objective function.
C
      ELSE
          DELSQ=DELTA*DELTA
          SCADEN=ZERO
          BIGLSQ=ZERO
          KNEW=0
          DO 350 K=1,NPT
          IF (K .EQ. KOPT) GOTO 350
          HDIAG=ZERO
          DO 330 JJ=1,NPTM
  330     HDIAG=HDIAG+ZMAT(K,JJ)**2
          DEN=BETA*HDIAG+VLAG(K)**2
          DISTSQ=ZERO
          DO 340 J=1,N
  340     DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
          TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
          IF (TEMP*DEN .GT. SCADEN) THEN
              SCADEN=TEMP*DEN
              KNEW=K
              DENOM=DEN
          END IF
          BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
  350     CONTINUE
          IF (SCADEN .LE. HALF*BIGLSQ) THEN
              IF (NF .GT. NRESC) GOTO 190
              IF (IPRINT .GT. 0) WRITE(IPRINT,320)
              GOTO 720
          END IF
      END IF
C
C     Put the variables for the next calculation of the objective function
C       in XNEW, with any adjustments for the bounds.
C
C
C     Calculate the value of the objective function at XBASE+XNEW, unless
C       the limit on the number of calculations of F has been reached.
C
  360 DO 380 I=1,N
      X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XNEW(I)),XU(I))
      IF (XNEW(I) .EQ. SL(I)) X(I)=XL(I)
      IF (XNEW(I) .EQ. SU(I)) X(I)=XU(I)
  380 CONTINUE
      IF (NF .GE. MAXFUN) THEN
          IF (IPRINT .GT. 0) WRITE(IPRINT,390)
  390     FORMAT (/4X,'Return from BOBYQA because CALFUN has been',
     1      ' called MAXFUN times.')
          GOTO 720
      END IF
      NF=NF+1
      CALL CALFUN (N,X,F)
      IF (IPRINT .GE. 3) THEN
          WRITE(IPRINT,400) NF,F,(X(I),I=1,N)
  400     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NTRITS .EQ. -1) THEN
          FSAVE=F
          GOTO 720
      END IF
C
C     Use the quadratic model to predict the change in F due to the step D,
C       and set DIFF to the error of this prediction.
C
      FOPT=FVAL(KOPT)
      VQUAD=ZERO
      IH=0
      DO 410 J=1,N
      VQUAD=VQUAD+D(J)*GOPT(J)
      DO 410 I=1,J
      IH=IH+1
      TEMP=D(I)*D(J)
      IF (I .EQ. J) TEMP=HALF*TEMP
  410 VQUAD=VQUAD+HQ(IH)*TEMP
      DO 420 K=1,NPT
  420 VQUAD=VQUAD+HALF*PQ(K)*W(NPT+K)**2
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM .GT. RHO) NFSAV=NF
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (NTRITS .GT. 0) THEN
          IF (VQUAD .GE. ZERO) THEN
              IF (IPRINT .GT. 0) WRITE(IPRINT,430)
  430         FORMAT (/4X,'Return from BOBYQA because a trust',
     1          ' region step has failed to reduce Q.')
              GOTO 720
          END IF
          RATIO=(F-FOPT)/VQUAD
          IF (RATIO .LE. TENTH) THEN
              DELTA=DMIN1(HALF*DELTA,DNORM)
          ELSE IF (RATIO .LE. 0.7D0) THEN
              DELTA=DMAX1(HALF*DELTA,DNORM)
          ELSE
              DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
          END IF
          IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
C
C     Recalculate KNEW and DENOM if the new F is less than FOPT.
C
          IF (F .LT. FOPT) THEN
              KSAV=KNEW
              DENSAV=DENOM
              DELSQ=DELTA*DELTA
              SCADEN=ZERO
              BIGLSQ=ZERO
              KNEW=0
              DO 460 K=1,NPT
              HDIAG=ZERO
              DO 440 JJ=1,NPTM
  440         HDIAG=HDIAG+ZMAT(K,JJ)**2
              DEN=BETA*HDIAG+VLAG(K)**2
              DISTSQ=ZERO
              DO 450 J=1,N
  450         DISTSQ=DISTSQ+(XPT(K,J)-XNEW(J))**2
              TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
              IF (TEMP*DEN .GT. SCADEN) THEN
                  SCADEN=TEMP*DEN
                  KNEW=K
                  DENOM=DEN
              END IF
  460         BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
              IF (SCADEN .LE. HALF*BIGLSQ) THEN
                  KNEW=KSAV
                  DENOM=DENSAV
              END IF
          END IF
      END IF
C
C     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
C     moved. Also update the second derivative terms of the model.
C
      CALL BOBUPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
      IH=0
      PQOLD=PQ(KNEW)
      PQ(KNEW)=ZERO
      DO 470 I=1,N
      TEMP=PQOLD*XPT(KNEW,I)
      DO 470 J=1,I
      IH=IH+1
  470 HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
      DO 480 JJ=1,NPTM
      TEMP=DIFF*ZMAT(KNEW,JJ)
      DO 480 K=1,NPT
  480 PQ(K)=PQ(K)+TEMP*ZMAT(K,JJ)
C
C     Include the new interpolation point, and make the changes to GOPT at
C     the old XOPT that are caused by the updating of the quadratic model.
C
      FVAL(KNEW)=F
      DO 490 I=1,N
      XPT(KNEW,I)=XNEW(I)
  490 W(I)=BMAT(KNEW,I)
      DO 520 K=1,NPT
      SUMA=ZERO
      DO 500 JJ=1,NPTM
  500 SUMA=SUMA+ZMAT(KNEW,JJ)*ZMAT(K,JJ)
      SUMB=ZERO
      DO 510 J=1,N
  510 SUMB=SUMB+XPT(K,J)*XOPT(J)
      TEMP=SUMA*SUMB
      DO 520 I=1,N
  520 W(I)=W(I)+TEMP*XPT(K,I)
      DO 530 I=1,N
  530 GOPT(I)=GOPT(I)+DIFF*W(I)
C
C     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
C
      IF (F .LT. FOPT) THEN
          KOPT=KNEW
          XOPTSQ=ZERO
          IH=0
          DO 540 J=1,N
          XOPT(J)=XNEW(J)
          XOPTSQ=XOPTSQ+XOPT(J)**2
          DO 540 I=1,J
          IH=IH+1
          IF (I .LT. J) GOPT(J)=GOPT(J)+HQ(IH)*D(I)
  540     GOPT(I)=GOPT(I)+HQ(IH)*D(J)
          DO 560 K=1,NPT
          TEMP=ZERO
          DO 550 J=1,N
  550     TEMP=TEMP+XPT(K,J)*D(J)
          TEMP=PQ(K)*TEMP
          DO 560 I=1,N
  560     GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
      END IF
C
C     Calculate the parameters of the least Frobenius norm interpolant to
C     the current data, the gradient of this interpolant at XOPT being put
C     into VLAG(NPT+I), I=1,2,...,N.
C
      IF (NTRITS .GT. 0) THEN
          DO 570 K=1,NPT
          VLAG(K)=FVAL(K)-FVAL(KOPT)
  570     W(K)=ZERO
          DO 590 J=1,NPTM
          SUM=ZERO
          DO 580 K=1,NPT
  580     SUM=SUM+ZMAT(K,J)*VLAG(K)
          DO 590 K=1,NPT
  590     W(K)=W(K)+SUM*ZMAT(K,J)
          DO 610 K=1,NPT
          SUM=ZERO
          DO 600 J=1,N
  600     SUM=SUM+XPT(K,J)*XOPT(J)
          W(K+NPT)=W(K)
  610     W(K)=SUM*W(K)
          GQSQ=ZERO
          GISQ=ZERO
          DO 630 I=1,N
          SUM=ZERO
          DO 620 K=1,NPT
  620     SUM=SUM+BMAT(K,I)*VLAG(K)+XPT(K,I)*W(K)
          IF (XOPT(I) .EQ. SL(I)) THEN
              GQSQ=GQSQ+DMIN1(ZERO,GOPT(I))**2
              GISQ=GISQ+DMIN1(ZERO,SUM)**2
          ELSE IF (XOPT(I) .EQ. SU(I)) THEN
              GQSQ=GQSQ+DMAX1(ZERO,GOPT(I))**2
              GISQ=GISQ+DMAX1(ZERO,SUM)**2
          ELSE
              GQSQ=GQSQ+GOPT(I)**2
              GISQ=GISQ+SUM*SUM
          END IF
  630     VLAG(NPT+I)=SUM
C
C     Test whether to replace the new quadratic model by the least Frobenius
C     norm interpolant, making the replacement if the test is satisfied.
C
          ITEST=ITEST+1
          IF (GQSQ .LT. TEN*GISQ) ITEST=0
          IF (ITEST .GE. 3) THEN
              DO 640 I=1,MAX0(NPT,NH)
              IF (I .LE. N) GOPT(I)=VLAG(NPT+I)
              IF (I .LE. NPT) PQ(I)=W(NPT+I)
              IF (I .LE. NH) HQ(I)=ZERO
              ITEST=0
  640         CONTINUE
          END IF
      END IF
C
C     If a trust region step has provided a sufficient decrease in F, then
C     branch for another trust region calculation. The case NTRITS=0 occurs
C     when the new interpolation point was reached by an alternative step.
C
      IF (NTRITS .EQ. 0) GOTO 60
      IF (F .LE. FOPT+TENTH*VQUAD) GOTO 60
C
C     Alternatively, find out if the interpolation points are close enough
C       to the best point so far.
C
      DISTSQ=DMAX1((TWO*DELTA)**2,(TEN*RHO)**2)
  650 KNEW=0
      DO 670 K=1,NPT
      SUM=ZERO
      DO 660 J=1,N
  660 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM .GT. DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  670 CONTINUE
C
C     If KNEW is positive, then ALTMOV finds alternative new positions for
C     the KNEW-th interpolation point within distance ADELT of XOPT. It is
C     reached via label 90. Otherwise, there is a branch to label 60 for
C     another trust region iteration, unless the calculations with the
C     current RHO are complete.
C
      IF (KNEW .GT. 0) THEN
          DIST=DSQRT(DISTSQ)
          IF (NTRITS .EQ. -1) THEN
              DELTA=DMIN1(TENTH*DELTA,HALF*DIST)
              IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
          END IF
          NTRITS=0
          ADELT=DMAX1(DMIN1(TENTH*DIST,DELTA),RHO)
          DSQ=ADELT*ADELT
          GOTO 90
      END IF
      IF (NTRITS .EQ. -1) GOTO 680
      IF (RATIO .GT. ZERO) GOTO 60
      IF (DMAX1(DELTA,DNORM) .GT. RHO) GOTO 60
C
C     The calculations with the current value of RHO are complete. Pick the
C       next values of RHO and DELTA.
C
  680 IF (RHO .GT. RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO .LE. 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO .LE. 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) WRITE(IPRINT,690)
  690         FORMAT (5X)
              WRITE(IPRINT,700) RHO,NF
  700         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              WRITE(IPRINT,710) FVAL(KOPT),(XBASE(I)+XOPT(I),I=1,N)
  710         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          NTRITS=0
          NFSAV=NF
          GOTO 60
      END IF
C
C     Return from the calculation, after another Newton-Raphson step, if
C       it is too short to have been tried before.
C
      IF (NTRITS .EQ. -1) GOTO 360
  720 IF (FVAL(KOPT) .LE. FSAVE) THEN
          DO 730 I=1,N
          X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XOPT(I)),XU(I))
          IF (XOPT(I) .EQ. SL(I)) X(I)=XL(I)
          IF (XOPT(I) .EQ. SU(I)) X(I)=XU(I)
  730     CONTINUE
          F=FVAL(KOPT)
      END IF
      IF (IPRINT .GE. 1) THEN
          WRITE(IPRINT,740) NF
  740     FORMAT (/4X,'At the return from BOBYQA',5X,
     1      'Number of function values =',I6)
          WRITE(IPRINT,710) F,(X(I),I=1,N)
      END IF
      RETURN
      END
      SUBROUTINE RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,
     1  FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,
     2  KOPT,VLAG,PTSAUX,PTSID,W,CALFUN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),XOPT(*),
     1  GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*),
     2  VLAG(*),PTSAUX(2,*),PTSID(*),W(*)
      EXTERNAL CALFUN
C
C     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
C       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
C       the corresponding arguments of BOBYQB on the entry to RESCUE.
C     NF is maintained as the number of calls of CALFUN so far, except that
C       NF is set to -1 if the value of MAXFUN prevents further progress.
C     KOPT is maintained so that FVAL(KOPT) is the least calculated function
C       value. Its correct value must be given on entry. It is updated if a
C       new least function value is found, but the corresponding changes to
C       XOPT and GOPT have to be made later by the calling program.
C     DELTA is the current trust region radius.
C     VLAG is a working space vector that will be used for the values of the
C       provisional Lagrange functions at each of the interpolation points.
C       They are part of a product that requires VLAG to be of length NDIM.
C     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
C       PTSAUX(2,J) specify the two positions of provisional interpolation
C       points when a nonzero step is taken along e_J (the J-th coordinate
C       direction) through XBASE+XOPT, as specified below. Usually these
C       steps have length DELTA, but other lengths are chosen if necessary
C       in order to satisfy the given bounds on the variables.
C     PTSID is also a working space array. It has NPT components that denote
C       provisional new positions of the original interpolation points, in
C       case changes are needed to restore the linear independence of the
C       interpolation conditions. The K-th point is a candidate for change
C       if and only if PTSID(K) is nonzero. In this case let p and q be the
C       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
C       and q are both positive, the step from XBASE+XOPT to the new K-th
C       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
C       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
C       p=0, respectively.
C     The first NDIM+NPT elements of the array W are used for working space. 
C     The final elements of BMAT and ZMAT are set in a well-conditioned way
C       to the values that are appropriate for the new interpolation points.
C     The elements of GOPT, HQ and PQ are also revised to the values that are
C       appropriate to the final quadratic model.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      NP=N+1
      SFRAC=HALF/DFLOAT(NP)
      NPTM=NPT-NP
C
C     Shift the interpolation points so that XOPT becomes the origin, and set
C     the elements of ZMAT to zero. The value of SUMPQ is required in the
C     updating of HQ below. The squares of the distances from XOPT to the
C     other interpolation points are set at the end of W. Increments of WINC
C     may be added later to these squares to balance the consideration of
C     the choice of point that is going to become current.
C
      SUMPQ=ZERO
      WINC=ZERO
      DO 20 K=1,NPT
      DISTSQ=ZERO
      DO 10 J=1,N
      XPT(K,J)=XPT(K,J)-XOPT(J)
   10 DISTSQ=DISTSQ+XPT(K,J)**2
      SUMPQ=SUMPQ+PQ(K)
      W(NDIM+K)=DISTSQ
      WINC=DMAX1(WINC,DISTSQ)
      DO 20 J=1,NPTM
   20 ZMAT(K,J)=ZERO
C
C     Update HQ so that HQ and PQ define the second derivatives of the model
C     after XBASE has been shifted to the trust region centre.
C
      IH=0
      DO 40 J=1,N
      W(J)=HALF*SUMPQ*XOPT(J)
      DO 30 K=1,NPT
   30 W(J)=W(J)+PQ(K)*XPT(K,J)
      DO 40 I=1,J
      IH=IH+1
   40 HQ(IH)=HQ(IH)+W(I)*XOPT(J)+W(J)*XOPT(I)
C
C     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
C     also set the elements of PTSAUX.
C
      DO 50 J=1,N
      XBASE(J)=XBASE(J)+XOPT(J)
      SL(J)=SL(J)-XOPT(J)
      SU(J)=SU(J)-XOPT(J)
      XOPT(J)=ZERO
      PTSAUX(1,J)=DMIN1(DELTA,SU(J))
      PTSAUX(2,J)=DMAX1(-DELTA,SL(J))
      IF (PTSAUX(1,J)+PTSAUX(2,J) .LT. ZERO) THEN
          TEMP=PTSAUX(1,J)
          PTSAUX(1,J)=PTSAUX(2,J)
          PTSAUX(2,J)=TEMP
      END IF
      IF (DABS(PTSAUX(2,J)) .LT. HALF*DABS(PTSAUX(1,J))) THEN
          PTSAUX(2,J)=HALF*PTSAUX(1,J)
      END IF
      DO 50 I=1,NDIM
   50 BMAT(I,J)=ZERO
      FBASE=FVAL(KOPT)
C
C     Set the identifiers of the artificial interpolation points that are
C     along a coordinate direction from XOPT, and set the corresponding
C     nonzero elements of BMAT and ZMAT.
C
      PTSID(1)=SFRAC
      DO 60 J=1,N
      JP=J+1
      JPN=JP+N
      PTSID(JP)=DFLOAT(J)+SFRAC
      IF (JPN .LE. NPT) THEN
          PTSID(JPN)=DFLOAT(J)/DFLOAT(NP)+SFRAC
          TEMP=ONE/(PTSAUX(1,J)-PTSAUX(2,J))
          BMAT(JP,J)=-TEMP+ONE/PTSAUX(1,J)
          BMAT(JPN,J)=TEMP+ONE/PTSAUX(2,J)
          BMAT(1,J)=-BMAT(JP,J)-BMAT(JPN,J)
          ZMAT(1,J)=DSQRT(2.0D0)/DABS(PTSAUX(1,J)*PTSAUX(2,J))
          ZMAT(JP,J)=ZMAT(1,J)*PTSAUX(2,J)*TEMP
          ZMAT(JPN,J)=-ZMAT(1,J)*PTSAUX(1,J)*TEMP
      ELSE
          BMAT(1,J)=-ONE/PTSAUX(1,J)
          BMAT(JP,J)=ONE/PTSAUX(1,J)
          BMAT(J+NPT,J)=-HALF*PTSAUX(1,J)**2
      END IF
   60 CONTINUE
C
C     Set any remaining identifiers with their nonzero elements of ZMAT.
C
      IF (NPT .GE. N+NP) THEN
          DO 70 K=2*NP,NPT
          IW=(DFLOAT(K-NP)-HALF)/DFLOAT(N)
          IP=K-NP-IW*N
          IQ=IP+IW
          IF (IQ .GT. N) IQ=IQ-N
          PTSID(K)=DFLOAT(IP)+DFLOAT(IQ)/DFLOAT(NP)+SFRAC
          TEMP=ONE/(PTSAUX(1,IP)*PTSAUX(1,IQ))
          ZMAT(1,K-NP)=TEMP
          ZMAT(IP+1,K-NP)=-TEMP
          ZMAT(IQ+1,K-NP)=-TEMP
   70     ZMAT(K,K-NP)=TEMP
      END IF
      NREM=NPT
      KOLD=1
      KNEW=KOPT
C
C     Reorder the provisional points in the way that exchanges PTSID(KOLD)
C     with PTSID(KNEW).
C
   80 DO 90 J=1,N
      TEMP=BMAT(KOLD,J)
      BMAT(KOLD,J)=BMAT(KNEW,J)
   90 BMAT(KNEW,J)=TEMP
      DO 100 J=1,NPTM
      TEMP=ZMAT(KOLD,J)
      ZMAT(KOLD,J)=ZMAT(KNEW,J)
  100 ZMAT(KNEW,J)=TEMP
      PTSID(KOLD)=PTSID(KNEW)
      PTSID(KNEW)=ZERO
      W(NDIM+KNEW)=ZERO
      NREM=NREM-1
      IF (KNEW .NE. KOPT) THEN
          TEMP=VLAG(KOLD)
          VLAG(KOLD)=VLAG(KNEW)
          VLAG(KNEW)=TEMP
C
C     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
C     interpolation point can be changed from provisional to original. The
C     branch to label 350 occurs if all the original points are reinstated.
C     The nonnegative values of W(NDIM+K) are required in the search below.
C
          CALL BOBUPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
          IF (NREM .EQ. 0) GOTO 350
          DO 110 K=1,NPT
  110     W(NDIM+K)=DABS(W(NDIM+K))
      END IF
C
C     Pick the index KNEW of an original interpolation point that has not
C     yet replaced one of the provisional interpolation points, giving
C     attention to the closeness to XOPT and to previous tries with KNEW.
C
  120 DSQMIN=ZERO
      DO 130 K=1,NPT
      IF (W(NDIM+K) .GT. ZERO) THEN
          IF (DSQMIN .EQ. ZERO .OR. W(NDIM+K) .LT. DSQMIN) THEN
              KNEW=K
              DSQMIN=W(NDIM+K)
          END IF
      END IF
  130 CONTINUE
      IF (DSQMIN .EQ. ZERO) GOTO 260
C
C     Form the W-vector of the chosen original interpolation point.
C
      DO 140 J=1,N
  140 W(NPT+J)=XPT(KNEW,J)
      DO 160 K=1,NPT
      SUM=ZERO
      IF (K .EQ. KOPT) THEN
          CONTINUE
      ELSE IF (PTSID(K) .EQ. ZERO) THEN
          DO 150 J=1,N
  150     SUM=SUM+W(NPT+J)*XPT(K,J)
      ELSE
          IP=PTSID(K)
          IF (IP .GT. 0) SUM=W(NPT+IP)*PTSAUX(1,IP)
          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
          IF (IQ .GT. 0) THEN
              IW=1
              IF (IP .EQ. 0) IW=2
              SUM=SUM+W(NPT+IQ)*PTSAUX(IW,IQ)
          END IF
      END IF
  160 W(K)=HALF*SUM*SUM
C
C     Calculate VLAG and BETA for the required updating of the H matrix if
C     XPT(KNEW,.) is reinstated in the set of interpolation points.
C
      DO 180 K=1,NPT
      SUM=ZERO
      DO 170 J=1,N
  170 SUM=SUM+BMAT(K,J)*W(NPT+J)
  180 VLAG(K)=SUM
      BETA=ZERO
      DO 200 J=1,NPTM
      SUM=ZERO
      DO 190 K=1,NPT
  190 SUM=SUM+ZMAT(K,J)*W(K)
      BETA=BETA-SUM*SUM
      DO 200 K=1,NPT
  200 VLAG(K)=VLAG(K)+SUM*ZMAT(K,J)
      BSUM=ZERO
      DISTSQ=ZERO
      DO 230 J=1,N
      SUM=ZERO
      DO 210 K=1,NPT
  210 SUM=SUM+BMAT(K,J)*W(K)
      JP=J+NPT
      BSUM=BSUM+SUM*W(JP)
      DO 220 IP=NPT+1,NDIM
  220 SUM=SUM+BMAT(IP,J)*W(IP)
      BSUM=BSUM+SUM*W(JP)
      VLAG(JP)=SUM
  230 DISTSQ=DISTSQ+XPT(KNEW,J)**2
      BETA=HALF*DISTSQ*DISTSQ+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     KOLD is set to the index of the provisional interpolation point that is
C     going to be deleted to make way for the KNEW-th original interpolation
C     point. The choice of KOLD is governed by the avoidance of a small value
C     of the denominator in the updating calculation of UPDATE.
C
      DENOM=ZERO
      VLMXSQ=ZERO
      DO 250 K=1,NPT
      IF (PTSID(K) .NE. ZERO) THEN
          HDIAG=ZERO
          DO 240 J=1,NPTM
  240     HDIAG=HDIAG+ZMAT(K,J)**2
          DEN=BETA*HDIAG+VLAG(K)**2
          IF (DEN .GT. DENOM) THEN
              KOLD=K
              DENOM=DEN
          END IF
      END IF
  250 VLMXSQ=DMAX1(VLMXSQ,VLAG(K)**2)
      IF (DENOM .LE. 1.0D-2*VLMXSQ) THEN
          W(NDIM+KNEW)=-W(NDIM+KNEW)-WINC
          GOTO 120
      END IF
      GOTO 80
C
C     When label 260 is reached, all the final positions of the interpolation
C     points have been chosen although any changes have not been included yet
C     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
C     from the shift of XBASE, the updating of the quadratic model remains to
C     be done. The following cycle through the new interpolation points begins
C     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
C     except that a RETURN occurs if MAXFUN prohibits another value of F.
C
  260 DO 340 KPT=1,NPT
      IF (PTSID(KPT) .EQ. ZERO) GOTO 340
      IF (NF .GE. MAXFUN) THEN
          NF=-1
          GOTO 350
      END IF
      IH=0
      DO 270 J=1,N
      W(J)=XPT(KPT,J)
      XPT(KPT,J)=ZERO
      TEMP=PQ(KPT)*W(J)
      DO 270 I=1,J
      IH=IH+1
  270 HQ(IH)=HQ(IH)+TEMP*W(I)
      PQ(KPT)=ZERO
      IP=PTSID(KPT)
      IQ=DFLOAT(NP)*PTSID(KPT)-DFLOAT(IP*NP)
      IF (IP .GT. 0) THEN
          XP=PTSAUX(1,IP)
          XPT(KPT,IP)=XP
      END IF
      IF (IQ .GT. 0) THEN
          XQ=PTSAUX(1,IQ)
          IF (IP .EQ. 0) XQ=PTSAUX(2,IQ)
          XPT(KPT,IQ)=XQ
      END IF
C
C     Set VQUAD to the value of the current model at the new point.
C
      VQUAD=FBASE
      IF (IP .GT. 0) THEN
          IHP=(IP+IP*IP)/2
          VQUAD=VQUAD+XP*(GOPT(IP)+HALF*XP*HQ(IHP))
      END IF
      IF (IQ .GT. 0) THEN
          IHQ=(IQ+IQ*IQ)/2
          VQUAD=VQUAD+XQ*(GOPT(IQ)+HALF*XQ*HQ(IHQ))
          IF (IP .GT. 0) THEN
              IW=MAX0(IHP,IHQ)-IABS(IP-IQ)
              VQUAD=VQUAD+XP*XQ*HQ(IW)
          END IF
      END IF
      DO 280 K=1,NPT
      TEMP=ZERO
      IF (IP .GT. 0) TEMP=TEMP+XP*XPT(K,IP)
      IF (IQ .GT. 0) TEMP=TEMP+XQ*XPT(K,IQ)
  280 VQUAD=VQUAD+HALF*PQ(K)*TEMP*TEMP
C
C     Calculate F at the new interpolation point, and set DIFF to the factor
C     that is going to multiply the KPT-th Lagrange function when the model
C     is updated to provide interpolation to the new function value.
C
      DO 290 I=1,N
      W(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XPT(KPT,I)),XU(I))
      IF (XPT(KPT,I) .EQ. SL(I)) W(I)=XL(I)
      IF (XPT(KPT,I) .EQ. SU(I)) W(I)=XU(I)
  290 CONTINUE
      NF=NF+1
      CALL CALFUN (N,W,F)
      IF (IPRINT .EQ. 3) THEN
          WRITE(IPRINT,300) NF,F,(W(I),I=1,N)
  300     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1      '    The corresponding X is:'/(2X,5D15.6))
      END IF
      FVAL(KPT)=F
      IF (F .LT. FVAL(KOPT)) KOPT=KPT
      DIFF=F-VQUAD
C
C     Update the quadratic model. The RETURN from the subroutine occurs when
C     all the new interpolation points are included in the model.
C
      DO 310 I=1,N
  310 GOPT(I)=GOPT(I)+DIFF*BMAT(KPT,I)
      DO 330 K=1,NPT
      SUM=ZERO
      DO 320 J=1,NPTM
  320 SUM=SUM+ZMAT(K,J)*ZMAT(KPT,J)
      TEMP=DIFF*SUM
      IF (PTSID(K) .EQ. ZERO) THEN
          PQ(K)=PQ(K)+TEMP
      ELSE
          IP=PTSID(K)
          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
          IHQ=(IQ*IQ+IQ)/2
          IF (IP .EQ. 0) THEN
              HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(2,IQ)**2
          ELSE
              IHP=(IP*IP+IP)/2
              HQ(IHP)=HQ(IHP)+TEMP*PTSAUX(1,IP)**2
              IF (IQ .GT. 0) THEN
                  HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(1,IQ)**2
                  IW=MAX0(IHP,IHQ)-IABS(IQ-IP)
                  HQ(IW)=HQ(IW)+TEMP*PTSAUX(1,IP)*PTSAUX(1,IQ)
              END IF
          END IF
      END IF
  330 CONTINUE
      PTSID(KPT)=ZERO
  340 CONTINUE
  350 RETURN
      END
      SUBROUTINE TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,
     1  XNEW,D,GNEW,XBDI,S,HS,HRED,DSQ,CRVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPT(NPT,*),XOPT(*),GOPT(*),HQ(*),PQ(*),SL(*),SU(*),
     1  XNEW(*),D(*),GNEW(*),XBDI(*),S(*),HS(*),HRED(*)
C
C     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
C       meanings as the corresponding arguments of BOBYQB.
C     DELTA is the trust region radius for the present calculation, which
C       seeks a small value of the quadratic model within distance DELTA of
C       XOPT subject to the bounds on the variables.
C     XNEW will be set to a new vector of variables that is approximately
C       the one that minimizes the quadratic model within the trust region
C       subject to the SL and SU constraints on the variables. It satisfies
C       as equations the bounds that become active during the calculation.
C     D is the calculated trial step from XOPT, generated iteratively from an
C       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
C     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
C       when D is updated.
C     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
C       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
C       I-th variable has become fixed at a bound, the bound being SL(I) or
C       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
C       information is accumulated during the construction of XNEW.
C     The arrays S, HS and HRED are also used for working space. They hold the
C       current search direction, and the changes in the gradient of Q along S
C       and the reduced D, respectively, where the reduced D is the same as D,
C       except that the components of the fixed variables are zero.
C     DSQ will be set to the square of the length of XNEW-XOPT.
C     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
C       it is set to the least curvature of H that occurs in the conjugate
C       gradient searches that are not restricted by any constraints. The
C       value CRVMIN=-1.0D0 is set, however, if all of these searches are
C       constrained.
C
C     A version of the truncated conjugate gradient is applied. If a line
C     search is restricted by a constraint, then the procedure is restarted,
C     the values of the variables that are at their bounds being fixed. If
C     the trust region boundary is reached, then further changes may be made
C     to D, each one being in the two dimensional space that is spanned
C     by the current D and the gradient of Q at XOPT+D, staying on the trust
C     region boundary. Termination occurs when the reduction in Q seems to
C     be close to the greatest reduction that can be achieved.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ONEMIN=-1.0D0
      ZERO=0.0D0
C
C     The sign of GOPT(I) gives the sign of the change to the I-th variable
C     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
C     or not to fix the I-th variable at one of its bounds initially, with
C     NACT being set to the number of fixed variables. D and GNEW are also
C     set for the first iteration. DELSQ is the upper bound on the sum of
C     squares of the free variables. QRED is the reduction in Q so far.
C
      ITERC=0
      NACT=0
      SQSTP=ZERO
      DO 10 I=1,N
      XBDI(I)=ZERO
      IF (XOPT(I) .LE. SL(I)) THEN
          IF (GOPT(I) .GE. ZERO) XBDI(I)=ONEMIN
      ELSE IF (XOPT(I) .GE. SU(I)) THEN
          IF (GOPT(I) .LE. ZERO) XBDI(I)=ONE
      END IF
      IF (XBDI(I) .NE. ZERO) NACT=NACT+1
      D(I)=ZERO
   10 GNEW(I)=GOPT(I)
      DELSQ=DELTA*DELTA
      QRED=ZERO
      CRVMIN=ONEMIN
C
C     Set the next search direction of the conjugate gradient method. It is
C     the steepest descent direction initially and when the iterations are
C     restarted because a variable has just been fixed by a bound, and of
C     course the components of the fixed variables are zero. ITERMAX is an
C     upper bound on the indices of the conjugate gradient iterations.
C
   20 BETA=ZERO
   30 STEPSQ=ZERO
      DO 40 I=1,N
      IF (XBDI(I) .NE. ZERO) THEN
          S(I)=ZERO
      ELSE IF (BETA .EQ. ZERO) THEN
          S(I)=-GNEW(I)
      ELSE
          S(I)=BETA*S(I)-GNEW(I)
      END IF
   40 STEPSQ=STEPSQ+S(I)**2
      IF (STEPSQ .EQ. ZERO) GOTO 190
      IF (BETA .EQ. ZERO) THEN
          GREDSQ=STEPSQ
          ITERMAX=ITERC+N-NACT
      END IF
      IF (GREDSQ*DELSQ .LE. 1.0D-4*QRED*QRED) GO TO 190
C
C     Multiply the search direction by the second derivative matrix of Q and
C     calculate some scalars for the choice of steplength. Then set BLEN to
C     the length of the the step to the trust region boundary and STPLEN to
C     the steplength, ignoring the simple bounds.
C
      GOTO 210
   50 RESID=DELSQ
      DS=ZERO
      SHS=ZERO
      DO 60 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          RESID=RESID-D(I)**2
          DS=DS+S(I)*D(I)
          SHS=SHS+S(I)*HS(I)
      END IF
   60 CONTINUE
      IF (RESID .LE. ZERO) GOTO 90
      TEMP=DSQRT(STEPSQ*RESID+DS*DS)
      IF (DS .LT. ZERO) THEN
          BLEN=(TEMP-DS)/STEPSQ
      ELSE
          BLEN=RESID/(TEMP+DS)
      END IF
      STPLEN=BLEN
      IF (SHS .GT. ZERO) THEN
          STPLEN=DMIN1(BLEN,GREDSQ/SHS)
      END IF
      
C
C     Reduce STPLEN if necessary in order to preserve the simple bounds,
C     letting IACT be the index of the new constrained variable.
C
      IACT=0
      DO 70 I=1,N
      IF (S(I) .NE. ZERO) THEN
          XSUM=XOPT(I)+D(I)
          IF (S(I) .GT. ZERO) THEN
              TEMP=(SU(I)-XSUM)/S(I)
          ELSE
              TEMP=(SL(I)-XSUM)/S(I)
          END IF
          IF (TEMP .LT. STPLEN) THEN
              STPLEN=TEMP
              IACT=I
          END IF
      END IF
   70 CONTINUE
C
C     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
C
      SDEC=ZERO
      IF (STPLEN .GT. ZERO) THEN
          ITERC=ITERC+1
          TEMP=SHS/STEPSQ
          IF (IACT .EQ. 0 .AND. TEMP .GT. ZERO) THEN
              CRVMIN=DMIN1(CRVMIN,TEMP)
              IF (CRVMIN .EQ. ONEMIN) CRVMIN=TEMP
          END IF 
          GGSAV=GREDSQ
          GREDSQ=ZERO
          DO 80 I=1,N
          GNEW(I)=GNEW(I)+STPLEN*HS(I)
          IF (XBDI(I) .EQ. ZERO) GREDSQ=GREDSQ+GNEW(I)**2
   80     D(I)=D(I)+STPLEN*S(I)
          SDEC=DMAX1(STPLEN*(GGSAV-HALF*STPLEN*SHS),ZERO)
          QRED=QRED+SDEC
      END IF
C
C     Restart the conjugate gradient method if it has hit a new bound.
C
      IF (IACT .GT. 0) THEN
          NACT=NACT+1
          XBDI(IACT)=ONE
          IF (S(IACT) .LT. ZERO) XBDI(IACT)=ONEMIN
          DELSQ=DELSQ-D(IACT)**2
          IF (DELSQ .LE. ZERO) GOTO 90
          GOTO 20
      END IF
C
C     If STPLEN is less than BLEN, then either apply another conjugate
C     gradient iteration or RETURN.
C
      IF (STPLEN .LT. BLEN) THEN
          IF (ITERC .EQ. ITERMAX) GOTO 190
          IF (SDEC .LE. 0.01D0*QRED) GOTO 190
          BETA=GREDSQ/GGSAV
          GOTO 30
      END IF
   90 CRVMIN=ZERO
C
C     Prepare for the alternative iteration by calculating some scalars and
C     by multiplying the reduced D by the second derivative matrix of Q.
C
  100 IF (NACT .GE. N-1) GOTO 190
      DREDSQ=ZERO
      DREDG=ZERO
      GREDSQ=ZERO
      DO 110 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          DREDSQ=DREDSQ+D(I)**2
          DREDG=DREDG+D(I)*GNEW(I)
          GREDSQ=GREDSQ+GNEW(I)**2
          S(I)=D(I)
      ELSE
          S(I)=ZERO
      END IF
  110 CONTINUE
      ITCSAV=ITERC
      GOTO 210
C
C     Let the search direction S be a linear combination of the reduced D
C     and the reduced G that is orthogonal to the reduced D.
C
  120 ITERC=ITERC+1
      TEMP=GREDSQ*DREDSQ-DREDG*DREDG
      IF (TEMP .LE. 1.0D-4*QRED*QRED) GOTO 190
      TEMP=DSQRT(TEMP)
      DO 130 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          S(I)=(DREDG*D(I)-DREDSQ*GNEW(I))/TEMP
      ELSE
          S(I)=ZERO
      END IF
  130 CONTINUE
      SREDG=-TEMP
C
C     By considering the simple bounds on the variables, calculate an upper
C     bound on the tangent of half the angle of the alternative iteration,
C     namely ANGBD, except that, if already a free variable has reached a
C     bound, there is a branch back to label 100 after fixing that variable.
C
      ANGBD=ONE
      IACT=0
      DO 140 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          TEMPA=XOPT(I)+D(I)-SL(I)
          TEMPB=SU(I)-XOPT(I)-D(I)
          IF (TEMPA .LE. ZERO) THEN
              NACT=NACT+1
              XBDI(I)=ONEMIN
              GOTO 100
          ELSE IF (TEMPB .LE. ZERO) THEN
              NACT=NACT+1
              XBDI(I)=ONE
              GOTO 100
          END IF
          RATIO=ONE
          SSQ=D(I)**2+S(I)**2
          TEMP=SSQ-(XOPT(I)-SL(I))**2
          IF (TEMP .GT. ZERO) THEN
              TEMP=DSQRT(TEMP)-S(I)
              IF (ANGBD*TEMP .GT. TEMPA) THEN
                  ANGBD=TEMPA/TEMP
                  IACT=I
                  XSAV=ONEMIN
              END IF
          END IF
          TEMP=SSQ-(SU(I)-XOPT(I))**2
          IF (TEMP .GT. ZERO) THEN
              TEMP=DSQRT(TEMP)+S(I)
              IF (ANGBD*TEMP .GT. TEMPB) THEN
                  ANGBD=TEMPB/TEMP
                  IACT=I
                  XSAV=ONE
              END IF
          END IF
      END IF
  140 CONTINUE
C
C     Calculate HHD and some curvatures for the alternative iteration.
C
      GOTO 210
  150 SHS=ZERO
      DHS=ZERO
      DHD=ZERO
      DO 160 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          SHS=SHS+S(I)*HS(I)
          DHS=DHS+D(I)*HS(I)
          DHD=DHD+D(I)*HRED(I)
      END IF
  160 CONTINUE
C
C     Seek the greatest reduction in Q for a range of equally spaced values
C     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
C     the alternative iteration.
C
      REDMAX=ZERO
      ISAV=0
      REDSAV=ZERO
      IU=17.0D0*ANGBD+3.1D0
      DO 170 I=1,IU
      ANGT=ANGBD*DFLOAT(I)/DFLOAT(IU)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      REDNEW=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      IF (REDNEW .GT. REDMAX) THEN
          REDMAX=REDNEW
          ISAV=I
          RDPREV=REDSAV
      ELSE IF (I .EQ. ISAV+1) THEN
          RDNEXT=REDNEW
      END IF
  170 REDSAV=REDNEW
C
C     Return if the reduction is zero. Otherwise, set the sine and cosine
C     of the angle of the alternative iteration, and calculate SDEC.
C
      IF (ISAV .EQ. 0) GOTO 190
      IF (ISAV .LT. IU) THEN
          TEMP=(RDNEXT-RDPREV)/(REDMAX+REDMAX-RDPREV-RDNEXT)
          ANGT=ANGBD*(DFLOAT(ISAV)+HALF*TEMP)/DFLOAT(IU)
      END IF
      CTH=(ONE-ANGT*ANGT)/(ONE+ANGT*ANGT)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      SDEC=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      IF (SDEC .LE. ZERO) GOTO 190
C
C     Update GNEW, D and HRED. If the angle of the alternative iteration
C     is restricted by a bound on a free variable, that variable is fixed
C     at the bound.
C
      DREDG=ZERO
      GREDSQ=ZERO
      DO 180 I=1,N
      GNEW(I)=GNEW(I)+(CTH-ONE)*HRED(I)+STH*HS(I)
      IF (XBDI(I) .EQ. ZERO) THEN
          D(I)=CTH*D(I)+STH*S(I)
          DREDG=DREDG+D(I)*GNEW(I)
          GREDSQ=GREDSQ+GNEW(I)**2
      END IF
  180 HRED(I)=CTH*HRED(I)+STH*HS(I)
      QRED=QRED+SDEC
      IF (IACT .GT. 0 .AND. ISAV .EQ. IU) THEN
          NACT=NACT+1
          XBDI(IACT)=XSAV
          GOTO 100
      END IF
C
C     If SDEC is sufficiently small, then RETURN after setting XNEW to
C     XOPT+D, giving careful attention to the bounds.
C
      IF (SDEC .GT. 0.01D0*QRED) GOTO 120
  190 DSQ=ZERO
      DO 200 I=1,N
      XNEW(I)=DMAX1(DMIN1(XOPT(I)+D(I),SU(I)),SL(I))
      IF (XBDI(I) .EQ. ONEMIN) XNEW(I)=SL(I)
      IF (XBDI(I) .EQ. ONE) XNEW(I)=SU(I)
      D(I)=XNEW(I)-XOPT(I)
  200 DSQ=DSQ+D(I)**2
      RETURN
 
C     The following instructions multiply the current S-vector by the second
C     derivative matrix of the quadratic model, putting the product in HS.
C     They are reached from three different parts of the software above and
C     they can be regarded as an external subroutine.
C
  210 IH=0
      DO 220 J=1,N
      HS(J)=ZERO
      DO 220 I=1,J
      IH=IH+1
      IF (I .LT. J) HS(J)=HS(J)+HQ(IH)*S(I)
  220 HS(I)=HS(I)+HQ(IH)*S(J)
      DO 250 K=1,NPT
      IF (PQ(K) .NE. ZERO) THEN
          TEMP=ZERO
          DO 230 J=1,N
  230     TEMP=TEMP+XPT(K,J)*S(J)
          TEMP=TEMP*PQ(K)
          DO 240 I=1,N
  240     HS(I)=HS(I)+TEMP*XPT(K,I)
      END IF
  250 CONTINUE
      IF (CRVMIN .NE. ZERO) GOTO 50
      IF (ITERC .GT. ITCSAV) GOTO 150
      DO 260 I=1,N
  260 HRED(I)=HS(I)
      GOTO 120
      END
      SUBROUTINE BOBUPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,
     1  KNEW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
C
C     The arrays BMAT and ZMAT are updated, as required by the new position
C     of the interpolation point that has the index KNEW. The vector VLAG has
C     N+NPT components, set on entry to the first NPT and last N components
C     of the product Hw in equation (4.11) of the Powell (2006) paper on
C     NEWUOA. Further, BETA is set on entry to the value of the parameter
C     with that name, and DENOM is set to the denominator of the updating
C     formula. Elements of ZMAT may be treated as zero if their moduli are
C     at most ZTEST. The first NDIM elements of W are used for working space.
C
C     Set some constants.
C
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
      ZTEST=ZERO
      DO 10 K=1,NPT
      DO 10 J=1,NPTM
   10 ZTEST=DMAX1(ZTEST,DABS(ZMAT(K,J)))
      ZTEST=1.0D-20*ZTEST
C
C     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
C
      JL=1
      DO 30 J=2,NPTM
      IF (DABS(ZMAT(KNEW,J)) .GT. ZTEST) THEN
          TEMP=DSQRT(ZMAT(KNEW,1)**2+ZMAT(KNEW,J)**2)
          TEMPA=ZMAT(KNEW,1)/TEMP
          TEMPB=ZMAT(KNEW,J)/TEMP
          DO 20 I=1,NPT
          TEMP=TEMPA*ZMAT(I,1)+TEMPB*ZMAT(I,J)
          ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,1)
   20     ZMAT(I,1)=TEMP
      END IF
      ZMAT(KNEW,J)=ZERO
   30 CONTINUE
C
C     Put the first NPT components of the KNEW-th column of HLAG into W,
C     and calculate the parameters of the updating formula.
C
      DO 40 I=1,NPT
      W(I)=ZMAT(KNEW,1)*ZMAT(I,1)
   40 CONTINUE
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      VLAG(KNEW)=VLAG(KNEW)-ONE
C
C     Complete the updating of ZMAT.
C
      TEMP=DSQRT(DENOM)
      TEMPB=ZMAT(KNEW,1)/TEMP
      TEMPA=TAU/TEMP
      DO 50 I=1,NPT
   50 ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
C
C     Finally, update the matrix BMAT.
C
      DO 60 J=1,N
      JP=NPT+J
      W(JP)=BMAT(KNEW,J)
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
      DO 60 I=1,JP
      BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
      IF (I .GT. NPT) BMAT(JP,I-NPT)=BMAT(I,J)
   60 CONTINUE
      RETURN
      END
      SUBROUTINE BOPRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,
     1  XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT,CALFUN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),GOPT(*),
     1  HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*)
      EXTERNAL CALFUN
C
C     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
C       same as the corresponding arguments in SUBROUTINE BOBYQA.
C     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
C       are the same as the corresponding arguments in BOBYQB, the elements
C       of SL and SU being set in BOBYQA.
C     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
C       it is set by BOPRELIM to the gradient of the quadratic model at XBASE.
C       If XOPT is nonzero, BOBYQB will change it to its usual value later.
C     NF is maintaned as the number of calls of CALFUN so far.
C     KOPT will be such that the least calculated value of F so far is at
C       the point XPT(KOPT,.)+XBASE in the space of the variables.
C
C     SUBROUTINE BOPRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
C     BMAT and ZMAT for the first iteration, and it maintains the values of
C     NF and KOPT. The vector X is also changed by BOPRELIM.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TWO=2.0D0
      ZERO=0.0D0
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      NP=N+1
C
C     Set XBASE to the initial vector of variables, and set the initial
C     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
C
      DO 20 J=1,N
      XBASE(J)=X(J)
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 IH=1,(N*NP)/2
   30 HQ(IH)=ZERO
      DO 40 K=1,NPT
      PQ(K)=ZERO
      DO 40 J=1,NPT-NP
   40 ZMAT(K,J)=ZERO
C
C     Begin the initialization procedure. NF becomes one more than the number
C     of function values so far. The coordinates of the displacement of the
C     next initial interpolation point from XBASE are set in XPT(NF+1,.).
C
      NF=0
   50 NFM=NF
      NFX=NF-N
      NF=NF+1
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              STEPA=RHOBEG
              IF (SU(NFM) .EQ. ZERO) STEPA=-STEPA
              XPT(NF,NFM)=STEPA
          ELSE IF (NFM .GT. N) THEN
              STEPA=XPT(NF-N,NFX)
              STEPB=-RHOBEG
              IF (SL(NFX) .EQ. ZERO) STEPB=DMIN1(TWO*RHOBEG,SU(NFX))
              IF (SU(NFX) .EQ. ZERO) STEPB=DMAX1(-TWO*RHOBEG,SL(NFX))
              XPT(NF,NFX)=STEPB
          END IF
      ELSE
          ITEMP=(NFM-NP)/N
          JPT=NFM-ITEMP*N-N
          IPT=JPT+ITEMP
          IF (IPT .GT. N) THEN
              ITEMP=JPT
              JPT=IPT-N
              IPT=ITEMP
          END IF
          XPT(NF,IPT)=XPT(IPT+1,IPT)
          XPT(NF,JPT)=XPT(JPT+1,JPT)
      END IF
C
C     Calculate the next value of F. The least function value so far and
C     its index are required.
C
      DO 60 J=1,N
      X(J)=DMIN1(DMAX1(XL(J),XBASE(J)+XPT(NF,J)),XU(J))
      IF (XPT(NF,J) .EQ. SL(J)) X(J)=XL(J)
      IF (XPT(NF,J) .EQ. SU(J)) X(J)=XU(J)
   60 CONTINUE
      CALL CALFUN (N,X,F)
      IF (IPRINT .GE. 3) THEN
          WRITE(IPRINT,70) NF,F,(X(I),I=1,N)
   70     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1       '    The corresponding X is:'/(2X,5D15.6))
      END IF
      FVAL(NF)=F
      IF (NF .EQ. 1) THEN
          FBEG=F
          KOPT=1
      ELSE IF (F .LT. FVAL(KOPT)) THEN
          KOPT=NF
      END IF
C
C     Set the nonzero initial elements of BMAT and the quadratic model in the
C     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
C     of the NF-th and (NF-N)-th interpolation points may be switched, in
C     order that the function value at the first of them contributes to the
C     off-diagonal second derivative terms of the initial quadratic model.
C
      IF (NF .LE. 2*N+1) THEN
          IF (NF .GE. 2 .AND. NF .LE. N+1) THEN
              GOPT(NFM)=(F-FBEG)/STEPA
              IF (NPT .LT. NF+N) THEN
                  BMAT(1,NFM)=-ONE/STEPA
                  BMAT(NF,NFM)=ONE/STEPA
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              END IF
          ELSE IF (NF .GE. N+2) THEN
              IH=(NFX*(NFX+1))/2
              TEMP=(F-FBEG)/STEPB
              DIFF=STEPB-STEPA
              HQ(IH)=TWO*(TEMP-GOPT(NFX))/DIFF
              GOPT(NFX)=(GOPT(NFX)*STEPB-TEMP*STEPA)/DIFF
              IF (STEPA*STEPB .LT. ZERO) THEN
                  IF (F .LT. FVAL(NF-N)) THEN
                      FVAL(NF)=FVAL(NF-N)
                      FVAL(NF-N)=F
                      IF (KOPT .EQ. NF) KOPT=NF-N
                      XPT(NF-N,NFX)=STEPB
                      XPT(NF,NFX)=STEPA
                  END IF
              END IF
              BMAT(1,NFX)=-(STEPA+STEPB)/(STEPA*STEPB)
              BMAT(NF,NFX)=-HALF/XPT(NF-N,NFX)
              BMAT(NF-N,NFX)=-BMAT(1,NFX)-BMAT(NF,NFX)
              ZMAT(1,NFX)=DSQRT(TWO)/(STEPA*STEPB)
              ZMAT(NF,NFX)=DSQRT(HALF)/RHOSQ
              ZMAT(NF-N,NFX)=-ZMAT(1,NFX)-ZMAT(NF,NFX)
          END IF
C
C     Set the off-diagonal second derivatives of the Lagrange functions and
C     the initial quadratic model.
C
      ELSE
          IH=(IPT*(IPT-1))/2+JPT
          ZMAT(1,NFX)=RECIP
          ZMAT(NF,NFX)=RECIP
          ZMAT(IPT+1,NFX)=-RECIP
          ZMAT(JPT+1,NFX)=-RECIP
          TEMP=XPT(NF,IPT)*XPT(NF,JPT)
          HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/TEMP
      END IF
      IF (NF .LT. NPT .AND. NF .LT. MAXFUN) GOTO 50
      RETURN
      END
      SUBROUTINE LINCOA (N,NPT,M,A,IA,B,X,RHOBEG,RHOEND,LSTDOUT,
     1  MAXFUN,W,CALFUN,OUTFILE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,*),B(*),X(*),W(*)
C     Definition of stream control variables 256 appears in calling codes
      EXTERNAL      CALFUN
      INTEGER       LSTDOUT
      CHARACTER*256 OUTFILE
      IF (LSTDOUT .gt. 1) THEN
          IPRINT = 12
          OPEN (IPRINT, FILE=OUTFILE, STATUS='OLD')
      ELSE
          IPRINT = 0
      ENDIF
C     End of it. KMK
C
C     This subroutine seeks the least value of a function of many variables,
C       subject to general linear inequality constraints, by a trust region
C       method that forms quadratic models by interpolation. Usually there
C       is much freedom in each new model after satisfying the interpolation
C       conditions, which is taken up by minimizing the Frobenius norm of
C       the change to the second derivative matrix of the model. One new
C       function value is calculated on each iteration, usually at a point
C       where the current model predicts a reduction in the least value so
C       far of the objective function subject to the linear constraints.
C       Alternatively, a new vector of variables may be chosen to replace
C       an interpolation point that may be too far away for reliability, and
C       then the new point does not have to satisfy the linear constraints.
C       The arguments of the subroutine are as follows.
C
C     N must be set to the number of variables and must be at least two.
C     NPT must be set to the number of interpolation conditions, which is
C       required to be in the interval [N+2,(N+1)(N+2)/2]. Typical choices
C       of the author are NPT=N+6 and NPT=2*N+1. Larger values tend to be
C       highly inefficent when the number of variables is substantial, due
C       to the amount of work and extra difficulty of adjusting more points.
C     M must be set to the number of linear inequality constraints.
C     A is a matrix whose columns are the constraint gradients, which are
C       required to be nonzero.
C     IA is the first dimension of the array A, which must be at least N.
C     B is the vector of right hand sides of the constraints, the J-th
C       constraint being that the scalar product of A(.,J) with X(.) is at
C       most B(J). The initial vector X(.) is made feasible by increasing
C       the value of B(J) if necessary.
C     X is the vector of variables. Initial values of X(1),X(2),...,X(N)
C       must be supplied. If they do not satisfy the constraints, then B
C       is increased as mentioned above. X contains on return the variables
C       that have given the least calculated F subject to the constraints.
C     RHOBEG and RHOEND must be set to the initial and final values of a
C       trust region radius, so both must be positive with RHOEND<=RHOBEG.
C       Typically, RHOBEG should be about one tenth of the greatest expected
C       change to a variable, and RHOEND should indicate the accuracy that
C       is required in the final values of the variables.
C     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
C       amount of printing. Specifically, there is no output if IPRINT=0 and
C       there is output only at the return if IPRINT=1. Otherwise, the best
C       feasible vector of variables so far and the corresponding value of
C       the objective function are printed whenever RHO is reduced, where
C       RHO is the current lower bound on the trust region radius. Further,
C       each new value of F with its variables are output if IPRINT=3.
C     MAXFUN must be set to an upper bound on the number of calls of CALFUN,
C       its value being at least NPT+1.
C     W is an array used for working space. Its length must be at least
C       M*(2+N) + NPT*(4+N+NPT) + N*(9+3*N) + MAX [ M+3*N, 2*M+N, 2*NPT ].
C       On return, W(1) is set to the final value of F, and W(2) is set to
C       the total number of function evaluations plus 0.5.
C
C     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
C       F to the value of the objective function for the variables X(1),
C       X(2),...,X(N). The value of the argument F is positive when CALFUN
C       is called if and only if the current X satisfies the constraints
C       to working accuracy.
C
C     Check that N, NPT and MAXFUN are acceptable.
C
      ZERO=0.0D0
      SMALLX=1.0D-6*RHOEND
      NP=N+1
      NPTM=NPT-NP
      IF (N .LE. 1) THEN
          IF (IPRINT .GT. 0) WRITE(IPRINT,10)
   10     FORMAT (/4X,'Return from LINCOA because N is less than 2.')
          GOTO 80
      END IF
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
          IF (IPRINT .GT. 0) WRITE(IPRINT,20)
   20     FORMAT (/4X,'Return from LINCOA because NPT is not in',
     1      ' the required interval.')
          GOTO 80
      END IF
      IF (MAXFUN .LE. NPT) THEN
          IF (IPRINT .GT. 0) WRITE(IPRINT,30)
   30     FORMAT (/4X,'Return from LINCOA because MAXFUN is less',
     1      ' than NPT+1.')
          GOTO 80
      END IF
C
C     Normalize the constraints, and copy the resultant constraint matrix
C       and right hand sides into working space, after increasing the right
C       hand sides if necessary so that the starting point is feasible.
C
      IAMAT=MAX0(M+3*N,2*M+N,2*NPT)+1
      IB=IAMAT+M*N
      IFLAG=0
      IF (M .GT. 0) THEN
          IW=IAMAT-1
          DO 60 J=1,M
          SUM=ZERO
          TEMP=ZERO
          DO 40 I=1,N
          SUM=SUM+A(I,J)*X(I)
   40     TEMP=TEMP+A(I,J)**2
          IF (TEMP .EQ. ZERO) THEN
              IF (IPRINT .GT. 0) WRITE(IPRINT,50)
   50         FORMAT (/4X,'Return from LINCOA because the gradient of',
     1          ' a constraint is zero.')
              GOTO 80
          END IF
          TEMP=DSQRT(TEMP)
          IF (SUM-B(J) .GT. SMALLX*TEMP) IFLAG=1
          W(IB+J-1)=DMAX1(B(J),SUM)/TEMP
          DO 60 I=1,N
          IW=IW+1
   60     W(IW)=A(I,J)/TEMP
      END IF
      IF (IFLAG .EQ. 1) THEN
          IF (IPRINT .GT. 0) WRITE(IPRINT,70)
   70     FORMAT (/4X,'LINCOA has made the initial X feasible by',
     1      ' increasing part(s) of B.')
      END IF
C
C     Partition the working space array, so that different parts of it can be
C     treated separately by the subroutine that performs the main calculation.
C
      NDIM=NPT+N
      IXB=IB+M
      IXP=IXB+N
      IFV=IXP+N*NPT
      IXS=IFV+NPT
      IXO=IXS+N
      IGO=IXO+N
      IHQ=IGO+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ISTP=IZMAT+NPT*NPTM
      ISP=ISTP+N
      IXN=ISP+NPT+NPT
      IAC=IXN+N
      IRC=IAC+N
      IQF=IRC+M
      IRF=IQF+N*N
      IPQW=IRF+(N*NP)/2
C
C     The above settings provide a partition of W for subroutine LINCOB.
C
      CALL LINCOB (N,NPT,M,W(IAMAT),W(IB),X,RHOBEG,RHOEND,IPRINT,
     1  MAXFUN,W(IXB),W(IXP),W(IFV),W(IXS),W(IXO),W(IGO),W(IHQ),
     2  W(IPQ),W(IBMAT),W(IZMAT),NDIM,W(ISTP),W(ISP),W(IXN),W(IAC),
     3  W(IRC),W(IQF),W(IRF),W(IPQW),W,CALFUN)
   80 RETURN
      END
      SUBROUTINE LINCOB (N,NPT,M,AMAT,B,X,RHOBEG,RHOEND,IPRINT,
     1  MAXFUN,XBASE,XPT,FVAL,XSAV,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,
     2  STEP,SP,XNEW,IACT,RESCON,QFAC,RFAC,PQW,W,CALFUN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AMAT(N,*),B(*),X(*),XBASE(*),XPT(NPT,*),FVAL(*),
     1  XSAV(*),XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),
     2  ZMAT(NPT,*),STEP(*),SP(*),XNEW(*),IACT(*),RESCON(*),
     3  QFAC(N,*),RFAC(*),PQW(*),W(*)
      EXTERNAL      CALFUN
C
C     The arguments N, NPT, M, X, RHOBEG, RHOEND, IPRINT and MAXFUN are
C       identical to the corresponding arguments in SUBROUTINE LINCOA.
C     AMAT is a matrix whose columns are the constraint gradients, scaled
C       so that they have unit length.
C     B contains on entry the right hand sides of the constraints, scaled
C       as above, but later B is modified for variables relative to XBASE.
C     XBASE holds a shift of origin that should reduce the contributions
C       from rounding errors to values of the model and Lagrange functions.
C     XPT contains the interpolation point coordinates relative to XBASE.
C     FVAL holds the values of F at the interpolation points.
C     XSAV holds the best feasible vector of variables so far, without any
C       shift of origin.
C     XOPT is set to XSAV-XBASE, which is the displacement from XBASE of
C       the feasible vector of variables that provides the least calculated
C       F so far, this vector being the current trust region centre.
C     GOPT holds the gradient of the quadratic model at XSAV = XBASE+XOPT.
C     HQ holds the explicit second derivatives of the quadratic model.
C     PQ contains the parameters of the implicit second derivatives of the
C       quadratic model.
C     BMAT holds the last N columns of the big inverse matrix H.
C     ZMAT holds the factorization of the leading NPT by NPT submatrix
C       of H, this factorization being ZMAT times Diag(DZ) times ZMAT^T,
C       where the elements of DZ are plus or minus one, as specified by IDZ.
C     NDIM is the first dimension of BMAT and has the value NPT+N.
C     STEP is employed for trial steps from XOPT. It is also used for working
C       space when XBASE is shifted and in PRELIM.
C     SP is reserved for the scalar products XOPT^T XPT(K,.), K=1,2,...,NPT,
C       followed by STEP^T XPT(K,.), K=1,2,...,NPT.
C     XNEW is the displacement from XBASE of the vector of variables for
C       the current calculation of F, except that SUBROUTINE TRSTEP uses it
C       for working space.
C     IACT is an integer array for the indices of the active constraints.
C     RESCON holds useful information about the constraint residuals. Every
C       nonnegative RESCON(J) is the residual of the J-th constraint at the
C       current trust region centre. Otherwise, if RESCON(J) is negative, the
C       J-th constraint holds as a strict inequality at the trust region
C       centre, its residual being at least |RESCON(J)|; further, the value
C       of |RESCON(J)| is at least the current trust region radius DELTA.
C     QFAC is the orthogonal part of the QR factorization of the matrix of
C       active constraint gradients, these gradients being ordered in
C       accordance with IACT. When NACT is less than N, columns are added
C       to QFAC to complete an N by N orthogonal matrix, which is important
C       for keeping calculated steps sufficiently close to the boundaries
C       of the active constraints.
C     RFAC is the upper triangular part of this QR factorization, beginning
C       with the first diagonal element, followed by the two elements in the
C       upper triangular part of the second column and so on.
C     PQW is used for working space, mainly for storing second derivative
C       coefficients of quadratic functions. Its length is NPT+N.
C     The array W is also used for working space. The required number of
C       elements, namely MAX[M+3*N,2*M+N,2*NPT], is set in LINCOA.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TENTH=0.1D0
      ZERO=0.0D0
      NP=N+1
      NH=(N*NP)/2
      NPTM=NPT-NP
C
C     Set the elements of XBASE, XPT, FVAL, XSAV, XOPT, GOPT, HQ, PQ, BMAT,
C       ZMAT and SP for the first iteration. An important feature is that,
C       if the interpolation point XPT(K,.) is not feasible, where K is any
C       integer from [1,NPT], then a change is made to XPT(K,.) if necessary
C       so that the constraint violation is at least 0.2*RHOBEG. Also KOPT
C       is set so that XPT(KOPT,.) is the initial trust region centre.
C
      CALL LINPRELIM (N,NPT,M,AMAT,B,X,RHOBEG,IPRINT,XBASE,XPT,FVAL,
     1  XSAV,XOPT,GOPT,KOPT,HQ,PQ,BMAT,ZMAT,IDZ,NDIM,SP,RESCON,
     2  STEP,PQW,W,CALFUN)
C
C     Begin the iterative procedure.
C
      NF=NPT
      FOPT=FVAL(KOPT)
      RHO=RHOBEG
      DELTA=RHO
      IFEAS=0
      NACT=0
      ITEST=3
   10 KNEW=0
      NVALA=0
      NVALB=0
C
C     Shift XBASE if XOPT may be too far from XBASE. First make the changes
C       to BMAT that do not depend on ZMAT.
C
   20 FSAVE=FOPT
      XOPTSQ=ZERO
      DO 30 I=1,N
   30 XOPTSQ=XOPTSQ+XOPT(I)**2
      IF (XOPTSQ .GE. 1.0D4*DELTA*DELTA) THEN
          QOPTSQ=0.25D0*XOPTSQ
          DO 50 K=1,NPT
          SUM=ZERO
          DO 40 I=1,N
   40     SUM=SUM+XPT(K,I)*XOPT(I)
          SUM=SUM-HALF*XOPTSQ
          W(NPT+K)=SUM
          SP(K)=ZERO
          DO 50 I=1,N
          XPT(K,I)=XPT(K,I)-HALF*XOPT(I)
          STEP(I)=BMAT(K,I)
          W(I)=SUM*XPT(K,I)+QOPTSQ*XOPT(I)
          IP=NPT+I
          DO 50 J=1,I
   50     BMAT(IP,J)=BMAT(IP,J)+STEP(I)*W(J)+W(I)*STEP(J)
C
C     Then the revisions of BMAT that depend on ZMAT are calculated.
C
          DO 90 K=1,NPTM
          SUMZ=ZERO
          DO 60 I=1,NPT
          SUMZ=SUMZ+ZMAT(I,K)
   60     W(I)=W(NPT+I)*ZMAT(I,K)
          DO 80 J=1,N
          SUM=QOPTSQ*SUMZ*XOPT(J)
          DO 70 I=1,NPT
   70     SUM=SUM+W(I)*XPT(I,J)
          STEP(J)=SUM
          IF (K .LT. IDZ) SUM=-SUM
          DO 80 I=1,NPT
   80     BMAT(I,J)=BMAT(I,J)+SUM*ZMAT(I,K)
          DO 90 I=1,N
          IP=I+NPT
          TEMP=STEP(I)
          IF (K .LT. IDZ) TEMP=-TEMP
          DO 90 J=1,I
   90     BMAT(IP,J)=BMAT(IP,J)+TEMP*STEP(J)
C
C     Update the right hand sides of the constraints.
C
          IF (M .GT. 0) THEN
              DO 110 J=1,M
              TEMP=ZERO
              DO 100 I=1,N
  100         TEMP=TEMP+AMAT(I,J)*XOPT(I)
  110         B(J)=B(J)-TEMP
          END IF
C
C     The following instructions complete the shift of XBASE, including the
C       changes to the parameters of the quadratic model.
C
          IH=0
          DO 130 J=1,N
          W(J)=ZERO
          DO 120 K=1,NPT
          W(J)=W(J)+PQ(K)*XPT(K,J)
  120     XPT(K,J)=XPT(K,J)-HALF*XOPT(J)
          DO 130 I=1,J
          IH=IH+1
          HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
  130     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 140 J=1,N
          XBASE(J)=XBASE(J)+XOPT(J)
          XOPT(J)=ZERO
  140     XPT(KOPT,J)=ZERO
      END IF
C
C     In the case KNEW=0, generate the next trust region step by calling
C       TRSTEP, where SNORM is the current trust region radius initially.
C       The final value of SNORM is the length of the calculated step,
C       except that SNORM is zero on return if the projected gradient is
C       unsuitable for starting the conjugate gradient iterations.
C
      DELSAV=DELTA
      KSAVE=KNEW
      IF (KNEW .EQ. 0) THEN
          SNORM=DELTA
          DO 150 I=1,N
  150     XNEW(I)=GOPT(I)
          CALL TRSTEP (N,NPT,M,AMAT,B,XPT,HQ,PQ,NACT,IACT,RESCON,
     1      QFAC,RFAC,SNORM,STEP,XNEW,W,W(M+1),PQW,PQW(NP),W(M+NP))
C
C     A trust region step is applied whenever its length, namely SNORM, is at
C       least HALF*DELTA. It is also applied if its length is at least 0.1999
C       times DELTA and if a line search of TRSTEP has caused a change to the
C       active set. Otherwise there is a branch below to label 530 or 560.
C
          TEMP=HALF*DELTA
          IF (XNEW(1) .GE. HALF) TEMP=0.1999D0*DELTA
          IF (SNORM .LE. TEMP) THEN
              DELTA=HALF*DELTA
              IF (DELTA .LE. 1.4D0*RHO) DELTA=RHO
              NVALA=NVALA+1
              NVALB=NVALB+1
              TEMP=SNORM/RHO
              IF (DELSAV .GT. RHO) TEMP=ONE
              IF (TEMP .GE. HALF) NVALA=ZERO
              IF (TEMP .GE. TENTH) NVALB=ZERO
              IF (DELSAV .GT. RHO) GOTO 530
              IF (NVALA .LT. 5 .AND. NVALB .LT. 3) GOTO 530
              IF (SNORM .GT. ZERO) KSAVE=-1
              GOTO 560
          END IF
          NVALA=ZERO
          NVALB=ZERO
C
C     Alternatively, KNEW is positive. Then the model step is calculated
C       within a trust region of radius DEL, after setting the gradient at
C       XBASE and the second derivative parameters of the KNEW-th Lagrange
C       function in W(1) to W(N) and in PQW(1) to PQW(NPT), respectively.
C
      ELSE
          DEL=DMAX1(TENTH*DELTA,RHO)
          DO 160 I=1,N
  160     W(I)=BMAT(KNEW,I)
          DO 170 K=1,NPT
  170     PQW(K)=ZERO
          DO 180 J=1,NPTM
          TEMP=ZMAT(KNEW,J)
          IF (J .LT. IDZ) TEMP=-TEMP
          DO 180 K=1,NPT
  180     PQW(K)=PQW(K)+TEMP*ZMAT(K,J)
          CALL QMSTEP (N,NPT,M,AMAT,B,XPT,XOPT,NACT,IACT,RESCON,
     1      QFAC,KOPT,KNEW,DEL,STEP,W,PQW,W(NP),W(NP+M),IFEAS)
      END IF
C
C     Set VQUAD to the change to the quadratic model when the move STEP is
C       made from XOPT. If STEP is a trust region step, then VQUAD should be
C       negative. If it is nonnegative due to rounding errors in this case,
C       there is a branch to label 530 to try to improve the model.
C
      VQUAD=ZERO
      IH=0
      DO 190 J=1,N
      VQUAD=VQUAD+STEP(J)*GOPT(J)
      DO 190 I=1,J
      IH=IH+1
      TEMP=STEP(I)*STEP(J)
      IF (I .EQ. J) TEMP=HALF*TEMP
  190 VQUAD=VQUAD+TEMP*HQ(IH)
      DO 210 K=1,NPT
      TEMP=ZERO
      DO 200 J=1,N
      TEMP=TEMP+XPT(K,J)*STEP(J)
  200 SP(NPT+K)=TEMP
  210 VQUAD=VQUAD+HALF*PQ(K)*TEMP*TEMP
      IF (KSAVE .EQ. 0 .AND. VQUAD .GE. ZERO) GOTO 530
C
C     Calculate the next value of the objective function. The difference
C       between the actual new value of F and the value predicted by the
C       model is recorded in DIFF.
C
  220 NF=NF+1
      IF (NF .GT. MAXFUN) THEN
          NF=NF-1
          IF (IPRINT .GT. 0) WRITE(IPRINT,230)
  230     FORMAT (/4X,'Return from LINCOA because CALFUN has been',
     1      ' called MAXFUN times.')
          GOTO 600
      END IF
      XDIFF=ZERO
      DO 240 I=1,N
      XNEW(I)=XOPT(I)+STEP(I)
      X(I)=XBASE(I)+XNEW(I)
  240 XDIFF=XDIFF+(X(I)-XSAV(I))**2
      XDIFF=DSQRT(XDIFF)
      IF (KSAVE .EQ. -1) XDIFF=RHO
      IF (XDIFF .LE. TENTH*RHO .OR. XDIFF .GE. DELTA+DELTA) THEN
          IFEAS=0
          IF (IPRINT .GT. 0) WRITE(IPRINT,250)
  250     FORMAT (/4X,'Return from LINCOA because rounding errors',
     1      ' prevent reasonable changes to X.')
          GOTO 600
      END IF
      IF (KSAVE .LE. 0) IFEAS=1
      F=DFLOAT(IFEAS)
      CALL CALFUN (N,X,F)
      IF (IPRINT .GE. 3) THEN
          WRITE(IPRINT,260) NF,F,(X(I),I=1,N)
  260     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1      '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (KSAVE .EQ. -1) GOTO 600
      DIFF=F-FOPT-VQUAD
C
C     If X is feasible, then set DFFALT to the difference between the new
C       value of F and the value predicted by the alternative model.
C
      IF (IFEAS .EQ. 1 .AND. ITEST .LT. 3) THEN
          DO 270 K=1,NPT
          PQW(K)=ZERO
  270     W(K)=FVAL(K)-FVAL(KOPT)
          DO 290 J=1,NPTM
          SUM=ZERO
          DO 280 I=1,NPT
  280     SUM=SUM+W(I)*ZMAT(I,J)
          IF (J .LT. IDZ) SUM=-SUM
          DO 290 K=1,NPT
  290     PQW(K)=PQW(K)+SUM*ZMAT(K,J)
          VQALT=ZERO
          DO 310 K=1,NPT
          SUM=ZERO
          DO 300 J=1,N
  300     SUM=SUM+BMAT(K,J)*STEP(J)
          VQALT=VQALT+SUM*W(K)
  310     VQALT=VQALT+PQW(K)*SP(NPT+K)*(HALF*SP(NPT+K)+SP(K))
          DFFALT=F-FOPT-VQALT
      END IF
      IF (ITEST .EQ. 3) THEN
          DFFALT=DIFF
          ITEST=0
      END IF
C
C     Pick the next value of DELTA after a trust region step.
C
      IF (KSAVE .EQ. 0) THEN
          RATIO=(F-FOPT)/VQUAD
          IF (RATIO .LE. TENTH) THEN
              DELTA=HALF*DELTA
          ELSE IF (RATIO .LE. 0.7D0) THEN
              DELTA=DMAX1(HALF*DELTA,SNORM)
          ELSE 
              TEMP=DSQRT(2.0D0)*DELTA
              DELTA=DMAX1(HALF*DELTA,SNORM+SNORM)
              DELTA=DMIN1(DELTA,TEMP)
          END IF
          IF (DELTA .LE. 1.4D0*RHO) DELTA=RHO
      END IF
C
C     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
C       can be moved. If STEP is a trust region step, then KNEW is zero at
C       present, but a positive value is picked by subroutine UPDATE.
C
      CALL LINUPDATE (N,NPT,XPT,BMAT,ZMAT,IDZ,NDIM,SP,STEP,KOPT,
     1  KNEW,PQW,W)
      IF (KNEW .EQ. 0) THEN
          IF (IPRINT .GT. 0) WRITE(IPRINT,320)
  320     FORMAT (/4X,'Return from LINCOA because the denominator'
     1      ' of the updating formula is zero.')
          GOTO 600
      END IF
C
C     If ITEST is increased to 3, then the next quadratic model is the
C       one whose second derivative matrix is least subject to the new
C       interpolation conditions. Otherwise the new model is constructed
C       by the symmetric Broyden method in the usual way.
C
      IF (IFEAS .EQ. 1) THEN
          ITEST=ITEST+1
          IF (DABS(DFFALT) .GE. TENTH*DABS(DIFF)) ITEST=0
      END IF
C
C     Update the second derivatives of the model by the symmetric Broyden
C       method, using PQW for the second derivative parameters of the new
C       KNEW-th Lagrange function. The contribution from the old parameter
C       PQ(KNEW) is included in the second derivative matrix HQ. W is used
C       later for the gradient of the new KNEW-th Lagrange function.       
C
      IF (ITEST .LT. 3) THEN
          DO 330 K=1,NPT
  330     PQW(K)=ZERO
          DO 350 J=1,NPTM
          TEMP=ZMAT(KNEW,J)
          IF (TEMP .NE. ZERO) THEN
              IF (J .LT. IDZ) TEMP=-TEMP
              DO 340 K=1,NPT
  340         PQW(K)=PQW(K)+TEMP*ZMAT(K,J)
          END IF
  350     CONTINUE
          IH=0
          DO 360 I=1,N
          W(I)=BMAT(KNEW,I)
          TEMP=PQ(KNEW)*XPT(KNEW,I)
          DO 360 J=1,I
          IH=IH+1
  360     HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
          PQ(KNEW)=ZERO
          DO 370 K=1,NPT
  370     PQ(K)=PQ(K)+DIFF*PQW(K)
      END IF
C
C     Include the new interpolation point with the corresponding updates of
C       SP. Also make the changes of the symmetric Broyden method to GOPT at
C       the old XOPT if ITEST is less than 3.
C
      FVAL(KNEW)=F
      SP(KNEW)=SP(KOPT)+SP(NPT+KOPT)
      SSQ=ZERO
      DO 380 I=1,N
      XPT(KNEW,I)=XNEW(I)
  380 SSQ=SSQ+STEP(I)**2
      SP(NPT+KNEW)=SP(NPT+KOPT)+SSQ
      IF (ITEST .LT. 3) THEN
          DO 390 K=1,NPT
          TEMP=PQW(K)*SP(K)
          DO 390 I=1,N
  390     W(I)=W(I)+TEMP*XPT(K,I)
          DO 400 I=1,N
  400     GOPT(I)=GOPT(I)+DIFF*W(I)
      END IF
C
C     Update FOPT, XSAV, XOPT, KOPT, RESCON and SP if the new F is the
C       least calculated value so far with a feasible vector of variables.
C
      IF (F .LT. FOPT .AND. IFEAS .EQ. 1) THEN
          FOPT=F
          DO 410 J=1,N
          XSAV(J)=X(J)
  410     XOPT(J)=XNEW(J)
          KOPT=KNEW
          SNORM=DSQRT(SSQ)
          DO 430 J=1,M
          IF (RESCON(J) .GE. DELTA+SNORM) THEN
              RESCON(J)=SNORM-RESCON(J)
          ELSE
              RESCON(J)=RESCON(J)+SNORM
              IF (RESCON(J)+DELTA .GT. ZERO) THEN
                  TEMP=B(J)
                  DO 420 I=1,N
  420             TEMP=TEMP-XOPT(I)*AMAT(I,J)
                  TEMP=DMAX1(TEMP,ZERO)
                  IF (TEMP .GE. DELTA) TEMP=-TEMP
                  RESCON(J)=TEMP
              END IF
          END IF
  430     CONTINUE
          DO 440 K=1,NPT
  440     SP(K)=SP(K)+SP(NPT+K)
C
C     Also revise GOPT when symmetric Broyden updating is applied.
C
          IF (ITEST .LT. 3) THEN
              IH=0
              DO 450 J=1,N
              DO 450 I=1,J
              IH=IH+1
              IF (I .LT. J) GOPT(J)=GOPT(J)+HQ(IH)*STEP(I)
  450         GOPT(I)=GOPT(I)+HQ(IH)*STEP(J)
              DO 460 K=1,NPT
              TEMP=PQ(K)*SP(NPT+K)
              DO 460 I=1,N
  460         GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
          END IF
      END IF
C
C     Replace the current model by the least Frobenius norm interpolant if
C       this interpolant gives substantial reductions in the predictions
C       of values of F at feasible points.
C
      IF (ITEST .EQ. 3) THEN
          DO 470 K=1,NPT
          PQ(K)=ZERO
  470     W(K)=FVAL(K)-FVAL(KOPT)
          DO 490 J=1,NPTM
          SUM=ZERO
          DO 480 I=1,NPT
  480     SUM=SUM+W(I)*ZMAT(I,J)
          IF (J .LT. IDZ) SUM=-SUM
          DO 490 K=1,NPT
  490     PQ(K)=PQ(K)+SUM*ZMAT(K,J)
          DO 500 J=1,N
          GOPT(J)=ZERO
          DO 500 I=1,NPT
  500     GOPT(J)=GOPT(J)+W(I)*BMAT(I,J)
          DO 510 K=1,NPT
          TEMP=PQ(K)*SP(K)
          DO 510 I=1,N
  510     GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
          DO 520 IH=1,NH
  520     HQ(IH)=ZERO
      END IF
C
C     If a trust region step has provided a sufficient decrease in F, then
C       branch for another trust region calculation. Every iteration that
C       takes a model step is followed by an attempt to take a trust region
C       step.
C
      KNEW=0
      IF (KSAVE .GT. 0) GOTO 20
      IF (RATIO .GE. TENTH) GOTO 20
C
C     Alternatively, find out if the interpolation points are close enough
C       to the best point so far.
C
  530 DISTSQ=DMAX1(DELTA*DELTA,4.0D0*RHO*RHO)
      DO 550 K=1,NPT
      SUM=ZERO
      DO 540 J=1,N
  540 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM .GT. DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  550 CONTINUE
C
C     If KNEW is positive, then branch back for the next iteration, which
C       will generate a "model step". Otherwise, if the current iteration
C       has reduced F, or if DELTA was above its lower bound when the last
C       trust region step was calculated, then try a "trust region" step
C       instead.
C
      IF (KNEW .GT. 0) GOTO 20
      KNEW=0
      IF (FOPT .LT. FSAVE) GOTO 20
      IF (DELSAV .GT. RHO) GOTO 20
C
C     The calculations with the current value of RHO are complete.
C       Pick the next value of RHO.
C
  560 IF (RHO .GT. RHOEND) THEN
          DELTA=HALF*RHO
          IF (RHO .GT. 250.0D0*RHOEND) THEN
              RHO=TENTH*RHO
          ELSE IF (RHO .LE. 16.0D0*RHOEND) THEN
              RHO=RHOEND
          ELSE
              RHO=DSQRT(RHO*RHOEND)
          END IF 
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) WRITE(IPRINT,570)
  570         FORMAT (5X)
              WRITE(IPRINT,580) RHO,NF
  580         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              WRITE(IPRINT,590) FOPT,(XBASE(I)+XOPT(I),I=1,N)
  590         FORMAT (4X,'Least value of F =',1PD23.15,9X,
     1          'The corresponding X is:'/(2X,5D15.6))
          END IF
          GOTO 10
      END IF
C
C     Return from the calculation, after branching to label 220 for another
C       Newton-Raphson step if it has not been tried before.
C
      IF (KSAVE .EQ. -1) GOTO 220
  600 IF (FOPT .LE. F .OR. IFEAS .EQ. 0) THEN
          DO 610 I=1,N
  610     X(I)=XSAV(I)
          F=FOPT
      END IF
      IF (IPRINT .GE. 1) THEN
          WRITE(IPRINT,620) NF
  620     FORMAT (/4X,'At the return from LINCOA',5X,
     1      'Number of function values =',I6)
          WRITE(IPRINT,590) F,(X(I),I=1,N)
      END IF
      W(1)=F
      W(2)=DFLOAT(NF)+HALF
      RETURN
      END
      SUBROUTINE TRSTEP (N,NPT,M,AMAT,B,XPT,HQ,PQ,NACT,IACT,RESCON,
     1  QFAC,RFAC,SNORM,STEP,G,RESNEW,RESACT,D,DW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AMAT(N,*),B(*),XPT(NPT,*),HQ(*),PQ(*),IACT(*),
     1  RESCON(*),QFAC(N,*),RFAC(*),STEP(*),G(*),RESNEW(*),RESACT(*),
     2  D(*),DW(*),W(*)
C
C     N, NPT, M, AMAT, B, XPT, HQ, PQ, NACT, IACT, RESCON, QFAC and RFAC
C       are the same as the terms with these names in LINCOB. If RESCON(J)
C       is negative, then |RESCON(J)| must be no less than the trust region
C       radius, so that the J-th constraint can be ignored.
C     SNORM is set to the trust region radius DELTA initially. On the
C       return, however, it is the length of the calculated STEP, which is
C       set to zero if the constraints do not allow a long enough step.
C     STEP is the total calculated step so far from the trust region centre,
C       its final value being given by the sequence of CG iterations, which
C       terminate if the trust region boundary is reached.
C     G must be set on entry to the gradient of the quadratic model at the
C       trust region centre. It is used as working space, however, and is
C       always the gradient of the model at the current STEP, except that
C       on return the value of G(1) is set to ONE instead of to ZERO if
C       and only if GETACT is called more than once.
C     RESNEW, RESACT, D, DW and W are used for working space. A negative
C       value of RESNEW(J) indicates that the J-th constraint does not
C       restrict the CG steps of the current trust region calculation, a
C       zero value of RESNEW(J) indicates that the J-th constraint is active,
C       and otherwise RESNEW(J) is set to the greater of TINY and the actual
C       residual of the J-th constraint for the current STEP. RESACT holds
C       the residuals of the active constraints, which may be positive.
C       D is the search direction of each line search. DW is either another
C       search direction or the change in gradient along D. The length of W
C       must be at least MAX[M,2*N].
C
C     Set some numbers for the conjugate gradient iterations.
C
      HALF=0.5D0
      ONE=1.0D0
      TINY=1.0D-60
      ZERO=0.0D0
      CTEST=0.01D0
      SNSQ=SNORM*SNORM
C
C     Set the initial elements of RESNEW, RESACT and STEP.
C
      IF (M .GT. 0) THEN
          DO 10 J=1,M
          RESNEW(J)=RESCON(J)
          IF (RESCON(J) .GE. SNORM) THEN
              RESNEW(J)=-ONE
          ELSE IF (RESCON(J) .GE. ZERO) THEN
              RESNEW(J)=DMAX1(RESNEW(J),TINY)
          END IF
   10     CONTINUE
          IF (NACT .GT. 0) THEN
              DO 20 K=1,NACT
              RESACT(K)=RESCON(IACT(K))
   20         RESNEW(IACT(K))=ZERO
          END IF
      END IF
      DO 30 I=1,N
   30 STEP(I)=ZERO
      SS=ZERO
      REDUCT=ZERO
      NCALL=0
C
C     GETACT picks the active set for the current STEP. It also sets DW to
C       the vector closest to -G that is orthogonal to the normals of the
C       active constraints. DW is scaled to have length 0.2*SNORM, as then
C       a move of DW from STEP is allowed by the linear constraints.
C
   40 NCALL=NCALL+1
      CALL GETACT (N,M,AMAT,B,NACT,IACT,QFAC,RFAC,SNORM,RESNEW,
     1  RESACT,G,DW,W,W(N+1))
      IF (W(N+1) .EQ. ZERO) GOTO 320
      SCALE=0.2D0*SNORM/DSQRT(W(N+1))
      DO 50 I=1,N
   50 DW(I)=SCALE*DW(I)
C
C     If the modulus of the residual of an active constraint is substantial,
C       then set D to the shortest move from STEP to the boundaries of the
C       active constraints.
C
      RESMAX=ZERO
      IF (NACT .GT. 0) THEN
          DO 60 K=1,NACT
   60     RESMAX=DMAX1(RESMAX,RESACT(K))
      END IF
      GAMMA=ZERO
      IF (RESMAX .GT. 1.0D-4*SNORM) THEN
          IR=0
          DO 80 K=1,NACT
          TEMP=RESACT(K)
          IF (K .GE. 2) THEN
              DO 70 I=1,K-1
              IR=IR+1
   70         TEMP=TEMP-RFAC(IR)*W(I)
          END IF
          IR=IR+1
   80     W(K)=TEMP/RFAC(IR)
          DO 90 I=1,N
          D(I)=ZERO
          DO 90 K=1,NACT
   90     D(I)=D(I)+W(K)*QFAC(I,K)
C
C     The vector D that has just been calculated is also the shortest move
C       from STEP+DW to the boundaries of the active constraints. Set GAMMA
C       to the greatest steplength of this move that satisfies the trust
C       region bound.
C
          RHS=SNSQ
          DS=ZERO
          DD=ZERO
          DO 100 I=1,N
          SUM=STEP(I)+DW(I)
          RHS=RHS-SUM*SUM
          DS=DS+D(I)*SUM
  100     DD=DD+D(I)**2
          IF (RHS .GT. ZERO) THEN
              TEMP=DSQRT(DS*DS+DD*RHS)
              IF (DS .LE. ZERO) THEN
                  GAMMA=(TEMP-DS)/DD
              ELSE
                  GAMMA=RHS/(TEMP+DS)
              END IF
          END IF
C
C     Reduce the steplength GAMMA if necessary so that the move along D
C       also satisfies the linear constraints.
C
          J=0
  110     IF (GAMMA .GT. ZERO) THEN
              J=J+1
              IF (RESNEW(J) .GT. ZERO) THEN
                  AD=ZERO
                  ADW=ZERO
                  DO 120 I=1,N
                  AD=AD+AMAT(I,J)*D(I)
  120             ADW=ADW+AMAT(I,J)*DW(I)
                  IF (AD .GT. ZERO) THEN
                      TEMP=DMAX1((RESNEW(J)-ADW)/AD,ZERO)
                      GAMMA=DMIN1(GAMMA,TEMP)
                  END IF
              END IF
              IF (J .LT. M) GOTO 110
          END IF
          GAMMA=DMIN1(GAMMA,ONE)
      END IF
C
C     Set the next direction for seeking a reduction in the model function
C       subject to the trust region bound and the linear constraints.
C
      IF (GAMMA .LE. ZERO) THEN
          DO 130 I=1,N
  130     D(I)=DW(I)
          ICOUNT=NACT
      ELSE
          DO 140 I=1,N
  140     D(I)=DW(I)+GAMMA*D(I)
          ICOUNT=NACT-1
      END IF
      ALPBD=ONE
C
C     Set ALPHA to the steplength from STEP along D to the trust region
C       boundary. Return if the first derivative term of this step is
C       sufficiently small or if no further progress is possible.
C
  150 ICOUNT=ICOUNT+1
      RHS=SNSQ-SS
      IF (RHS .LE. ZERO) GOTO 320
      DG=ZERO
      DS=ZERO
      DD=ZERO
      DO 160 I=1,N
      DG=DG+D(I)*G(I)
      DS=DS+D(I)*STEP(I)
  160 DD=DD+D(I)**2
      IF (DG .GE. ZERO) GOTO 320
      TEMP=DSQRT(RHS*DD+DS*DS)
      IF (DS .LE. ZERO) THEN
          ALPHA=(TEMP-DS)/DD
      ELSE
          ALPHA=RHS/(TEMP+DS)
      END IF
      IF (-ALPHA*DG .LE. CTEST*REDUCT) GOTO 320
C
C     Set DW to the change in gradient along D.
C
      IH=0
      DO 170 J=1,N
      DW(J)=ZERO
      DO 170 I=1,J
      IH=IH+1
      IF (I .LT. J) DW(J)=DW(J)+HQ(IH)*D(I)
  170 DW(I)=DW(I)+HQ(IH)*D(J)
      DO 190 K=1,NPT
      TEMP=ZERO
      DO 180 J=1,N
  180 TEMP=TEMP+XPT(K,J)*D(J)
      TEMP=PQ(K)*TEMP
      DO 190 I=1,N
  190 DW(I)=DW(I)+TEMP*XPT(K,I)
C
C     Set DGD to the curvature of the model along D. Then reduce ALPHA if
C       necessary to the value that minimizes the model.
C
      DGD=ZERO
      DO 200 I=1,N
  200 DGD=DGD+D(I)*DW(I)
      ALPHT=ALPHA
      IF (DG+ALPHA*DGD .GT. ZERO) THEN
          ALPHA=-DG/DGD
      END IF
C
C     Make a further reduction in ALPHA if necessary to preserve feasibility,
C       and put some scalar products of D with constraint gradients in W.
C
      ALPHM=ALPHA
      JSAV=0
      IF (M .GT. 0) THEN
          DO 220 J=1,M
          AD=ZERO
          IF (RESNEW(J) .GT. ZERO) THEN
              DO 210 I=1,N
  210         AD=AD+AMAT(I,J)*D(I)
              IF (ALPHA*AD .GT. RESNEW(J)) THEN
                  ALPHA=RESNEW(J)/AD
                  JSAV=J
              END IF
          END IF
  220     W(J)=AD
      END IF
      ALPHA=DMAX1(ALPHA,ALPBD)
      ALPHA=DMIN1(ALPHA,ALPHM)
      IF (ICOUNT .EQ. NACT) ALPHA=DMIN1(ALPHA,ONE)
C
C     Update STEP, G, RESNEW, RESACT and REDUCT.
C
      SS=ZERO
      DO 230 I=1,N
      STEP(I)=STEP(I)+ALPHA*D(I)
      SS=SS+STEP(I)**2
  230 G(I)=G(I)+ALPHA*DW(I)
      IF (M .GT. 0) THEN
          DO 240 J=1,M
          IF (RESNEW(J) .GT. ZERO) THEN
              RESNEW(J)=DMAX1(RESNEW(J)-ALPHA*W(J),TINY)
          END IF
  240     CONTINUE
      END IF
      IF (ICOUNT .EQ. NACT .AND. NACT .GT. 0) THEN
          DO 250 K=1,NACT
  250     RESACT(K)=(ONE-GAMMA)*RESACT(K)
      END IF
      REDUCT=REDUCT-ALPHA*(DG+HALF*ALPHA*DGD)
C
C     Test for termination. Branch to label 40 if there is a new active
C       constraint and if the distance from STEP to the trust region
C       boundary is at least 0.2*SNORM.
C
      IF (ALPHA .EQ. ALPHT) GOTO 320
      TEMP=-ALPHM*(DG+HALF*ALPHM*DGD)
      IF (TEMP .LE. CTEST*REDUCT) GOTO 320
      IF (JSAV .GT. 0) THEN
          IF (SS .LE. 0.64D0*SNSQ) GOTO 40
          GOTO 320
      END IF
      IF (ICOUNT .EQ. N) GOTO 320
C
C     Calculate the next search direction, which is conjugate to the
C       previous one except in the case ICOUNT=NACT.
C
      IF (NACT .GT. 0) THEN
          DO 260 J=NACT+1,N
          W(J)=ZERO
          DO 260 I=1,N
  260     W(J)=W(J)+G(I)*QFAC(I,J)
          DO 280 I=1,N
          TEMP=ZERO
          DO 270 J=NACT+1,N
  270     TEMP=TEMP+QFAC(I,J)*W(J)
  280     W(N+I)=TEMP
      ELSE
          DO 290 I=1,N
  290     W(N+I)=G(I)
      END IF
      IF (ICOUNT .EQ. NACT) THEN
          BETA=ZERO
      ELSE
          WGD=ZERO
          DO 300 I=1,N
  300     WGD=WGD+W(N+I)*DW(I)
          BETA=WGD/DGD
      END IF
      DO 310 I=1,N
  310 D(I)=-W(N+I)+BETA*D(I)
      ALPBD=ZERO
      GOTO 150
C
C     Return from the subroutine.
C
  320 SNORM=ZERO
      IF (REDUCT .GT. ZERO) SNORM=DSQRT(SS)
      G(1)=ZERO
      IF (NCALL .GT. 1) G(1)=ONE
      RETURN
      END
      SUBROUTINE GETACT (N,M,AMAT,B,NACT,IACT,QFAC,RFAC,SNORM,
     1  RESNEW,RESACT,G,DW,VLAM,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AMAT(N,*),B(*),IACT(*),QFAC(N,*),RFAC(*),
     1  RESNEW(*),RESACT(*),G(*),DW(*),VLAM(*),W(*)
C
C     N, M, AMAT, B, NACT, IACT, QFAC and RFAC are the same as the terms
C       with these names in SUBROUTINE LINCOB. The current values must be
C       set on entry. NACT, IACT, QFAC and RFAC are kept up to date when
C       GETACT changes the current active set.
C     SNORM, RESNEW, RESACT, G and DW are the same as the terms with these
C       names in SUBROUTINE TRSTEP. The elements of RESNEW and RESACT are
C       also kept up to date.
C     VLAM and W are used for working space, the vector VLAM being reserved
C       for the Lagrange multipliers of the calculation. Their lengths must
C       be at least N.
C     The main purpose of GETACT is to pick the current active set. It is
C       defined by the property that the projection of -G into the space
C       orthogonal to the active constraint normals is as large as possible,
C       subject to this projected steepest descent direction moving no closer
C       to the boundary of every constraint whose current residual is at most
C       0.2*SNORM. On return, the settings in NACT, IACT, QFAC and RFAC are
C       all appropriate to this choice of active set.
C     Occasionally this projected direction is zero, and then the final value
C       of W(1) is set to zero. Otherwise, the direction itself is returned
C       in DW, and W(1) is set to the square of the length of the direction.
C
C     Set some constants and a temporary VLAM.
C
      ONE=1.0D0
      TINY=1.0D-60
      ZERO=0.0D0
      TDEL=0.2D0*SNORM
      DDSAV=ZERO
      DO 10 I=1,N
      DDSAV=DDSAV+G(I)**2
   10 VLAM(I)=ZERO
      DDSAV=DDSAV+DDSAV
C
C     Set the initial QFAC to the identity matrix in the case NACT=0.
C
      IF (NACT .EQ. 0) THEN
          DO 30 I=1,N
          DO 20 J=1,N
   20     QFAC(I,J)=ZERO
   30     QFAC(I,I)=ONE
          GOTO 100
      END IF
C
C     Remove any constraints from the initial active set whose residuals
C       exceed TDEL.
C
      IFLAG=1
      IC=NACT
   40 IF (RESACT(IC) .GT. TDEL) GOTO 800
   50 IC=IC-1
      IF (IC .GT. 0) GOTO 40
C
C     Remove any constraints from the initial active set whose Lagrange
C       multipliers are nonnegative, and set the surviving multipliers.
C
      IFLAG=2
   60 IF (NACT .EQ. 0) GOTO 100
      IC=NACT
   70 TEMP=ZERO
      DO 80 I=1,N
   80 TEMP=TEMP+QFAC(I,IC)*G(I)
      IDIAG=(IC*IC+IC)/2
      IF (IC .LT. NACT) THEN
          JW=IDIAG+IC
          DO 90 J=IC+1,NACT
          TEMP=TEMP-RFAC(JW)*VLAM(J)
   90     JW=JW+J
      END IF
      IF (TEMP .GE. ZERO) GOTO 800
      VLAM(IC)=TEMP/RFAC(IDIAG)
      IC=IC-1
      IF (IC .GT. 0) GOTO 70
C
C     Set the new search direction D. Terminate if the 2-norm of D is zero
C       or does not decrease, or if NACT=N holds. The situation NACT=N
C       occurs for sufficiently large SNORM if the origin is in the convex
C       hull of the constraint gradients.
C
  100 IF (NACT .EQ. N) GOTO 290
      DO 110 J=NACT+1,N
      W(J)=ZERO
      DO 110 I=1,N
  110 W(J)=W(J)+QFAC(I,J)*G(I)
      DD=ZERO
      DO 130 I=1,N
      DW(I)=ZERO
      DO 120 J=NACT+1,N
  120 DW(I)=DW(I)-W(J)*QFAC(I,J)
  130 DD=DD+DW(I)**2
      IF (DD .GE. DDSAV) GOTO 290
      IF (DD .EQ. ZERO) GOTO 300
      DDSAV=DD
      DNORM=DSQRT(DD)
C
C     Pick the next integer L or terminate, a positive value of L being
C       the index of the most violated constraint. The purpose of CTOL
C       below is to estimate whether a positive value of VIOLMX may be
C       due to computer rounding errors.
C
      L=0
      IF (M .GT. 0) THEN
          TEST=DNORM/SNORM
          VIOLMX=ZERO
          DO 150 J=1,M
          IF (RESNEW(J) .GT. ZERO .AND. RESNEW(J) .LE. TDEL) THEN
              SUM=ZERO
              DO 140 I=1,N
  140         SUM=SUM+AMAT(I,J)*DW(I)
              IF (SUM .GT. TEST*RESNEW(J)) THEN
                  IF (SUM .GT. VIOLMX) THEN
                      L=J
                      VIOLMX=SUM
                  END IF
              END IF
          END IF
  150     CONTINUE
          CTOL=ZERO
          TEMP=0.01D0*DNORM
          IF (VIOLMX .GT. ZERO .AND. VIOLMX .LT. TEMP) THEN
              IF (NACT .GT. 0) THEN
                  DO 170 K=1,NACT
                  J=IACT(K)
                  SUM=ZERO
                  DO 160 I=1,N
  160             SUM=SUM+DW(I)*AMAT(I,J)
  170             CTOL=DMAX1(CTOL,DABS(SUM))
              END IF
          END IF
      END IF
      W(1)=ONE
      IF (L .EQ. 0) GOTO 300
      IF (VIOLMX .LE. 10.0D0*CTOL) GOTO 300
C
C     Apply Givens rotations to the last (N-NACT) columns of QFAC so that
C       the first (NACT+1) columns of QFAC are the ones required for the
C       addition of the L-th constraint, and add the appropriate column
C       to RFAC.
C
      NACTP=NACT+1
      IDIAG=(NACTP*NACTP-NACTP)/2
      RDIAG=ZERO
      DO 200 J=N,1,-1
      SPROD=ZERO
      DO 180 I=1,N
  180 SPROD=SPROD+QFAC(I,J)*AMAT(I,L)
      IF (J .LE. NACT) THEN
          RFAC(IDIAG+J)=SPROD
      ELSE
          IF (DABS(RDIAG) .LE. 1.0D-20*DABS(SPROD)) THEN
              RDIAG=SPROD
          ELSE
              TEMP=DSQRT(SPROD*SPROD+RDIAG*RDIAG)
              COSV=SPROD/TEMP
              SINV=RDIAG/TEMP
              RDIAG=TEMP
              DO 190 I=1,N
              TEMP=COSV*QFAC(I,J)+SINV*QFAC(I,J+1)
              QFAC(I,J+1)=-SINV*QFAC(I,J)+COSV*QFAC(I,J+1)
  190         QFAC(I,J)=TEMP
          END IF
      END IF
  200 CONTINUE
      IF (RDIAG .LT. ZERO) THEN
          DO 210 I=1,N
  210     QFAC(I,NACTP)=-QFAC(I,NACTP)
      END IF
      RFAC(IDIAG+NACTP)=DABS(RDIAG)
      NACT=NACTP
      IACT(NACT)=L
      RESACT(NACT)=RESNEW(L)
      VLAM(NACT)=ZERO
      RESNEW(L)=ZERO
C
C     Set the components of the vector VMU in W.
C
  220 W(NACT)=ONE/RFAC((NACT*NACT+NACT)/2)**2
      IF (NACT .GT. 1) THEN
          DO 240 I=NACT-1,1,-1
          IDIAG=(I*I+I)/2
          JW=IDIAG+I
          SUM=ZERO
          DO 230 J=I+1,NACT
          SUM=SUM-RFAC(JW)*W(J)
  230     JW=JW+J
  240     W(I)=SUM/RFAC(IDIAG)
      END IF
C
C     Calculate the multiple of VMU to subtract from VLAM, and update VLAM.
C
      VMULT=VIOLMX
      IC=0
      J=1
  250 IF (J .LT. NACT) THEN
          IF (VLAM(J) .GE. VMULT*W(J)) THEN
              IC=J
              VMULT=VLAM(J)/W(J)
          END IF
          J=J+1
          GOTO 250
      END IF
      DO 260 J=1,NACT
  260 VLAM(J)=VLAM(J)-VMULT*W(J)
      IF (IC .GT. 0) VLAM(IC)=ZERO
      VIOLMX=DMAX1(VIOLMX-VMULT,ZERO)
      IF (IC .EQ. 0) VIOLMX=ZERO
C
C     Reduce the active set if necessary, so that all components of the
C       new VLAM are negative, with resetting of the residuals of the
C       constraints that become inactive.
C
      IFLAG=3
      IC=NACT
  270 IF (VLAM(IC) .LT. ZERO) GOTO 280
      RESNEW(IACT(IC))=DMAX1(RESACT(IC),TINY)
      GOTO 800
  280 IC=IC-1
      IF (IC .GT. 0) GOTO 270
C
C     Calculate the next VMU if VIOLMX is positive. Return if NACT=N holds,
C       as then the active constraints imply D=0. Otherwise, go to label
C       100, to calculate the new D and to test for termination.
C
      IF (VIOLMX .GT. ZERO) GOTO 220
      IF (NACT .LT. N) GOTO 100
  290 DD=ZERO
  300 W(1)=DD
      RETURN
C
C     These instructions rearrange the active constraints so that the new
C       value of IACT(NACT) is the old value of IACT(IC). A sequence of
C       Givens rotations is applied to the current QFAC and RFAC. Then NACT
C       is reduced by one.
C
  800 RESNEW(IACT(IC))=DMAX1(RESACT(IC),TINY)
      JC=IC
  810 IF (JC .LT. NACT) THEN
          JCP=JC+1
          IDIAG=JC*JCP/2
          JW=IDIAG+JCP
          TEMP=DSQRT(RFAC(JW-1)**2+RFAC(JW)**2)
          CVAL=RFAC(JW)/TEMP
          SVAL=RFAC(JW-1)/TEMP
          RFAC(JW-1)=SVAL*RFAC(IDIAG)
          RFAC(JW)=CVAL*RFAC(IDIAG)
          RFAC(IDIAG)=TEMP
          IF (JCP .LT. NACT) THEN
              DO 820 J=JCP+1,NACT
              TEMP=SVAL*RFAC(JW+JC)+CVAL*RFAC(JW+JCP)
              RFAC(JW+JCP)=CVAL*RFAC(JW+JC)-SVAL*RFAC(JW+JCP)
              RFAC(JW+JC)=TEMP
  820         JW=JW+J
          END IF
          JDIAG=IDIAG-JC
          DO 830 I=1,N
          IF (I .LT. JC) THEN
              TEMP=RFAC(IDIAG+I)
              RFAC(IDIAG+I)=RFAC(JDIAG+I)
              RFAC(JDIAG+I)=TEMP
          END IF
          TEMP=SVAL*QFAC(I,JC)+CVAL*QFAC(I,JCP)
          QFAC(I,JCP)=CVAL*QFAC(I,JC)-SVAL*QFAC(I,JCP)
  830     QFAC(I,JC)=TEMP
          IACT(JC)=IACT(JCP)
          RESACT(JC)=RESACT(JCP)
          VLAM(JC)=VLAM(JCP)
          JC=JCP
          GOTO 810
      END IF
      NACT=NACT-1
      GOTO (50,60,280),IFLAG
      END
      SUBROUTINE QMSTEP (N,NPT,M,AMAT,B,XPT,XOPT,NACT,IACT,
     1  RESCON,QFAC,KOPT,KNEW,DEL,STEP,GL,PQW,RSTAT,W,IFEAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AMAT(N,*),B(*),XPT(NPT,*),XOPT(*),IACT(*),
     1  RESCON(*),QFAC(N,*),STEP(*),GL(*),PQW(*),RSTAT(*),W(*)
C
C     N, NPT, M, AMAT, B, XPT, XOPT, NACT, IACT, RESCON, QFAC, KOPT are the
C       same as the terms with these names in SUBROUTINE LINCOB.
C     KNEW is the index of the interpolation point that is going to be moved.
C     DEL is the current restriction on the length of STEP, which is never
C       greater than the current trust region radius DELTA.
C     STEP will be set to the required step from XOPT to the new point.
C     GL must be set on entry to the gradient of LFUNC at XBASE, where LFUNC
C       is the KNEW-th Lagrange function. It is used also for some other
C       gradients of LFUNC.
C     PQW provides the second derivative parameters of LFUNC.
C     RSTAT and W are used for working space. Their lengths must be at least
C       M and N, respectively. RSTAT(J) is set to -1.0, 0.0, or 1.0 if the
C       J-th constraint is irrelevant, active, or both inactive and relevant,
C       respectively.
C     IFEAS will be set to 0 or 1 if XOPT+STEP is infeasible or feasible.
C
C     STEP is chosen to provide a relatively large value of the modulus of
C       LFUNC(XOPT+STEP), subject to ||STEP|| .LE. DEL. A projected STEP is
C       calculated too, within the trust region, that does not alter the
C       residuals of the active constraints. The projected step is preferred
C       if its value of | LFUNC(XOPT+STEP) | is at least one fifth of the
C       original one, but the greatest violation of a linear constraint must
C       be at least 0.2*DEL, in order to keep the interpolation points apart.
C       The remedy when the maximum constraint violation is too small is to
C       restore the original step, which is perturbed if necessary so that
C       its maximum constraint violation becomes 0.2*DEL.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      TENTH=0.1D0
      ZERO=0.0D0
      TEST=0.2D0*DEL
C
C     Replace GL by the gradient of LFUNC at the trust region centre, and
C       set the elements of RSTAT.
C
      DO 20 K=1,NPT
      TEMP=ZERO
      DO 10 J=1,N
   10 TEMP=TEMP+XPT(K,J)*XOPT(J)
      TEMP=PQW(K)*TEMP
      DO 20 I=1,N
   20 GL(I)=GL(I)+TEMP*XPT(K,I)
      IF (M .GT. 0) THEN
          DO 30 J=1,M
          RSTAT(J)=ONE
   30     IF (DABS(RESCON(J)) .GE. DEL) RSTAT(J)=-ONE
          DO 40 K=1,NACT
   40     RSTAT(IACT(K))=ZERO
      END IF
C
C     Find the greatest modulus of LFUNC on a line through XOPT and
C       another interpolation point within the trust region.
C
      IFLAG=0
      VBIG=ZERO
      DO 60 K=1,NPT
      IF (K .EQ. KOPT) GOTO 60
      SS=ZERO
      SP=ZERO
      DO 50 I=1,N
      TEMP=XPT(K,I)-XOPT(I)
      SS=SS+TEMP*TEMP
   50 SP=SP+GL(I)*TEMP
      STP=-DEL/DSQRT(SS)
      IF (K .EQ. KNEW) THEN
          IF (SP*(SP-ONE) .LT. ZERO) STP=-STP
          VLAG=DABS(STP*SP)+STP*STP*DABS(SP-ONE)
      ELSE
          VLAG=DABS(STP*(ONE-STP)*SP)
      END IF
      IF (VLAG .GT. VBIG) THEN
          KSAV=K
          STPSAV=STP
          VBIG=VLAG
      END IF
   60 CONTINUE
C
C     Set STEP to the move that gives the greatest modulus calculated above.
C       This move may be replaced by a steepest ascent step from XOPT.
C
      GG=ZERO
      DO 70 I=1,N
      GG=GG+GL(I)**2
   70 STEP(I)=STPSAV*(XPT(KSAV,I)-XOPT(I))
      VGRAD=DEL*DSQRT(GG)
      IF (VGRAD .LE. TENTH*VBIG) GOTO 220
C
C     Make the replacement if it provides a larger value of VBIG.
C
      GHG=ZERO
      DO 90 K=1,NPT
      TEMP=ZERO
      DO 80 J=1,N
   80 TEMP=TEMP+XPT(K,J)*GL(J)
   90 GHG=GHG+PQW(K)*TEMP*TEMP
      VNEW=VGRAD+DABS(HALF*DEL*DEL*GHG/GG)
      IF (VNEW .GT. VBIG) THEN
          VBIG=VNEW
          STP=DEL/DSQRT(GG)
          IF (GHG .LT. ZERO) STP=-STP
          DO 100 I=1,N
  100     STEP(I)=STP*GL(I)
      END IF
      IF (NACT .EQ. 0 .OR. NACT .EQ. N) GOTO 220
C
C     Overwrite GL by its projection. Then set VNEW to the greatest
C       value of |LFUNC| on the projected gradient from XOPT subject to
C       the trust region bound. If VNEW is sufficiently large, then STEP
C       may be changed to a move along the projected gradient.
C
      DO 110 K=NACT+1,N
      W(K)=ZERO
      DO 110 I=1,N
  110 W(K)=W(K)+GL(I)*QFAC(I,K)
      GG=ZERO
      DO 130 I=1,N
      GL(I)=ZERO
      DO 120 K=NACT+1,N
  120 GL(I)=GL(I)+QFAC(I,K)*W(K)
  130 GG=GG+GL(I)**2
      VGRAD=DEL*DSQRT(GG)
      IF (VGRAD .LE. TENTH*VBIG) GOTO 220
      GHG=ZERO
      DO 150 K=1,NPT
      TEMP=ZERO
      DO 140 J=1,N
  140 TEMP=TEMP+XPT(K,J)*GL(J)
  150 GHG=GHG+PQW(K)*TEMP*TEMP
      VNEW=VGRAD+DABS(HALF*DEL*DEL*GHG/GG)
C
C     Set W to the possible move along the projected gradient.
C
      STP=DEL/DSQRT(GG)
      IF (GHG .LT. ZERO) STP=-STP
      WW=ZERO
      DO 160 I=1,N
      W(I)=STP*GL(I)
  160 WW=WW+W(I)**2
C
C     Set STEP to W if W gives a sufficiently large value of the modulus
C       of the Lagrange function, and if W either preserves feasibility
C       or gives a constraint violation of at least 0.2*DEL. The purpose
C       of CTOL below is to provide a check on feasibility that includes
C       a tolerance for contributions from computer rounding errors.
C
      IF (VNEW/VBIG .GE. 0.2D0) THEN
          IFEAS=1
          BIGV=ZERO
          J=0
  170     J=J+1
          IF (J .LE. M) THEN
              IF (RSTAT(J) .EQ. ONE) THEN
                  TEMP=-RESCON(J)
                  DO 180 I=1,N
  180             TEMP=TEMP+W(I)*AMAT(I,J)
                  BIGV=DMAX1(BIGV,TEMP)
              END IF
              IF (BIGV .LT. TEST) GOTO 170
              IFEAS=0
          END IF
          CTOL=ZERO
          TEMP=0.01D0*DSQRT(WW)
          IF (BIGV .GT. ZERO .AND. BIGV .LT. TEMP) THEN
              DO 200 K=1,NACT
              J=IACT(K)
              SUM=ZERO
              DO 190 I=1,N
  190         SUM=SUM+W(I)*AMAT(I,J)
  200         CTOL=DMAX1(CTOL,DABS(SUM))
          END IF
          IF (BIGV .LE. 10.0D0*CTOL .OR. BIGV .GE. TEST) THEN
              DO 210 I=1,N
  210         STEP(I)=W(I)
              GOTO 260
          END IF
      END IF
C
C     Calculate the greatest constraint violation at XOPT+STEP with STEP at
C       its original value. Modify STEP if this violation is unacceptable.
C
  220 IFEAS=1
      BIGV=ZERO
      RESMAX=ZERO
      J=0
  230 J=J+1
      IF (J .LE. M) THEN
          IF (RSTAT(J) .LT. ZERO) GOTO 230
          TEMP=-RESCON(J)
          DO 240 I=1,N
  240     TEMP=TEMP+STEP(I)*AMAT(I,J)
          RESMAX=DMAX1(RESMAX,TEMP)
          IF (TEMP .LT. TEST) THEN
              IF (TEMP .LE. BIGV) GOTO 230
              BIGV=TEMP
              JSAV=J
              IFEAS=-1
              GOTO 230
          END IF
          IFEAS=0
      END IF
      IF (IFEAS .EQ. -1) THEN
          DO 250 I=1,N
  250     STEP(I)=STEP(I)+(TEST-BIGV)*AMAT(I,JSAV)
          IFEAS=0
      END IF
C
C     Return the calculated STEP and the value of IFEAS.
C
  260 RETURN
      END
      SUBROUTINE LINUPDATE (N,NPT,XPT,BMAT,ZMAT,IDZ,NDIM,SP,STEP,
     1  KOPT,KNEW,VLAG,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),SP(*),STEP(*),
     2  VLAG(*),W(*)
C
C     The arguments N, NPT, XPT, BMAT, ZMAT, IDZ, NDIM ,SP and STEP are
C       identical to the corresponding arguments in SUBROUTINE LINCOB.
C     KOPT is such that XPT(KOPT,.) is the current trust region centre.
C     KNEW on exit is usually positive, and then it is the index of an
C       interpolation point to be moved to the position XPT(KOPT,.)+STEP(.).
C       It is set on entry either to its final value or to 0. In the latter
C       case, the final value of KNEW is chosen to maximize the denominator
C       of the matrix updating formula times a weighting factor.
C     VLAG and W are used for working space, the first NPT+N elements of
C       both of these vectors being required.
C
C     The arrays BMAT and ZMAT with IDZ are updated, the new matrices being
C       the ones that are suitable after the shift of the KNEW-th point to
C       the new position XPT(KOPT,.)+STEP(.). A return with KNEW set to zero
C       occurs if the calculation fails due to a zero denominator in the
C       updating formula, which should never happen.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
C
C     Calculate VLAG and BETA for the current choice of STEP. The first NPT
C       elements of VLAG are set to the values of the Lagrange functions at
C       XPT(KOPT,.)+STEP(.). The first NPT components of W_check are held
C       in W, where W_check is defined in a paper on the updating method.
C
      DO 20 K=1,NPT
      W(K)=SP(NPT+K)*(HALF*SP(NPT+K)+SP(K))
      SUM=ZERO
      DO 10 J=1,N
   10 SUM=SUM+BMAT(K,J)*STEP(J)
   20 VLAG(K)=SUM
      BETA=ZERO
      DO 40 K=1,NPTM
      SUM=ZERO
      DO 30 I=1,NPT
   30 SUM=SUM+ZMAT(I,K)*W(I)
      IF (K .LT. IDZ) THEN
          BETA=BETA+SUM*SUM
          SUM=-SUM
      ELSE
          BETA=BETA-SUM*SUM
      END IF
      DO 40 I=1,NPT
   40 VLAG(I)=VLAG(I)+SUM*ZMAT(I,K)
      BSUM=ZERO
      DX=ZERO
      SSQ=ZERO
      DO 70 J=1,N
      SUM=ZERO
      DO 50 I=1,NPT
   50 SUM=SUM+W(I)*BMAT(I,J)
      BSUM=BSUM+SUM*STEP(J)
      JP=NPT+J
      DO 60 K=1,N
   60 SUM=SUM+BMAT(JP,K)*STEP(K)
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*STEP(J)
      DX=DX+STEP(J)*XPT(KOPT,J)
   70 SSQ=SSQ+STEP(J)**2
      BETA=DX*DX+SSQ*(SP(KOPT)+DX+DX+HALF*SSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
C
C     If KNEW is zero initially, then pick the index of the interpolation
C       point to be deleted, by maximizing the absolute value of the
C       denominator of the updating formula times a weighting factor.
C       
C
      IF (KNEW .EQ. 0) THEN
          DENMAX=ZERO
          DO 100 K=1,NPT
          HDIAG=ZERO
          DO 80 J=1,NPTM
          TEMP=ONE
          IF (J .LT. IDZ) TEMP=-ONE
   80     HDIAG=HDIAG+TEMP*ZMAT(K,J)**2
          DENABS=DABS(BETA*HDIAG+VLAG(K)**2)
          DISTSQ=ZERO
          DO 90 J=1,N
   90     DISTSQ=DISTSQ+(XPT(K,J)-XPT(KOPT,J))**2
          TEMP=DENABS*DISTSQ*DISTSQ
          IF (TEMP .GT. DENMAX) THEN
              DENMAX=TEMP
              KNEW=K
          END IF
  100     CONTINUE
      END IF
C
C     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
C
      JL=1
      IF (NPTM .GE. 2) THEN
          DO 120 J=2,NPTM
          IF (J .EQ. IDZ) THEN
              JL=IDZ
          ELSE IF (ZMAT(KNEW,J) .NE. ZERO) THEN
              TEMP=DSQRT(ZMAT(KNEW,JL)**2+ZMAT(KNEW,J)**2)
              TEMPA=ZMAT(KNEW,JL)/TEMP
              TEMPB=ZMAT(KNEW,J)/TEMP
              DO 110 I=1,NPT
              TEMP=TEMPA*ZMAT(I,JL)+TEMPB*ZMAT(I,J)
              ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,JL)
  110         ZMAT(I,JL)=TEMP
              ZMAT(KNEW,J)=ZERO
          END IF
  120     CONTINUE
      END IF
C
C     Put the first NPT components of the KNEW-th column of the Z Z^T matrix
C       into W, and calculate the parameters of the updating formula.
C
      TEMPA=ZMAT(KNEW,1)
      IF (IDZ .GE. 2) TEMPA=-TEMPA
      IF (JL .GT. 1) TEMPB=ZMAT(KNEW,JL)
      DO 130 I=1,NPT
      W(I)=TEMPA*ZMAT(I,1)
      IF (JL .GT. 1) W(I)=W(I)+TEMPB*ZMAT(I,JL)
  130 CONTINUE
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      TAUSQ=TAU*TAU
      DENOM=ALPHA*BETA+TAUSQ
      VLAG(KNEW)=VLAG(KNEW)-ONE
      IF (DENOM .EQ. ZERO) THEN
          KNEW=0
          GOTO 180
      END IF
      SQRTDN=DSQRT(DABS(DENOM))
C
C     Complete the updating of ZMAT when there is only one nonzero element
C       in the KNEW-th row of the new matrix ZMAT. IFLAG is set to one when
C       the value of IDZ is going to be reduced.
C
      IFLAG=0
      IF (JL .EQ. 1) THEN
          TEMPA=TAU/SQRTDN
          TEMPB=ZMAT(KNEW,1)/SQRTDN
          DO 140 I=1,NPT
  140     ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
          IF (DENOM .LT. ZERO) THEN
              IF (IDZ .EQ. 1) THEN
                  IDZ=2
              ELSE
                  IFLAG=1
              END IF
          END IF
      ELSE
C
C     Complete the updating of ZMAT in the alternative case.
C
          JA=1
          IF (BETA .GE. ZERO) JA=JL
          JB=JL+1-JA
          TEMP=ZMAT(KNEW,JB)/DENOM
          TEMPA=TEMP*BETA
          TEMPB=TEMP*TAU
          TEMP=ZMAT(KNEW,JA)
          SCALA=ONE/DSQRT(DABS(BETA)*TEMP*TEMP+TAUSQ)
          SCALB=SCALA*SQRTDN
          DO 150 I=1,NPT
          ZMAT(I,JA)=SCALA*(TAU*ZMAT(I,JA)-TEMP*VLAG(I))
  150     ZMAT(I,JB)=SCALB*(ZMAT(I,JB)-TEMPA*W(I)-TEMPB*VLAG(I))
          IF (DENOM .LE. ZERO) THEN
              IF (BETA .LT. ZERO) THEN
                  IDZ=IDZ+1
              ELSE
                  IFLAG=1
              END IF
          END IF
      END IF
C
C     Reduce IDZ when the diagonal part of the ZMAT times Diag(DZ) times
C       ZMAT^T factorization gains another positive element. Then exchange
C       the first and IDZ-th columns of ZMAT.
C
      IF (IFLAG .EQ. 1) THEN
          IDZ=IDZ-1
          DO 160 I=1,NPT
          TEMP=ZMAT(I,1)
          ZMAT(I,1)=ZMAT(I,IDZ)
  160     ZMAT(I,IDZ)=TEMP
      END IF
C
C     Finally, update the matrix BMAT.
C
      DO 170 J=1,N
      JP=NPT+J
      W(JP)=BMAT(KNEW,J)
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
      DO 170 I=1,JP
      BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
      IF (I .GT. NPT) BMAT(JP,I-NPT)=BMAT(I,J)
  170 CONTINUE
  180 RETURN
      END
      SUBROUTINE LINPRELIM (N,NPT,M,AMAT,B,X,RHOBEG,IPRINT,XBASE,
     1  XPT,FVAL,XSAV,XOPT,GOPT,KOPT,HQ,PQ,BMAT,ZMAT,IDZ,NDIM,
     2  SP,RESCON,STEP,PQW,W,CALFUN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AMAT(N,*),B(*),X(*),XBASE(*),XPT(NPT,*),FVAL(*),
     1  XSAV(*),XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),
     2  SP(*),RESCON(*),STEP(*),PQW(*),W(*)
      EXTERNAL CALFUN
C
C     The arguments N, NPT, M, AMAT, B, X, RHOBEG, IPRINT, XBASE, XPT, FVAL,
C       XSAV, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SP and RESCON are the
C       same as the corresponding arguments in SUBROUTINE LINCOB.
C     KOPT is set to the integer such that XPT(KOPT,.) is the initial trust
C       region centre.
C     IDZ is going to be set to one, so that every element of Diag(DZ) is
C       one in the product ZMAT times Diag(DZ) times ZMAT^T, which is the
C       factorization of the leading NPT by NPT submatrix of H.
C     STEP, PQW and W are used for working space, the arrays STEP and PQW
C       being taken from LINCOB. The length of W must be at least N+NPT.
C
C     SUBROUTINE PRELIM provides the elements of XBASE, XPT, BMAT and ZMAT
C       for the first iteration, an important feature being that, if any of
C       of the columns of XPT is an infeasible point, then the largest of
C       the constraint violations there is at least 0.2*RHOBEG. It also sets
C       the initial elements of FVAL, XOPT, GOPT, HQ, PQ, SP and RESCON.
C
C     Set some constants.
C
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      RECIQ=DSQRT(HALF)/RHOSQ
      TEST=0.2D0*RHOBEG
      IDZ=1
      KBASE=1
C
C     Set the initial elements of XPT, BMAT, SP and ZMAT to zero. 
C
      DO 20 J=1,N
      XBASE(J)=X(J)
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 K=1,NPT
      SP(K)=ZERO
      DO 30 J=1,NPT-N-1
   30 ZMAT(K,J)=ZERO
C
C     Set the nonzero coordinates of XPT(K,.), K=1,2,...,min[2*N+1,NPT],
C       but they may be altered later to make a constraint violation
C       sufficiently large. The initial nonzero elements of BMAT and of
C       the first min[N,NPT-N-1] columns of ZMAT are set also.
C
      DO 40 J=1,N
      XPT(J+1,J)=RHOBEG
      IF (J .LT. NPT-N) THEN
          JP=N+J+1
          XPT(JP,J)=-RHOBEG
          BMAT(J+1,J)=HALF/RHOBEG
          BMAT(JP,J)=-HALF/RHOBEG
          ZMAT(1,J)=-RECIQ-RECIQ
          ZMAT(J+1,J)=RECIQ
          ZMAT(JP,J)=RECIQ
      ELSE
          BMAT(1,J)=-ONE/RHOBEG
          BMAT(J+1,J)=ONE/RHOBEG
          BMAT(NPT+J,J)=-HALF*RHOSQ
      END IF
   40 CONTINUE
C
C     Set the remaining initial nonzero elements of XPT and ZMAT when the
C       number of interpolation points exceeds 2*N+1.
C
      IF (NPT .GT. 2*N+1) THEN
          DO 50 K=N+1,NPT-N-1
          ITEMP=(K-1)/N
          IPT=K-ITEMP*N
          JPT=IPT+ITEMP
          IF (JPT .GT. N) JPT=JPT-N
          XPT(N+K+1,IPT)=RHOBEG
          XPT(N+K+1,JPT)=RHOBEG
          ZMAT(1,K)=RECIP
          ZMAT(IPT+1,K)=-RECIP
          ZMAT(JPT+1,K)=-RECIP
   50     ZMAT(N+K+1,K)=RECIP
      END IF
C
C     Update the constraint right hand sides to allow for the shift XBASE.
C
      IF (M .GT. 0) THEN
          DO 70 J=1,M
          TEMP=ZERO
          DO 60 I=1,N
   60     TEMP=TEMP+AMAT(I,J)*XBASE(I)
   70     B(J)=B(J)-TEMP
      END IF
C
C     Go through the initial points, shifting every infeasible point if
C       necessary so that its constraint violation is at least 0.2*RHOBEG.
C
      DO 150 NF=1,NPT
      FEAS=ONE
      BIGV=ZERO
      J=0
   80 J=J+1
      IF (J .LE. M .AND. NF .GE. 2) THEN
          RESID=-B(J)
          DO 90 I=1,N
   90     RESID=RESID+XPT(NF,I)*AMAT(I,J)
          IF (RESID .LE. BIGV) GOTO 80
          BIGV=RESID
          JSAV=J
          IF (RESID .LE. TEST) THEN
              FEAS=-ONE
              GOTO 80
          END IF
          FEAS=ZERO
      END IF
      IF (FEAS .LT. ZERO) THEN
          DO 100 I=1,N
  100     STEP(I)=XPT(NF,I)+(TEST-BIGV)*AMAT(I,JSAV)
          DO 110 K=1,NPT
          SP(NPT+K)=ZERO
          DO 110 J=1,N
  110     SP(NPT+K)=SP(NPT+K)+XPT(K,J)*STEP(J)
          CALL LINUPDATE (N,NPT,XPT,BMAT,ZMAT,IDZ,NDIM,SP,STEP,
     1      KBASE,NF,PQW,W)
          DO 120 I=1,N
  120     XPT(NF,I)=STEP(I)
      END IF
C
C     Calculate the objective function at the current interpolation point,
C       and set KOPT to the index of the first trust region centre.
C
      DO 130 J=1,N
  130 X(J)=XBASE(J)+XPT(NF,J)
      F=FEAS
      CALL CALFUN (N,X,F)
      IF (IPRINT .EQ. 3) THEN
          WRITE (IPRINT,140) NF,F,(X(I),I=1,N)
  140     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
     1      '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NF .EQ. 1) THEN
          KOPT=1
      ELSE IF (F .LT. FVAL(KOPT) .AND. FEAS .GT. ZERO) THEN
          KOPT=NF
      END IF
  150 FVAL(NF)=F
C
C     Set PQ for the first quadratic model.
C
      DO 160 J=1,NPTM
      W(J)=ZERO
      DO 160 K=1,NPT
  160 W(J)=W(J)+ZMAT(K,J)*FVAL(K)
      DO 170 K=1,NPT
      PQ(K)=ZERO
      DO 170 J=1,NPTM
  170 PQ(K)=PQ(K)+ZMAT(K,J)*W(J)
C
C     Set XOPT, SP, GOPT and HQ for the first quadratic model.
C
      DO 180 J=1,N
      XOPT(J)=XPT(KOPT,J)
      XSAV(J)=XBASE(J)+XOPT(J)
  180 GOPT(J)=ZERO
      DO 200 K=1,NPT
      SP(K)=ZERO
      DO 190 J=1,N
  190 SP(K)=SP(K)+XPT(K,J)*XOPT(J)
      TEMP=PQ(K)*SP(K)
      DO 200 J=1,N
  200 GOPT(J)=GOPT(J)+FVAL(K)*BMAT(K,J)+TEMP*XPT(K,J)
      DO 210 I=1,(N*N+N)/2
  210 HQ(I)=ZERO
C
C     Set the initial elements of RESCON.
C
      DO 230 J=1,M
      TEMP=B(J)
      DO 220 I=1,N
  220 TEMP=TEMP-XOPT(I)*AMAT(I,J)
      TEMP=DMAX1(TEMP,ZERO)
      IF (TEMP .GE. RHOBEG) TEMP=-TEMP
  230 RESCON(J)=TEMP  
      RETURN
      END
      SUBROUTINE COBYLA (N,M,X,RHOBEG,RHOEND,LSTDOUT,MAXFUN,W,IACT,
     1  CALCFC,OUTFILE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*),IACT(*)
C
C     This subroutine minimizes an objective function F(X) subject to M
C     inequality constraints on X, where X is a vector of variables that has
C     N components. The algorithm employs linear approximations to the
C     objective and constraint functions, the approximations being formed by
C     linear interpolation at N+1 points in the space of the variables.
C     We regard these interpolation points as vertices of a simplex. The
C     parameter RHO controls the size of the simplex and it is reduced
C     automatically from RHOBEG to RHOEND. For each RHO the subroutine tries
C     to achieve a good vector of variables for the current size, and then
C     RHO is reduced until the value RHOEND is reached. Therefore RHOBEG and
C     RHOEND should be set to reasonable initial changes to and the required   
C     accuracy in the variables respectively, but this accuracy should be
C     viewed as a subject for experimentation because it is not guaranteed.
C     The subroutine has an advantage over many of its competitors, however,
C     which is that it treats each constraint individually when calculating
C     a change to the variables, instead of lumping the constraints together
C     into a single penalty function. The name of the subroutine is derived
C     from the phrase Constrained Optimization BY Linear Approximations.
C
C     The user must set the values of N, M, RHOBEG and RHOEND, and must
C     provide an initial vector of variables in X. Further, the value of
C     IPRINT should be set to 0, 1, 2 or 3, which controls the amount of
C     printing during the calculation. Specifically, there is no output if
C     IPRINT=0 and there is output only at the end of the calculation if
C     IPRINT=1. Otherwise each new value of RHO and SIGMA is printed.
C     Further, the vector of variables and some function information are
C     given either when RHO is reduced or when each new value of F(X) is
C     computed in the cases IPRINT=2 or IPRINT=3 respectively. Here SIGMA
C     is a penalty parameter, it being assumed that a change to X is an
C     improvement if it reduces the merit function
C                F(X)+SIGMA*MAX(0.0,-C1(X),-C2(X),...,-CM(X)),
C     where C1,C2,...,CM denote the constraint functions that should become
C     nonnegative eventually, at least to the precision of RHOEND. In the
C     printed output the displayed term that is multiplied by SIGMA is
C     called MAXCV, which stands for 'MAXimum Constraint Violation'. The
C     argument MAXFUN is an integer variable that must be set by the user to a
C     limit on the number of calls of CALCFC, the purpose of this routine being
C     given below. The value of MAXFUN will be altered to the number of calls
C     of CALCFC that are made. The arguments W and IACT provide real and
C     integer arrays that are used as working space. Their lengths must be at
C     least N*(3*N+2*M+11)+4*M+6 and M+1 respectively.
C
C     In order to define the objective and constraint functions, we require
C     a subroutine that has the name and arguments
C                SUBROUTINE CALCFC (N,M,X,F,CON)
C                DIMENSION X(*),CON(*)  .
C     The values of N and M are fixed and have been defined already, while
C     X is now the current vector of variables. The subroutine should return
C     the objective and constraint functions at X in F and CON(1),CON(2),
C     ...,CON(M). Note that we are trying to adjust X so that F(X) is as
C     small as possible subject to the constraint functions being nonnegative.
C
C     Partition the working space array W to provide the storage that is needed
C     for the main calculation.
C
C     Definition of stream control variables 256 appears in calling codes
      EXTERNAL      CALCFC
      INTEGER       LSTDOUT
      CHARACTER*256 OUTFILE
      IF (LSTDOUT .gt. 1) THEN
          IPRINT = 12
          OPEN (IPRINT, FILE=OUTFILE, STATUS='OLD')
      ELSE
          IPRINT = 0
      ENDIF
C     End of it. KMK
      MPP=M+2
      ICON=1
      ISIM=ICON+MPP
      ISIMI=ISIM+N*N+N
      IDATM=ISIMI+N*N
      IA=IDATM+N*MPP+MPP
      IVSIG=IA+M*N+N
      IVETA=IVSIG+N
      ISIGB=IVETA+N
      IDX=ISIGB+N
      IWORK=IDX+N
      CALL COBYLB (N,M,MPP,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(ICON),
     1  W(ISIM),W(ISIMI),W(IDATM),W(IA),W(IVSIG),W(IVETA),W(ISIGB),
     2  W(IDX),W(IWORK),IACT,CALCFC)
      RETURN
      END
      SUBROUTINE COBYLB (N,M,MPP,X,RHOBEG,RHOEND,IPRINT,MAXFUN,
     1  CON,SIM,SIMI,DATMAT,A,VSIG,VETA,SIGBAR,DX,W,IACT,CALCFC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),CON(*),SIM(N,*),SIMI(N,*),DATMAT(MPP,*),
     1  A(N,*),VSIG(*),VETA(*),SIGBAR(*),DX(*),W(*),IACT(*)
      EXTERNAL      CALCFC
C
C     Set the initial values of some parameters. The last column of SIM holds
C     the optimal vertex of the current simplex, and the preceding N columns
C     hold the displacements from the optimal vertex to the other vertices.
C     Further, SIMI holds the inverse of the matrix that is contained in the
C     first N columns of SIM.
C
      IPTEM=MIN0(N,5)
      IPTEMP=IPTEM+1
      NP=N+1
      MP=M+1
      ALPHA=0.25
      BETA=2.1
      GAMMA=0.5
      DELTA=1.1
      RHO=RHOBEG
      PARMU=0.0
      IF (IPRINT .GE. 2) WRITE(IPRINT,10) RHO
   10 FORMAT (/3X,'The initial value of RHO is',1PE13.6,2X,
     1  'and PARMU is set to zero.')
      NFVALS=0
      TEMP=1.0/RHO
      DO 30 I=1,N
      SIM(I,NP)=X(I)
      DO 20 J=1,N
      SIM(I,J)=0.0
   20 SIMI(I,J)=0.0
      SIM(I,I)=RHO
   30 SIMI(I,I)=TEMP
      JDROP=NP
      IBRNCH=0
C
C     Make the next call of the user-supplied subroutine CALCFC. These
C     instructions are also used for calling CALCFC during the iterations of
C     the algorithm.
C
   40 IF (NFVALS .GE. MAXFUN .AND. NFVALS .GT. 0) THEN
          IF (IPRINT .GE. 1) WRITE(IPRINT,50)
   50     FORMAT (/3X,'Return from subroutine COBYLA because the ',
     1      'MAXFUN limit has been reached.')
          GOTO 600
      END IF
      NFVALS=NFVALS+1
      CALL CALCFC (N,M,X,F,CON)
      RESMAX=0.0
      IF (M .GT. 0) THEN
          DO 60 K=1,M
   60     RESMAX=AMAX1(RESMAX,-CON(K))
      END IF
      IF (NFVALS .EQ. IPRINT-1 .OR. IPRINT .GE. 3) THEN
          WRITE(IPRINT,70) NFVALS,F,RESMAX,(X(I),I=1,IPTEM)
   70     FORMAT (/3X,'NFVALS =',I5,3X,'F =',1PE13.6,4X,'MAXCV =',
     1      1PE13.6/3X,'X =',1PE13.6,1P4E15.6)
          IF (IPTEM .LT. N) WRITE(IPRINT,80) (X(I),I=IPTEMP,N)
   80     FORMAT (1PE19.6,1P4E15.6)
      END IF
      CON(MP)=F
      CON(MPP)=RESMAX
      IF (IBRNCH .EQ. 1) GOTO 440
C
C     Set the recently calculated function values in a column of DATMAT. This
C     array has a column for each vertex of the current simplex, the entries of
C     each column being the values of the constraint functions (if any)
C     followed by the objective function and the greatest constraint violation
C     at the vertex.
C
      DO 90 K=1,MPP
   90 DATMAT(K,JDROP)=CON(K)
      IF (NFVALS .GT. NP) GOTO 130
C
C     Exchange the new vertex of the initial simplex with the optimal vertex if
C     necessary. Then, if the initial simplex is not complete, pick its next
C     vertex and calculate the function values there.
C
      IF (JDROP .LE. N) THEN
          IF (DATMAT(MP,NP) .LE. F) THEN
              X(JDROP)=SIM(JDROP,NP)
          ELSE
              SIM(JDROP,NP)=X(JDROP)
              DO 100 K=1,MPP
              DATMAT(K,JDROP)=DATMAT(K,NP)
  100         DATMAT(K,NP)=CON(K)
              DO 120 K=1,JDROP
              SIM(JDROP,K)=-RHO
              TEMP=0.0
              DO 110 I=K,JDROP
  110         TEMP=TEMP-SIMI(I,K)
  120         SIMI(JDROP,K)=TEMP
          END IF
      END IF
      IF (NFVALS .LE. N) THEN
          JDROP=NFVALS
          X(JDROP)=X(JDROP)+RHO
          GOTO 40
      END IF
  130 IBRNCH=1
C
C     Identify the optimal vertex of the current simplex.
C
  140 PHIMIN=DATMAT(MP,NP)+PARMU*DATMAT(MPP,NP)
      NBEST=NP
      DO 150 J=1,N
      TEMP=DATMAT(MP,J)+PARMU*DATMAT(MPP,J)
      IF (TEMP .LT. PHIMIN) THEN
          NBEST=J
          PHIMIN=TEMP
      ELSE IF (TEMP .EQ. PHIMIN .AND. PARMU .EQ. 0.0) THEN
          IF (DATMAT(MPP,J) .LT. DATMAT(MPP,NBEST)) NBEST=J
      END IF
  150 CONTINUE
C
C     Switch the best vertex into pole position if it is not there already,
C     and also update SIM, SIMI and DATMAT.
C
      IF (NBEST .LE. N) THEN
          DO 160 I=1,MPP
          TEMP=DATMAT(I,NP)
          DATMAT(I,NP)=DATMAT(I,NBEST)
  160     DATMAT(I,NBEST)=TEMP
          DO 180 I=1,N
          TEMP=SIM(I,NBEST)
          SIM(I,NBEST)=0.0
          SIM(I,NP)=SIM(I,NP)+TEMP
          TEMPA=0.0
          DO 170 K=1,N
          SIM(I,K)=SIM(I,K)-TEMP
  170     TEMPA=TEMPA-SIMI(K,I)
  180     SIMI(NBEST,I)=TEMPA
      END IF
C
C     Make an error return if SIGI is a poor approximation to the inverse of
C     the leading N by N submatrix of SIG.
C
      ERROR=0.0
      DO 200 I=1,N
      DO 200 J=1,N
      TEMP=0.0
      IF (I .EQ. J) TEMP=TEMP-1.0
      DO 190 K=1,N
  190 TEMP=TEMP+SIMI(I,K)*SIM(K,J)
  200 ERROR=AMAX1(ERROR,ABS(TEMP))
      IF (ERROR .GT. 0.1) THEN
          IF (IPRINT .GE. 1) WRITE(IPRINT,210)
  210     FORMAT (/3X,'Return from subroutine COBYLA because ',
     1      'rounding errors are becoming damaging.')
          GOTO 600
      END IF
C
C     Calculate the coefficients of the linear approximations to the objective
C     and constraint functions, placing minus the objective function gradient
C     after the constraint gradients in the array A. The vector W is used for
C     working space.
C
      DO 240 K=1,MP
      CON(K)=-DATMAT(K,NP)
      DO 220 J=1,N
  220 W(J)=DATMAT(K,J)+CON(K)
      DO 240 I=1,N
      TEMP=0.0
      DO 230 J=1,N
  230 TEMP=TEMP+W(J)*SIMI(J,I)
      IF (K .EQ. MP) TEMP=-TEMP
  240 A(I,K)=TEMP
C
C     Calculate the values of sigma and eta, and set IFLAG=0 if the current
C     simplex is not acceptable.
C
      IFLAG=1
      PARSIG=ALPHA*RHO
      PARETA=BETA*RHO
      DO 260 J=1,N
      WSIG=0.0
      WETA=0.0
      DO 250 I=1,N
      WSIG=WSIG+SIMI(J,I)**2
  250 WETA=WETA+SIM(I,J)**2
      VSIG(J)=1.0/SQRT(WSIG)
      VETA(J)=SQRT(WETA)
      IF (VSIG(J) .LT. PARSIG .OR. VETA(J) .GT. PARETA) IFLAG=0
  260 CONTINUE
C
C     If a new vertex is needed to improve acceptability, then decide which
C     vertex to drop from the simplex.
C
      IF (IBRNCH .EQ. 1 .OR. IFLAG .EQ. 1) GOTO 370
      JDROP=0
      TEMP=PARETA
      DO 270 J=1,N
      IF (VETA(J) .GT. TEMP) THEN
          JDROP=J
          TEMP=VETA(J)
      END IF
  270 CONTINUE
      IF (JDROP .EQ. 0) THEN
          DO 280 J=1,N
          IF (VSIG(J) .LT. TEMP) THEN
              JDROP=J
              TEMP=VSIG(J)
          END IF
  280     CONTINUE
      END IF
C
C     Calculate the step to the new vertex and its sign.
C
      TEMP=GAMMA*RHO*VSIG(JDROP)
      DO 290 I=1,N
  290 DX(I)=TEMP*SIMI(JDROP,I)
      CVMAXP=0.0
      CVMAXM=0.0
      DO 310 K=1,MP
      SUM=0.0
      DO 300 I=1,N
  300 SUM=SUM+A(I,K)*DX(I)
      IF (K .LT. MP) THEN
          TEMP=DATMAT(K,NP)
          CVMAXP=AMAX1(CVMAXP,-SUM-TEMP)
          CVMAXM=AMAX1(CVMAXM,SUM-TEMP)
      END IF
  310 CONTINUE
      DXSIGN=1.0
      IF (PARMU*(CVMAXP-CVMAXM) .GT. SUM+SUM) DXSIGN=-1.0
C
C     Update the elements of SIM and SIMI, and set the next X.
C
      TEMP=0.0
      DO 320 I=1,N
      DX(I)=DXSIGN*DX(I)
      SIM(I,JDROP)=DX(I)
  320 TEMP=TEMP+SIMI(JDROP,I)*DX(I)
      DO 330 I=1,N
  330 SIMI(JDROP,I)=SIMI(JDROP,I)/TEMP
      DO 360 J=1,N
      IF (J .NE. JDROP) THEN
          TEMP=0.0
          DO 340 I=1,N
  340     TEMP=TEMP+SIMI(J,I)*DX(I)
          DO 350 I=1,N
  350     SIMI(J,I)=SIMI(J,I)-TEMP*SIMI(JDROP,I)
      END IF
  360 X(J)=SIM(J,NP)+DX(J)
      GOTO 40
C
C     Calculate DX=x(*)-x(0). Branch if the length of DX is less than 0.5*RHO.
C
  370 IZ=1
      IZDOTA=IZ+N*N
      IVMC=IZDOTA+N
      ISDIRN=IVMC+MP
      IDXNEW=ISDIRN+N
      IVMD=IDXNEW+N
      CALL COBTRSTLP (N,M,A,CON,RHO,DX,IFULL,IACT,W(IZ),W(IZDOTA),
     1  W(IVMC),W(ISDIRN),W(IDXNEW),W(IVMD))
      IF (IFULL .EQ. 0) THEN
          TEMP=0.0
          DO 380 I=1,N
  380     TEMP=TEMP+DX(I)**2
          IF (TEMP .LT. 0.25*RHO*RHO) THEN
              IBRNCH=1
              GOTO 550
          END IF
      END IF
C
C     Predict the change to F and the new maximum constraint violation if the
C     variables are altered from x(0) to x(0)+DX.
C
      RESNEW=0.0
      CON(MP)=0.0
      DO 400 K=1,MP
      SUM=CON(K)
      DO 390 I=1,N
  390 SUM=SUM-A(I,K)*DX(I)
      IF (K .LT. MP) RESNEW=AMAX1(RESNEW,SUM)
  400 CONTINUE
C
C     Increase PARMU if necessary and branch back if this change alters the
C     optimal vertex. Otherwise PREREM and PREREC will be set to the predicted
C     reductions in the merit function and the maximum constraint violation
C     respectively.
C
      BARMU=0.0
      PREREC=DATMAT(MPP,NP)-RESNEW
      IF (PREREC .GT. 0.0) BARMU=SUM/PREREC
      IF (PARMU .LT. 1.5*BARMU) THEN
          PARMU=2.0*BARMU
          IF (IPRINT .GE. 2) WRITE(IPRINT,410) PARMU
  410     FORMAT (/3X,'Increase in PARMU to',1PE13.6)
          PHI=DATMAT(MP,NP)+PARMU*DATMAT(MPP,NP)
          DO 420 J=1,N
          TEMP=DATMAT(MP,J)+PARMU*DATMAT(MPP,J)
          IF (TEMP .LT. PHI) GOTO 140
          IF (TEMP .EQ. PHI .AND. PARMU .EQ. 0.0) THEN
              IF (DATMAT(MPP,J) .LT. DATMAT(MPP,NP)) GOTO 140
          END IF
  420     CONTINUE
      END IF
      PREREM=PARMU*PREREC-SUM
C
C     Calculate the constraint and objective functions at x(*). Then find the
C     actual reduction in the merit function.
C
      DO 430 I=1,N
  430 X(I)=SIM(I,NP)+DX(I)
      IBRNCH=1
      GOTO 40
  440 VMOLD=DATMAT(MP,NP)+PARMU*DATMAT(MPP,NP)
      VMNEW=F+PARMU*RESMAX
      TRURED=VMOLD-VMNEW
      IF (PARMU .EQ. 0.0 .AND. F .EQ. DATMAT(MP,NP)) THEN
          PREREM=PREREC
          TRURED=DATMAT(MPP,NP)-RESMAX
      END IF
C
C     Begin the operations that decide whether x(*) should replace one of the
C     vertices of the current simplex, the change being mandatory if TRURED is
C     positive. Firstly, JDROP is set to the index of the vertex that is to be
C     replaced.
C
      RATIO=0.0
      IF (TRURED .LE. 0.0) RATIO=1.0
      JDROP=0
      DO 460 J=1,N
      TEMP=0.0
      DO 450 I=1,N
  450 TEMP=TEMP+SIMI(J,I)*DX(I)
      TEMP=ABS(TEMP)
      IF (TEMP .GT. RATIO) THEN
          JDROP=J
          RATIO=TEMP
      END IF
  460 SIGBAR(J)=TEMP*VSIG(J)
C
C     Calculate the value of ell.
C
      EDGMAX=DELTA*RHO
      L=0
      DO 480 J=1,N
      IF (SIGBAR(J) .GE. PARSIG .OR. SIGBAR(J) .GE. VSIG(J)) THEN
          TEMP=VETA(J)
          IF (TRURED .GT. 0.0) THEN
              TEMP=0.0
              DO 470 I=1,N
  470         TEMP=TEMP+(DX(I)-SIM(I,J))**2
              TEMP=SQRT(TEMP)
          END IF
          IF (TEMP .GT. EDGMAX) THEN
              L=J
              EDGMAX=TEMP
          END IF
      END IF
  480 CONTINUE
      IF (L .GT. 0) JDROP=L
      IF (JDROP .EQ. 0) GOTO 550
C
C     Revise the simplex by updating the elements of SIM, SIMI and DATMAT.
C
      TEMP=0.0
      DO 490 I=1,N
      SIM(I,JDROP)=DX(I)
  490 TEMP=TEMP+SIMI(JDROP,I)*DX(I)
      DO 500 I=1,N
  500 SIMI(JDROP,I)=SIMI(JDROP,I)/TEMP
      DO 530 J=1,N
      IF (J .NE. JDROP) THEN
          TEMP=0.0
          DO 510 I=1,N
  510     TEMP=TEMP+SIMI(J,I)*DX(I)
          DO 520 I=1,N
  520     SIMI(J,I)=SIMI(J,I)-TEMP*SIMI(JDROP,I)
      END IF
  530 CONTINUE
      DO 540 K=1,MPP
  540 DATMAT(K,JDROP)=CON(K)
C
C     Branch back for further iterations with the current RHO.
C
      IF (TRURED .GT. 0.0 .AND. TRURED .GE. 0.1*PREREM) GOTO 140
  550 IF (IFLAG .EQ. 0) THEN
          IBRNCH=0
          GOTO 140
      END IF
C
C     Otherwise reduce RHO if it is not at its least value and reset PARMU.
C
      IF (RHO .GT. RHOEND) THEN
          RHO=0.5*RHO
          IF (RHO .LE. 1.5*RHOEND) RHO=RHOEND
          IF (PARMU .GT. 0.0) THEN
              DENOM=0.0
              DO 570 K=1,MP
              CMIN=DATMAT(K,NP)
              CMAX=CMIN
              DO 560 I=1,N
              CMIN=AMIN1(CMIN,DATMAT(K,I))
  560         CMAX=AMAX1(CMAX,DATMAT(K,I))
              IF (K .LE. M .AND. CMIN .LT. 0.5*CMAX) THEN
                  TEMP=AMAX1(CMAX,0.0)-CMIN
                  IF (DENOM .LE. 0.0) THEN
                      DENOM=TEMP
                  ELSE
                      DENOM=AMIN1(DENOM,TEMP)
                  END IF
              END IF
  570         CONTINUE
              IF (DENOM .EQ. 0.0) THEN
                  PARMU=0.0
              ELSE IF (CMAX-CMIN .LT. PARMU*DENOM) THEN
                  PARMU=(CMAX-CMIN)/DENOM
              END IF
          END IF
          IF (IPRINT .GE. 2) WRITE(IPRINT,580) RHO,PARMU
  580     FORMAT (/3X,'Reduction in RHO to',1PE13.6,'  and PARMU =',
     1      1PE13.6)
          WRITE(IPRINT,70) NFVALS,DATMAT(MP,NP),DATMAT(MPP,NP),
     1          (SIM(I,NP),I=1,IPTEM)
          IF (IPTEM .LT. N) WRITE(IPRINT,80) (X(I),I=IPTEMP,N)
          GOTO 140
      END IF
C
C     Return the best calculated values of the variables.
C
      IF (IPRINT .GE. 1) WRITE(IPRINT,590)
  590 FORMAT (/3X,'Normal return from subroutine COBYLA')
      IF (IFULL .EQ. 1) GOTO 620
  600 DO 610 I=1,N
  610 X(I)=SIM(I,NP)
      F=DATMAT(MP,NP)
      RESMAX=DATMAT(MPP,NP)
  620 IF (IPRINT .GE. 1) THEN
          WRITE(IPRINT,70) NFVALS,F,RESMAX,(X(I),I=1,IPTEM)
          IF (IPTEM .LT. N) WRITE(IPRINT,80) (X(I),I=IPTEMP,N)
      END IF
      MAXFUN=NFVALS
      RETURN
      END
      SUBROUTINE COBTRSTLP (N,M,A,B,RHO,DX,IFULL,IACT,Z,ZDOTA,
     1  VMULTC,SDIRN,DXNEW,VMULTD)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,*),B(*),DX(*),IACT(*),Z(N,*),ZDOTA(*),
     1  VMULTC(*),SDIRN(*),DXNEW(*),VMULTD(*)
C
C     This subroutine calculates an N-component vector DX by applying the
C     following two stages. In the first stage, DX is set to the shortest
C     vector that minimizes the greatest violation of the constraints
C       A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K=2,3,...,M,
C     subject to the Euclidean length of DX being at most RHO. If its length is
C     strictly less than RHO, then we use the resultant freedom in DX to
C     minimize the objective function
C              -A(1,M+1)*DX(1)-A(2,M+1)*DX(2)-...-A(N,M+1)*DX(N)
C     subject to no increase in any greatest constraint violation. This
C     notation allows the gradient of the objective function to be regarded as
C     the gradient of a constraint. Therefore the two stages are distinguished
C     by MCON .EQ. M and MCON .GT. M respectively. It is possible that a
C     degeneracy may prevent DX from attaining the target length RHO. Then the
C     value IFULL=0 would be set, but usually IFULL=1 on return.
C
C     In general NACT is the number of constraints in the active set and
C     IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT
C     contains a permutation of the remaining constraint indices. Further, Z is
C     an orthogonal matrix whose first NACT columns can be regarded as the
C     result of Gram-Schmidt applied to the active constraint gradients. For
C     J=1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th
C     column of Z with the gradient of the J-th active constraint. DX is the
C     current vector of variables and here the residuals of the active
C     constraints should be zero. Further, the active constraints have
C     nonnegative Lagrange multipliers that are held at the beginning of
C     VMULTC. The remainder of this vector holds the residuals of the inactive
C     constraints at DX, the ordering of the components of VMULTC being in
C     agreement with the permutation of the indices of the constraints that is
C     in IACT. All these residuals are nonnegative, which is achieved by the
C     shift RESMAX that makes the least residual zero.
C
C     Initialize Z and some other variables. The value of RESMAX will be
C     appropriate to DX=0, while ICON will be the index of a most violated
C     constraint if RESMAX is positive. Usually during the first stage the
C     vector SDIRN gives a search direction that reduces all the active
C     constraint violations by one simultaneously.
C
      IFULL=1
      MCON=M
      NACT=0
      RESMAX=0.0
      DO 20 I=1,N
      DO 10 J=1,N
   10 Z(I,J)=0.0
      Z(I,I)=1.0
   20 DX(I)=0.0
      IF (M .GE. 1) THEN
          DO 30 K=1,M
          IF (B(K) .GT. RESMAX) THEN
              RESMAX=B(K)
              ICON=K
          END IF
   30     CONTINUE
          DO 40 K=1,M
          IACT(K)=K
   40     VMULTC(K)=RESMAX-B(K)
      END IF
      IF (RESMAX .EQ. 0.0) GOTO 480
      DO 50 I=1,N
   50 SDIRN(I)=0.0
C
C     End the current stage of the calculation if 3 consecutive iterations
C     have either failed to reduce the best calculated value of the objective
C     function or to increase the number of active constraints since the best
C     value was calculated. This strategy prevents cycling, but there is a
C     remote possibility that it will cause premature termination.
C
   60 OPTOLD=0.0
      ICOUNT=0
   70 IF (MCON .EQ. M) THEN
          OPTNEW=RESMAX
      ELSE
          OPTNEW=0.0
          DO 80 I=1,N
   80     OPTNEW=OPTNEW-DX(I)*A(I,MCON)
      END IF
      IF (ICOUNT .EQ. 0 .OR. OPTNEW .LT. OPTOLD) THEN
          OPTOLD=OPTNEW
          NACTX=NACT
          ICOUNT=3
      ELSE IF (NACT .GT. NACTX) THEN
          NACTX=NACT
          ICOUNT=3
      ELSE
          ICOUNT=ICOUNT-1
          IF (ICOUNT .EQ. 0) GOTO 490
      END IF
C
C     If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
C     the active set. Apply Givens rotations so that the last N-NACT-1 columns
C     of Z are orthogonal to the gradient of the new constraint, a scalar
C     product being set to zero if its nonzero value could be due to computer
C     rounding errors. The array DXNEW is used for working space.
C
      IF (ICON .LE. NACT) GOTO 260
      KK=IACT(ICON)
      DO 90 I=1,N
   90 DXNEW(I)=A(I,KK)
      TOT=0.0
      K=N
  100 IF (K .GT. NACT) THEN
          SP=0.0
          SPABS=0.0
          DO 110 I=1,N
          TEMP=Z(I,K)*DXNEW(I)
          SP=SP+TEMP
  110     SPABS=SPABS+ABS(TEMP)
          ACCA=SPABS+0.1*ABS(SP)
          ACCB=SPABS+0.2*ABS(SP)
          IF (SPABS .GE. ACCA .OR. ACCA .GE. ACCB) SP=0.0
          IF (TOT .EQ. 0.0) THEN
              TOT=SP
          ELSE
              KP=K+1
              TEMP=SQRT(SP*SP+TOT*TOT)
              ALPHA=SP/TEMP
              BETA=TOT/TEMP
              TOT=TEMP
              DO 120 I=1,N
              TEMP=ALPHA*Z(I,K)+BETA*Z(I,KP)
              Z(I,KP)=ALPHA*Z(I,KP)-BETA*Z(I,K)
  120         Z(I,K)=TEMP
          END IF
          K=K-1
          GOTO 100
      END IF
C
C     Add the new constraint if this can be done without a deletion from the
C     active set.
C
      IF (TOT .NE. 0.0) THEN
          NACT=NACT+1
          ZDOTA(NACT)=TOT
          VMULTC(ICON)=VMULTC(NACT)
          VMULTC(NACT)=0.0
          GOTO 210
      END IF
C
C     The next instruction is reached if a deletion has to be made from the
C     active set in order to make room for the new active constraint, because
C     the new constraint gradient is a linear combination of the gradients of
C     the old active constraints. Set the elements of VMULTD to the multipliers
C     of the linear combination. Further, set IOUT to the index of the
C     constraint to be deleted, but branch if no suitable index can be found.
C
      RATIO=-1.0
      K=NACT
  130 ZDOTV=0.0
      ZDVABS=0.0
      DO 140 I=1,N
      TEMP=Z(I,K)*DXNEW(I)
      ZDOTV=ZDOTV+TEMP
  140 ZDVABS=ZDVABS+ABS(TEMP)
      ACCA=ZDVABS+0.1*ABS(ZDOTV)
      ACCB=ZDVABS+0.2*ABS(ZDOTV)
      IF (ZDVABS .LT. ACCA .AND. ACCA .LT. ACCB) THEN
          TEMP=ZDOTV/ZDOTA(K)
          IF (TEMP .GT. 0.0 .AND. IACT(K) .LE. M) THEN
              TEMPA=VMULTC(K)/TEMP
              IF (RATIO .LT. 0.0 .OR. TEMPA .LT. RATIO) THEN
                  RATIO=TEMPA
                  IOUT=K
              END IF
           END IF
          IF (K .GE. 2) THEN
              KW=IACT(K)
              DO 150 I=1,N
  150         DXNEW(I)=DXNEW(I)-TEMP*A(I,KW)
          END IF
          VMULTD(K)=TEMP
      ELSE
          VMULTD(K)=0.0
      END IF
      K=K-1
      IF (K .GT. 0) GOTO 130
      IF (RATIO .LT. 0.0) GOTO 490
C
C     Revise the Lagrange multipliers and reorder the active constraints so
C     that the one to be replaced is at the end of the list. Also calculate the
C     new value of ZDOTA(NACT) and branch if it is not acceptable.
C
      DO 160 K=1,NACT
  160 VMULTC(K)=AMAX1(0.0,VMULTC(K)-RATIO*VMULTD(K))
      IF (ICON .LT. NACT) THEN
          ISAVE=IACT(ICON)
          VSAVE=VMULTC(ICON)
          K=ICON
  170     KP=K+1
          KW=IACT(KP)
          SP=0.0
          DO 180 I=1,N
  180     SP=SP+Z(I,K)*A(I,KW)
          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
          ALPHA=ZDOTA(KP)/TEMP
          BETA=SP/TEMP
          ZDOTA(KP)=ALPHA*ZDOTA(K)
          ZDOTA(K)=TEMP
          DO 190 I=1,N
          TEMP=ALPHA*Z(I,KP)+BETA*Z(I,K)
          Z(I,KP)=ALPHA*Z(I,K)-BETA*Z(I,KP)
  190     Z(I,K)=TEMP
          IACT(K)=KW
          VMULTC(K)=VMULTC(KP)
          K=KP
          IF (K .LT. NACT) GOTO 170
          IACT(K)=ISAVE
          VMULTC(K)=VSAVE
      END IF
      TEMP=0.0
      DO 200 I=1,N
  200 TEMP=TEMP+Z(I,NACT)*A(I,KK)
      IF (TEMP .EQ. 0.0) GOTO 490
      ZDOTA(NACT)=TEMP
      VMULTC(ICON)=0.0
      VMULTC(NACT)=RATIO
C
C     Update IACT and ensure that the objective function continues to be
C     treated as the last active constraint when MCON>M.
C
  210 IACT(ICON)=IACT(NACT)
      IACT(NACT)=KK
      IF (MCON .GT. M .AND. KK .NE. MCON) THEN
          K=NACT-1
          SP=0.0
          DO 220 I=1,N
  220     SP=SP+Z(I,K)*A(I,KK)
          TEMP=SQRT(SP*SP+ZDOTA(NACT)**2)
          ALPHA=ZDOTA(NACT)/TEMP
          BETA=SP/TEMP
          ZDOTA(NACT)=ALPHA*ZDOTA(K)
          ZDOTA(K)=TEMP
          DO 230 I=1,N
          TEMP=ALPHA*Z(I,NACT)+BETA*Z(I,K)
          Z(I,NACT)=ALPHA*Z(I,K)-BETA*Z(I,NACT)
  230     Z(I,K)=TEMP
          IACT(NACT)=IACT(K)
          IACT(K)=KK
          TEMP=VMULTC(K)
          VMULTC(K)=VMULTC(NACT)
          VMULTC(NACT)=TEMP
      END IF
C
C     If stage one is in progress, then set SDIRN to the direction of the next
C     change to the current vector of variables.
C
      IF (MCON .GT. M) GOTO 320
      KK=IACT(NACT)
      TEMP=0.0
      DO 240 I=1,N
  240 TEMP=TEMP+SDIRN(I)*A(I,KK)
      TEMP=TEMP-1.0
      TEMP=TEMP/ZDOTA(NACT)
      DO 250 I=1,N
  250 SDIRN(I)=SDIRN(I)-TEMP*Z(I,NACT)
      GOTO 340
C
C     Delete the constraint that has the index IACT(ICON) from the active set.
C
  260 IF (ICON .LT. NACT) THEN
          ISAVE=IACT(ICON)
          VSAVE=VMULTC(ICON)
          K=ICON
  270     KP=K+1
          KK=IACT(KP)
          SP=0.0
          DO 280 I=1,N
  280     SP=SP+Z(I,K)*A(I,KK)
          TEMP=SQRT(SP*SP+ZDOTA(KP)**2)
          ALPHA=ZDOTA(KP)/TEMP
          BETA=SP/TEMP
          ZDOTA(KP)=ALPHA*ZDOTA(K)
          ZDOTA(K)=TEMP
          DO 290 I=1,N
          TEMP=ALPHA*Z(I,KP)+BETA*Z(I,K)
          Z(I,KP)=ALPHA*Z(I,K)-BETA*Z(I,KP)
  290     Z(I,K)=TEMP
          IACT(K)=KK
          VMULTC(K)=VMULTC(KP)
          K=KP
          IF (K .LT. NACT) GOTO 270
          IACT(K)=ISAVE
          VMULTC(K)=VSAVE
      END IF
      NACT=NACT-1
C
C     If stage one is in progress, then set SDIRN to the direction of the next
C     change to the current vector of variables.
C
      IF (MCON .GT. M) GOTO 320
      TEMP=0.0
      DO 300 I=1,N
  300 TEMP=TEMP+SDIRN(I)*Z(I,NACT+1)
      DO 310 I=1,N
  310 SDIRN(I)=SDIRN(I)-TEMP*Z(I,NACT+1)
      GO TO 340
C
C     Pick the next search direction of stage two.
C
  320 TEMP=1.0/ZDOTA(NACT)
      DO 330 I=1,N
  330 SDIRN(I)=TEMP*Z(I,NACT)
C
C     Calculate the step to the boundary of the trust region or take the step
C     that reduces RESMAX to zero. The two statements below that include the
C     factor 1.0E-6 prevent some harmless underflows that occurred in a test
C     calculation. Further, we skip the step if it could be zero within a
C     reasonable tolerance for computer rounding errors.
C
  340 DD=RHO*RHO
      SD=0.0
      SS=0.0
      DO 350 I=1,N
      IF (ABS(DX(I)) .GE. 1.0E-6*RHO) DD=DD-DX(I)**2
      SD=SD+DX(I)*SDIRN(I)
  350 SS=SS+SDIRN(I)**2
      IF (DD .LE. 0.0) GOTO 490
      TEMP=SQRT(SS*DD)
      IF (ABS(SD) .GE. 1.0E-6*TEMP) TEMP=SQRT(SS*DD+SD*SD)
      STPFUL=DD/(TEMP+SD)
      STEP=STPFUL
      IF (MCON .EQ. M) THEN
          ACCA=STEP+0.1*RESMAX
          ACCB=STEP+0.2*RESMAX
          IF (STEP .GE. ACCA .OR. ACCA .GE. ACCB) GOTO 480
          STEP=AMIN1(STEP,RESMAX)
      END IF
C
C     Set DXNEW to the new variables if STEP is the steplength, and reduce
C     RESMAX to the corresponding maximum residual if stage one is being done.
C     Because DXNEW will be changed during the calculation of some Lagrange
C     multipliers, it will be restored to the following value later.
C
      DO 360 I=1,N
  360 DXNEW(I)=DX(I)+STEP*SDIRN(I)
      IF (MCON .EQ. M) THEN
          RESOLD=RESMAX
          RESMAX=0.0
          DO 380 K=1,NACT
          KK=IACT(K)
          TEMP=B(KK)
          DO 370 I=1,N
  370     TEMP=TEMP-A(I,KK)*DXNEW(I)
          RESMAX=AMAX1(RESMAX,TEMP)
  380     CONTINUE
      END IF
C
C     Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
C     device is included to force VMULTD(K)=0.0 if deviations from this value
C     can be attributed to computer rounding errors. First calculate the new
C     Lagrange multipliers.
C
      K=NACT
  390 ZDOTW=0.0
      ZDWABS=0.0
      DO 400 I=1,N
      TEMP=Z(I,K)*DXNEW(I)
      ZDOTW=ZDOTW+TEMP
  400 ZDWABS=ZDWABS+ABS(TEMP)
      ACCA=ZDWABS+0.1*ABS(ZDOTW)
      ACCB=ZDWABS+0.2*ABS(ZDOTW)
      IF (ZDWABS .GE. ACCA .OR. ACCA .GE. ACCB) ZDOTW=0.0
      VMULTD(K)=ZDOTW/ZDOTA(K)
      IF (K .GE. 2) THEN
          KK=IACT(K)
          DO 410 I=1,N
  410     DXNEW(I)=DXNEW(I)-VMULTD(K)*A(I,KK)
          K=K-1
          GOTO 390
      END IF
      IF (MCON .GT. M) VMULTD(NACT)=AMAX1(0.0,VMULTD(NACT))
C
C     Complete VMULTC by finding the new constraint residuals.
C
      DO 420 I=1,N
  420 DXNEW(I)=DX(I)+STEP*SDIRN(I)
      IF (MCON .GT. NACT) THEN
          KL=NACT+1
          DO 440 K=KL,MCON
          KK=IACT(K)
          SUM=RESMAX-B(KK)
          SUMABS=RESMAX+ABS(B(KK))
          DO 430 I=1,N
          TEMP=A(I,KK)*DXNEW(I)
          SUM=SUM+TEMP
  430     SUMABS=SUMABS+ABS(TEMP)
          ACCA=SUMABS+0.1*ABS(SUM)
          ACCB=SUMABS+0.2*ABS(SUM)
          IF (SUMABS .GE. ACCA .OR. ACCA .GE. ACCB) SUM=0.0
  440     VMULTD(K)=SUM
      END IF
C
C     Calculate the fraction of the step from DX to DXNEW that will be taken.
C
      RATIO=1.0
      ICON=0
      DO 450 K=1,MCON
      IF (VMULTD(K) .LT. 0.0) THEN
          TEMP=VMULTC(K)/(VMULTC(K)-VMULTD(K))
          IF (TEMP .LT. RATIO) THEN
              RATIO=TEMP
              ICON=K
          END IF
      END IF
  450 CONTINUE
C
C     Update DX, VMULTC and RESMAX.
C
      TEMP=1.0-RATIO
      DO 460 I=1,N
  460 DX(I)=TEMP*DX(I)+RATIO*DXNEW(I)
      DO 470 K=1,MCON
  470 VMULTC(K)=AMAX1(0.0,TEMP*VMULTC(K)+RATIO*VMULTD(K))
      IF (MCON .EQ. M) RESMAX=RESOLD+RATIO*(RESMAX-RESOLD)
C
C     If the full step is not acceptable then begin another iteration.
C     Otherwise switch to stage two or end the calculation.
C
      IF (ICON .GT. 0) GOTO 70
      IF (STEP .EQ. STPFUL) GOTO 500
  480 MCON=M+1
      ICON=MCON
      IACT(MCON)=MCON
      VMULTC(MCON)=0.0
      GOTO 60
C
C     We employ any freedom that may be available to reduce the objective
C     function before returning a DX whose length is less than RHO.
C
  490 IF (MCON .EQ. M) GOTO 480
      IFULL=0
  500 RETURN
      END



