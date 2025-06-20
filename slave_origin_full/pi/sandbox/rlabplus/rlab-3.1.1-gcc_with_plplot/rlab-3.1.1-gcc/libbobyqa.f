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
          PRINT 10
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
          PRINT 20
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
   40 RETURN
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
          IF (IPRINT .GT. 0) PRINT 390
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
          IF (IPRINT .GT. 0) PRINT 390
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
              IF (IPRINT .GT. 0) PRINT 320
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
              IF (IPRINT .GT. 0) PRINT 320
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
          IF (IPRINT .GT. 0) PRINT 390
  390     FORMAT (/4X,'Return from BOBYQA because CALFUN has been',
     1      ' called MAXFUN times.')
          GOTO 720
      END IF
      NF=NF+1
      CALL CALFUN (N,X,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 400, NF,F,(X(I),I=1,N)
  400      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
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
              IF (IPRINT .GT. 0) PRINT 430
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
      CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
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
              IF (IPRINT .GE. 3) PRINT 690
  690         FORMAT (5X)
              PRINT 700, RHO,NF
  700         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
     1          ' function values =',I6)
              PRINT 710, FVAL(KOPT),(XBASE(I)+XOPT(I),I=1,N)
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
          PRINT 740, NF
  740     FORMAT (/4X,'At the return from BOBYQA',5X,
     1      'Number of function values =',I6)
          PRINT 710, F,(X(I),I=1,N)
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
          CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
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
          PRINT 300, NF,F,(W(I),I=1,N)
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
      SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,
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
      IF (IPRINT .EQ. 3) THEN
          PRINT 70, NF,F,(X(I),I=1,N)
   70      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
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



