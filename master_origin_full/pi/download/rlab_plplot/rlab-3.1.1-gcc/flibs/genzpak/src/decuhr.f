      SUBROUTINE DEADHR (NDIM,NUMFUN,A,B,MAXSUB,FUNSUB,SINGUL,ALPHA,
     +LOGF,EPSABS,EPSREL,KEY,RESTAR,NUM,LENW,WTLENG,EMAX,MINPTS,MAXPTS,
     +NSUB,RESULT,ABSERR,NEVAL,IFAIL,VALUES,ERRORS,CENTRS,HWIDTS,GREATE,
     +DIR,WORK,Q,U,E,QOLD,QNEW,UNEW,UERR,EXTERR,NE,BETA,T,DIFF,G,W,
     +RULPTS,CENTER,HWIDTH,X,SCALES,NORMS,LIST,UPOINT,ORDER)
C***BEGIN PROLOGUE DEADHR
C***KEYWORDS automatic multidimensional integrator,
C            singular integrands,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive with extrapolation.
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals I, over a hyper-rectangular
C            region hopefully satisfying for each component of I the
C            following claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C***LAST MODIFICATION 92-09-16
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions with singularities.
C            DECUHR is a driver for the integration routine
C            DEADHR, which repeatedly subdivides the region
C            of integration. DEADHR estimates the integrals and the
C            errors over the subregions and subdivides the subregion
C            with greatest estimated errors until the error request
C            is met or MAXSUB subregions are stored.
C            The subdivision is done in such a way that we at any
C            stage have only one region containing the singularity.
C            When the singular region is picked it is subdivided in
C            SINGUL + 1 new regions (cutting each of the SINGUL
C            first directions once), creating one new singular region
C            and SINGUL non-singular regions.
C            The non-singular regions are divided in two equally
C            sized parts along the direction with greatest absolute
C            fourth order difference.
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MAXSUB Integer.
C            The computations proceed until there are at most
C            MAXSUB subregions in the data structure.
C
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand in the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     SINGUL Integer.
C            Dimension of the singularity
C     ALPHA  Real
C            Degree of homogeneous function. ALPHA > -SINGUL.
C     LOGF   Integer
C            Indicating power of logarithmic term in all components.
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            (In this case the output parameters from DEADHR
C            must not be changed since the last
C            exit from DEADHR.)
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     LENW   Integer.
C            Defines the length of the working array WORK.
C            LENW should be greater or equal to
C            16*NUMFUN.
C     WTLENG Integer.
C            The number of weights in the basic integration rule.
C     EMAX
C            The maximum number of extrapolation steps.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS
C            Maximum number of function evaluations.
C     NSUB   Integer.
C            If RESTAR = 1, then NSUB must specify the number
C            of subregions stored in the previous call to DEADHR.
C
C   ON RETURN
C
C     NSUB   Integer.
C            Number of stored subregions.
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute accuracies.
C     NEVAL  Integer.
C            Number of function evaluations used by DEADHR.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXSUB or less
C              subregions processed for all values of K,
C              1 <=  K <=  NUMFUN.
C            IFAIL = 1 if MAXSUB was too small for DEADHR
C              to obtain the required accuracy. In this case DEADHR
C              returns values of RESULT with estimated absolute
C              accuracies ABSERR.
C     VALUES Real array of dimension (NUMFUN,MAXSUB+1).
C            Used to store estimated values of the integrals
C            over the subregions.
C     ERRORS Real array of dimension (NUMFUN,MAXSUB+1).
C            Used to store the corresponding estimated errors.
C     CENTRS Real array of dimension (NDIM,MAXSUB+1).
C            Used to store the centers of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,MAXSUB+1).
C            Used to store the half widths of the stored subregions.
C     GREATE Real array of dimension MAXSUB+1.
C            Used to store the greatest estimated errors in
C            all subregions.
C     DIR    Real array of dimension MAXSUB+1.
C            DIR is used to store the directions for
C            further subdivision.
C     WORK   Real array of dimension LENW.
C            Used  in DERLHR and DETRHR.
C     Q      Real array of dimension (NUMFUN,0:MAXSUB)
C            contains the estimates over the singular regions
C            (Tail-estimates).
C     U      Real array of dimension (NUMFUN,MAXSUB)
C            contains the terms in series.
C     E      Real array of dimension (NUMFUN,MAXSUB)
C            contains the estimated errors in each U-term.
C     QOLD   Real array of dimension NUMFUN.
C            The last tail corrections for all functions in the vector.
C     QNEW   Real array of dimension NUMFUN.
C            The new tail correction for all functions in the vector.
C     UNEW   Real array of dimension NUMFUN.
C            This gives the next terms in the series (new extrapolation
C            step) else it is the correction to the u-values (updating).
C     UERR   Real array of dimension NUMFUN.
C            The estimated errors of all U-terms in the series.
C     EXTERR Real array of dimension NUMFUN.
C            These errors are associated with the singular region and
C            they are the pure extrapolation errors.
C     NE     Real array of dimension (0:EMAX)
C            Dummy parameter
C     BETA   Real array of dimension (0:EMAX,0:EMAX)
C            Dummy parameter
C     T      Real array of dimension (NUMFUN,0:EMAX)
C            contains the last rows in the extrapolation tableaus.
C     DIFF   Real array of dimension NDIM.
C            Work array.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the
C            points associated with the Jth weights.
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ..., W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ..., W(I,WTLENG) , for I > 1 are null rule weights.
C     RULPTS Real array of dimension WTLENG.
C            Work array used in DEINHR.
C     CENTER Real array of dimension NDIM.
C            Work array used in DETRHR.
C     HWIDTH Real array of dimension NDIM.
C            Work array used in DETRHR.
C     X      Real array of dimension NDIM.
C            Work array used in DERLHR.
C     SCALES Real array of dimension (3,WTLENG).
C            Work array used by DEINHR and DERLHR.
C     NORMS  Real array of dimension (3,WTLENG).
C            Work array used by DEINHR and DERLHR.
C     LIST   Integer array used in DETRHR of dimension MAXSUB.
C            Is a partially sorted list, where LIST(1) is the top
C            element in a heap of subregions.
C     UPOINT Integer array of dimension MAXSUB
C            Is an array of pointers to where in the U-sequence
C            a region belongs. This information is used when updating
C            the corresponding U-term after a subdivision.
C     ORDER  Integer array of dimension NDIM.
C            Work array used to give the order in which the singular
C            region will be cut.
C***REFERENCES
C
C   T.O.Espelid and A.Genz, DECUHR: An Algorithm for Automatic
C   Integration of Singular Functions over a Hyperrectangular Region.
C   Numerical Algorithms 8(1994), PP. 201-220.
C
C   T.O.Espelid, On integrating Vertex Singularities using
C   Extrapolation,  BIT 34,1(1994) 62-79.
C
C   T.O.Espelid, On integrating Singularities using non-Uniform
C   Subdivision and Extrapolation,  in Numerical Integration IV, eds.
C   Brass and Hammerlin, Birhauser, ISNM Vol. 112(1993) 77-89.
C
C   P. van Dooren and L. de Ridder, Algorithm 6, An adaptive algorithm
C   for numerical integration over an n-dimensional cube, J.Comput.Appl.
C   Math. 2(1976)207-217.
C
C   A.C.Genz and A.A.Malik, Algorithm 019. Remarks on algorithm 006:
C   An adaptive algorithm for numerical integration over an
C   N-dimensional rectangular region, J.Comput.Appl.Math.6(1980)295-302.
C
C***ROUTINES CALLED DEINHR,DERLHR,DETRHR,DESBHR,DEXTHR.
C***END PROLOGUE DEADHR
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN,MAXSUB,KEY,LENW,RESTAR,MINPTS,MAXPTS
      INTEGER NUM,NEVAL,NSUB,IFAIL,WTLENG,SINGUL
      INTEGER UPOINT(MAXSUB),LIST(MAXSUB),LOGF,EMAX,ORDER(NDIM)
      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL,ALPHA
      DOUBLE PRECISION RESULT(NUMFUN),ABSERR(NUMFUN),DIFF(NDIM)
      DOUBLE PRECISION VALUES(NUMFUN,0:MAXSUB),ERRORS(NUMFUN,0:MAXSUB)
      DOUBLE PRECISION CENTRS(NDIM,0:MAXSUB)
      DOUBLE PRECISION HWIDTS(NDIM,0:MAXSUB),T(NUMFUN,0:EMAX)
      DOUBLE PRECISION GREATE(0:MAXSUB),DIR(0:MAXSUB)
      DOUBLE PRECISION BETA(0:EMAX,0:EMAX)
      DOUBLE PRECISION WORK(LENW),RULPTS(WTLENG),EXTERR(NUMFUN)
      DOUBLE PRECISION G(NDIM,WTLENG),W(5,WTLENG),NE(EMAX)
      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM),X(NDIM)
      DOUBLE PRECISION SCALES(3,WTLENG),NORMS(3,WTLENG),QOLD(NUMFUN)
      DOUBLE PRECISION U(NUMFUN,MAXSUB),E(NUMFUN,MAXSUB),QNEW(NUMFUN)
      DOUBLE PRECISION Q(NUMFUN,0:MAXSUB),UNEW(NUMFUN),UERR(NUMFUN)
C
C   Local variables.
C
C   SBRGNS is the number of stored subregions.
C   NDIV   The number of subregions to be divided in each main step.
C   POINTR Pointer to the position in the data structure where
C          the new subregions are to be stored.
C   TOP    is a pointer to the top element in the heap of subregions.
C   VACANT Pointer to the vacant element.
C   DIRECT Direction of subdivision.
C   ERRCOF Heuristic error coefficient defined in DEINHR and used by
C          DERLHR and DEADHR.
C   UPDATE Pointer: either to a new element or to an old one to be updat
C          Value = 0 used to indicate a new extrapolation step.
C   FIRST  Logical: .TRUE. indicates that this is the first time
C          DEXTHR is called.
C   EXSTEP The exponent sequence's step size in the error expansion.
C          Default value = 1. Can be set in a parameter statement to a
C          different integer > 1 if that is appropriate. Advantage:
C          avoid to eliminate terms that are not present for the given
C          function. Warning: If the step size is set wrong (too big)
C          then this may result in poor performance.
C
      LOGICAL FIRST
      INTEGER I,J,K,SBRGNS,TOP,UPDATE,EXSTEP,NUMU
      INTEGER NDIV,POINTR,DIRECT,INDEX,VACANT
      DOUBLE PRECISION ERRCOF(6)
      PARAMETER (EXSTEP=1)
C
C***FIRST EXECUTABLE STATEMENT DEADHR
C
C   Call DEINHR to compute the weights and abscissas of
C   the function evaluation points.
C
      CALL DEINHR (NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C
C   If RESTAR = 1, then this is a restart run.
C
      IF (RESTAR.EQ.1) THEN
         SBRGNS=NSUB
         GO TO 50
      END IF
C
C   Initialize the SBRGNS, CENTRS,  HWIDTS and ORDER
C   Note the singular region will stay in position 0 and will not be
C   a part of the heap at any time.
C
      SBRGNS=0
      DO 10 J=1,NDIM
         CENTRS(J,0)=(A(J)+B(J))/2
         HWIDTS(J,0)=ABS(B(J)-A(J))/2
         ORDER(J)=J
 10   CONTINUE
C
C   Apply DERLHR over the whole region.
C   DERLHR will only change DIR(0) if the input value is non-negative.
C
      DIR(0)=-SINGUL
      CALL DERLHR (NDIM,CENTRS(1,0),HWIDTS(1,0),WTLENG,G,W,ERRCOF,
     +     NUMFUN,FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,0),ERRORS(1,0),
     +     DIR(0),GREATE(0),DIFF,ORDER)
C
C   Initialize RESULT, ABSERR, Q(*,0), NEVAL, NUMU and FIRST.
C
      DO 20 J=1,NUMFUN
         RESULT(J)=VALUES(J,0)
         ABSERR(J)=ERRORS(J,0)
         Q(J,0)=VALUES(J,0)
 20   CONTINUE
      NEVAL=NUM
      FIRST=.TRUE.
      NUMU=0
C
C   Initialize E and U
C
      DO 40 J=1,NUMFUN
         DO 30 I=1,MAXSUB
            E(J,I)=0
            U(J,I)=0
 30      CONTINUE
 40   CONTINUE
C
C   Initialize the pointer LIST(1) to point on the singular region.
C
      LIST(1)=0
C
C***End initialization.
C
C***Begin loop while the error is too large,
C
C     First we determine if we will subdivide the singular region
C     or a regular region, and then the number, NDIV, of new subregions.
C
 50   TOP=LIST(1)
C
C     First we determine if we will subdivide the singular region or a
C     non-singular region, and then the number, NDIV, of new subregions.
C
      IF (GREATE(TOP).LT.(GREATE(0))) THEN
         VACANT=0
      ELSE
         VACANT=TOP
      END IF
      DIRECT=DIR(VACANT)
      NDIV=MAX(2,-DIRECT+1)
C
C     Check if NEVAL+NDIV*NUM is less than or equal to MAXPTS:
C     MAXPTS is the maximum number of function evaluations that are
C     allowed to be computed.
C
      IF (NEVAL+NDIV*NUM.LE.MAXPTS) THEN
C
C     We are allowed to divide further: prepare to remove the region
C     with greatest error estimate from the heap or replace the singular
C     region with a smaller one.
C
C     Let POINTR point to the first free position in the data structure.
C
         POINTR=SBRGNS+1
C
C     Adjust, if necessary, the heap. (Reduce the size by one element.
C     Note: the information about the region we will subdivide
C     remains in the vacant position.)
C
         IF (VACANT.GT.0) THEN
            CALL DETRHR (1,SBRGNS,GREATE,LIST,K)
         END IF
C
C
C     Determine the new subregions.
C
         CALL DESBHR (NDIM,CENTRS,HWIDTS,DIR,DIRECT,POINTR,VACANT,ORDER)
C
C     Determine if this is a new extrapolation step or an update.
C     UPDATE will point to the element in the
C     U-series to be corrected or created.
C
         IF (VACANT.EQ.0) THEN
            NUMU=NUMU+1
            UPDATE=NUMU
         ELSE
            UPDATE=UPOINT(VACANT)
         END IF
C
C     Apply the basic rule to the new regions and let each new region
C     (except the singular region) point to its U-element.
C     First the region in the VACANT position. This region is the only
C     one that may be the singular and thus give a new Q-element.
C
         INDEX=VACANT
         IF (VACANT.EQ.0) THEN
            DO 60 J=1,NUMFUN
               UERR(J)=0
               UNEW(J)=0
 60         CONTINUE
         ELSE
            DO 70 J=1,NUMFUN
               UERR(J)=-ERRORS(J,INDEX)
               UNEW(J)=-VALUES(J,INDEX)
 70         CONTINUE
         END IF
         CALL DERLHR (NDIM,CENTRS(1,INDEX),HWIDTS(1,INDEX),WTLENG,G,W,
     +    ERRCOF,NUMFUN,FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,INDEX),
     +    ERRORS(1,INDEX),DIR(INDEX),GREATE(INDEX),DIFF,ORDER)
         IF (VACANT.EQ.0) THEN
            DO 80 J=1,NUMFUN
               Q(J,NUMU)=VALUES(J,0)
 80         CONTINUE
         ELSE
            UPOINT(INDEX)=UPDATE
            DO 90 J=1,NUMFUN
               UERR(J)=UERR(J)+ERRORS(J,INDEX)
               UNEW(J)=UNEW(J)+VALUES(J,INDEX)
 90         CONTINUE
         END IF
C
C     Then the rest of the regions
C
         DO 110 I=2,NDIV
            INDEX=POINTR+I-2
            CALL DERLHR (NDIM,CENTRS(1,INDEX),HWIDTS(1,INDEX),WTLENG,G,
     +       W,ERRCOF,NUMFUN,FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,INDEX),
     +       ERRORS(1,INDEX),DIR(INDEX),GREATE(INDEX),DIFF,ORDER)
            UPOINT(INDEX)=UPDATE
            DO 100 J=1,NUMFUN
               UERR(J)=UERR(J)+ERRORS(J,INDEX)
               UNEW(J)=UNEW(J)+VALUES(J,INDEX)
 100        CONTINUE
 110     CONTINUE
C
C     Compute the E and U terms (These may be new terms or terms that
C     have to be updated), QNEW and QOLD have no influence in a true
C     update.
C
         DO 120 J=1,NUMFUN
            QNEW(J)=Q(J,NUMU)
            QOLD(J)=Q(J,NUMU-1)
            U(J,UPDATE)=U(J,UPDATE)+UNEW(J)
            E(J,UPDATE)=E(J,UPDATE)+UERR(J)
 120     CONTINUE
C
C     Do the extrapolation and compute the global results and errors.
C     UPDATE is used to signal an extrapolation step.
C
         IF (VACANT.EQ.0) THEN
            UPDATE=0
         END IF
         CALL DEXTHR (NUMFUN,ALPHA,LOGF,SINGUL,EXSTEP,NUMU,T,UPDATE,
     +    UNEW,QNEW,QOLD,FIRST,EMAX,E,RESULT,ABSERR,EXTERR,NE,BETA)
         FIRST=.FALSE.
         NEVAL=NEVAL+NDIV*NUM
C
C     Change the error estimates in position 0.
C     We define the pure extrapolation error to be the new error of
C     the singular region.
C
         GREATE(0)=0
         DO 130 J=1,NUMFUN
            GREATE(0)=MAX(GREATE(0),EXTERR(J))
            ERRORS(J,0)=EXTERR(J)
 130     CONTINUE
C
C     Store results in heap.
C
         IF (VACANT.GT.0) THEN
            CALL DETRHR (2,POINTR-1,GREATE,LIST,VACANT)
         END IF
         DO 140 I=2,NDIV
            INDEX=POINTR+I-2
            CALL DETRHR (2,INDEX,GREATE,LIST,INDEX)
 140     CONTINUE
         SBRGNS=POINTR+NDIV-2
C
C     Check for termination.
C
         IF (NEVAL.LT.MINPTS) THEN
            GO TO 50
         END IF
         DO 150 J=1,NUMFUN
            IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)).AND.ABSERR(J).GT.
     +       EPSABS) THEN
               GO TO 50
            END IF
 150     CONTINUE
         IFAIL=0
C
C     Else we did not succeed with the
C     given value of MAXSUB.
C
      ELSE
         IFAIL=1
      END IF
      NSUB=SBRGNS
C
C***END DEADHR
C
      END
      SUBROUTINE DECALP ( NDIM, A, B, NUMFUN, ALPHA, SINGUL, LOGF,
     +     FUNSUB, X, FUNS )
C***BEGIN PROLOGUE DECALP
C***PURPOSE  DECALP computes  ALPHA and LOGF.
C***LAST MODIFICATION 94-05-04
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     NUMFUN Integer.
C            Number of components of the integral.
C     ALPHA  Real.
C            Degree of homogeneous function.
C     SINGUL Integer
C            Dimension of the singularity.
C     LOGF   Integer.
C            LOGF = 0, then no logarithmic singularity
C            LOGF = 1, then a logarithmic singularity of order 1.
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C     X       Real work array of length at least NDIM.
C     FUNS    Real work array of length at least NUMFUN.
C
C   ON RETURN
C
C     ALPHA  Real.
C            Degree of homogeneous function. The value will be estimated
C            by this routine.
C      LOGF  Integer.
C            LOGF = 0, then no logarithmic singularity
C            LOGF = 1, then a logarithmic singularity of order 1.
C            The value will be set by this routine.
C
C***END PROLOGUE DECALP
C
C   Global variables.
C
      EXTERNAL FUNSUB
      DOUBLE PRECISION A(*), B(*), X(*), FUNS(*), ALPHA
      INTEGER NDIM, NUMFUN, LOGF, SINGUL
C
C   Local variables.
C

      INTEGER I, N, J, ITW, BDG, L, NB, INDEX, SIGNAL
      DOUBLE PRECISION H, LGTWO, TERM1, TERM2, KEST
      PARAMETER ( BDG = 10, L = 10 )
      DOUBLE PRECISION T(0:2*BDG+L), SUMX, SUMY
C
C***FIRST EXECUTABLE STATEMENT DECALP
C

      LGTWO = LOG(2D0)
      H = 2

C
C   Estimate ALPHA and LOGF
C
C   Step 1: Create the ratio table
C
C
C   Compute points along a line by halving the singular directions.
C
 10   H = H/2
      DO 20 I = 1,NDIM
         X(I) = A(I) + H*(B(I)-A(I))
 20   CONTINUE
      CALL FUNSUB ( NDIM, X, NUMFUN, FUNS )
      SUMY = 0
      DO 30 I = 1,NUMFUN
         SUMY = SUMY + ABS( FUNS(I) )
 30   CONTINUE
      IF ( SUMY .LE. 0 ) GO TO 10
      DO 60 J = 0,2*BDG+L
         SUMX = SUMY
         H = H/2
         DO 40 I = 1,SINGUL
            X(I) = A(I) + H*(B(I)-A(I))
 40      CONTINUE
         CALL FUNSUB ( NDIM, X, NUMFUN, FUNS )
         SUMY = 0
         DO 50 I = 1,NUMFUN
            SUMY = SUMY + ABS( FUNS(I) )
 50      CONTINUE
C
C     Check if we can compute the next ratio
C
         IF ( SUMY .LE. 0 ) GO TO 10
         T(J) = LOG( SUMX/SUMY )/LGTWO
 60   CONTINUE
C
C    Step 2: Perform linear extrapolation
C
      DO 80 I = 1,2*BDG+L
         N = 0
         DO 70 J = I-1,MAX(0,I-L),-1
            N = 2*N+1
            T(J) = T(J+1) + (T(J+1) - T(J))/N
 70      CONTINUE
 80   CONTINUE
C
C     Now T(0), T(1), ..., T(2*BDG) are all approximations to ALPHA.
C
C     Check for logarithmic term:
C
      IF ( ABS( T(2*BDG)-T(2*BDG-1) ) .GT. 5D-6 ) THEN
         LOGF = 1
C
C     Step 3: Eliminate the effect of the logarithmic term on the
C     alpha estimates:
C
C     SIGNAL is the index at which we have to stop the
C     extrapolation. This is based on the paper by Bjorstad,
C     Grosse and Dahlquist, BIT, 21, 1981. Based on our experience
C     we allow one step more than the stopping criteria indicates.
C
         INDEX = 0
         SIGNAL = -1
         NB = 2*BDG
 90      NB = NB - 2
         ITW = 2*BDG - NB
         IF ( NB .GT. 0) THEN
            IF ( ABS( T(NB+1)+T(NB-1)-2*T(NB) ) .GT. 0 ) THEN
               TERM1 = ( T(NB+1)-T(NB) )/( T(NB+1)+T(NB-1)-2*T(NB) )
            ELSE
               TERM1 = 0
            END IF
            IF ( ABS( T(NB+2)+T(NB)-2*T(NB+1) ) .GT. 0) THEN
               TERM2 = ( T(NB+2)-T(NB+1) )/( T(NB+2)+T(NB)-2*T(NB+1) )
            ELSE
               TERM2 = 0
            END IF
            KEST = -1 - 1/( TERM2-TERM1 )
         ELSE
            KEST = 0
         END IF
         DO  100 J = 0,NB
            IF ( ABS( T(J+2)+T(J)-2*T(J+1) ) .GT. 0 ) THEN
               T(J) = T(J+1) - ITW*( T(J+2)-T(J+1) )*( T(J+1)-T(J) )
     +              /( ( T(J+2)+T(J)-2*T(J+1) )*( ITW-1 ) )
            ELSE
               T(J) = T(J+1)
            END IF
 100     CONTINUE
         IF (ABS(KEST-ITW+1).GT.1 .AND. SIGNAL.LT.0) SIGNAL=MAX(0,NB-2)
         IF ( NB .NE. SIGNAL .AND. NB .GT. 0 ) GO TO 90
         INDEX = SIGNAL
         ALPHA = T(INDEX)
      ELSE
         LOGF = 0
         ALPHA = T(L)
      END IF
C
C***END DECALP
C
      END
      SUBROUTINE DECHHR (MAXDIM,NDIM,NUMFUN,A,B,ALPHA,SINGUL,LOGF,
     +     MINPTS,MAXPTS,EPSABS,EPSREL,KEY,NW,RESTAR,EMAX,WRKSUB,
     +     NUM,MAXSUB,KEYF,IFAIL,WTLENG)
C***BEGIN PROLOGUE DECHHR
C***PURPOSE  DECHHR checks the validity of input parameters to DECUHR.
C
C***LAST MODIFICATION 94-08-02
C***DESCRIPTION
C            DECHHR computes NUM, MAXSUB, KEYF, WTLENG and
C            IFAIL as functions of the input parameters to DECUHR,
C            and checks the validity the input parameters to DECUHR.
C
C   ON ENTRY
C
C     MAXDIM Integer.
C            The maximum allowed number of dimensions.
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     ALPHA  Real.
C            Degree of homogeneous function.
C     SINGUL Integer.
C            Dimension of the singularity.
C     LOGF   Integer.
C            LOGF = 0, then no logarithmic singularity
C            LOGF = 1, then a logarithmic singularity of order 1.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
C              NUM = 65
C            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
C              NUM = 127
C            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
C              NUM = 1 + 2*NDIM + 6*NDIM*NDIM +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            Elseif (KEY = 4) Then
C              NUM = 1 + 2*NDIM*(NDIM+2) + 2**NDIM
C            MAXPTS >= (SINGUL+2)*NUM and MAXPTS >= MINPTS
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require the most adaptivity.
C     NW     Integer.
C            Defines the length of the working array WORK.
C            Let WRKSUB denote the maximum allowed number of subregions
C            NW should then be greater or equal to
C                3+17*NUMFUN+WRKSUB*(NDIM+NUMFUN+1)*2+EMAX+
C                NUMFUN*(3*WRKSUB+9+EMAX)+(EMAX+1)**2 +3*NDIM
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C     EMAX   Integer
C            The maximum number of extrapolation steps.
C     WRKSUB Integer.
C            Maximum size of MAXSUB.
C
C   ON RETURN
C
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     MAXSUB Integer.
C            The maximum possible number of subregions for the
C            given values of MAXPTS, KEY and NDIM. (See the code)
C            MAXSUB is not allowed to be greater than WRKSUB.
C            The fact that MAXSUB may be smaller allows a restart
C            option with MAXPTS increased.
C     KEYF   Integer.
C            Key to selected integration rule.
C     IFAIL  Integer.
C            IFAIL =  0 for normal exit.
C            IFAIL =  2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL =  3 if NDIM is less than 2 or NDIM greater than
C                       MAXDIM.
C            IFAIL =  4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL =  5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL =  6 if NUMFUN less than 1.
C            IFAIL =  7 if A(I) >= B(I), for a value of 0 < I < NDIM+1.
C            IFAIL =  8 if MAXPTS is less than MINPTS.
C            IFAIL =  9 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 10 if RESTAR > 1 or RESTAR < 0.
C            IFAIL = 11 if SINGUL > NDIM or SINGUL < 1.
C            IFAIL = 12 if LOG < 0 or LOG > 1.
C            IFAIL = 13 if MAXPTS is less than (SINGUL+2)*NUM.
C            IFAIL = 14 if ALPHA <=  -SINGUL.(Set outside this routine).
C            IFAIL = 15 if EMAX < 1.
C            IFAIL = 16 if NW is too small.
C            IFAIL = 17 if WRKSUB is too small.
C     WTLENG Integer.
C            The number of generators of the chosen integration rule.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE DECHHR
C
C   Global variables.
C
      INTEGER NDIM,NUMFUN,SINGUL,MINPTS,MAXPTS,KEY,NW,MAXSUB,LOGF
      INTEGER RESTAR,NUM,KEYF,IFAIL,MAXDIM,WTLENG,EMAX,WRKSUB
      DOUBLE PRECISION A(NDIM),B(NDIM),ALPHA,EPSABS,EPSREL
C
C   Local variables.
C
      INTEGER LIMIT,J,C1,CSING
C
C***FIRST EXECUTABLE STATEMENT DECHHR
C
      IFAIL=0
C
C   Check KEY.
C
      IF (KEY.LT.0.OR.KEY.GT.4) THEN
         IFAIL=2
         RETURN
      END IF
C
C   Check NDIM.
C
      IF (NDIM.LT.2.OR.NDIM.GT.MAXDIM) THEN
         IFAIL=3
         RETURN
      END IF
C
C   For KEY = 1, NDIM must be equal to 2.
C
      IF (KEY.EQ.1.AND.NDIM.NE.2) THEN
         IFAIL=4
         RETURN
      END IF
C
C   For KEY = 2, NDIM must be equal to 3.
C
      IF (KEY.EQ.2.AND.NDIM.NE.3) THEN
         IFAIL=5
         RETURN
      END IF
C
C   For KEY = 0, we point at the selected integration rule.
C
      IF (KEY.EQ.0) THEN
         IF (NDIM.EQ.2) THEN
            KEYF=1
         ELSE IF (NDIM.EQ.3) THEN
            KEYF=2
         ELSE
            KEYF=3
         END IF
      ELSE
         KEYF = KEY
      END IF
C
C   Compute NUM and WTLENG as a function of KEYF and NDIM.
C
      IF (KEYF.EQ.1) THEN
         NUM = 65
         WTLENG = 14
      ELSE IF (KEYF.EQ.2) THEN
         NUM = 127
         WTLENG = 13
      ELSE IF (KEYF.EQ.3) THEN
         NUM = 1 + 2*NDIM + 6*NDIM*NDIM
     +        + 4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
         WTLENG = 9
         IF (NDIM.EQ.2) WTLENG = 8
      ELSE IF (KEYF.EQ.4) THEN
         NUM = 1 + 2*NDIM*(NDIM+2) + 2**NDIM
         WTLENG = 6
      END IF
C
C   Check NUMFUN.
C
      IF (NUMFUN.LT.1) THEN
         IFAIL=6
         RETURN
      END IF
C
C   Check upper and lower limits of integration.
C
      DO 10 J=1,NDIM
         IF (A(J).GE.B(J)) THEN
            IFAIL=7
            RETURN
         END IF
 10   CONTINUE
C
C   Check MAXPTS >= MINPTS.
C
      IF (MAXPTS.LT.MINPTS) THEN
         IFAIL=8
         RETURN
      END IF
C
C   Check accuracy requests.
C
      IF (EPSABS.LT.0.AND.EPSREL.LT.0) THEN
         IFAIL=9
         RETURN
      END IF
C
C    Check RESTAR.
C
      IF (RESTAR.NE.0.AND.RESTAR.NE.1) THEN
         IFAIL=10
         RETURN
      END IF
C
C    Check SINGUL.
C
      IF (SINGUL.LT.1.OR.SINGUL.GT.NDIM) THEN
         IFAIL=11
         RETURN
      END IF
C
C    Check LOGF.
C
      IF (LOGF.LT.0.OR.LOGF.GT.1) THEN
         IFAIL=12
         RETURN
      END IF
C
C   Check that MAXPTS allows at least one subdivision of the
C   singular region in SINGUL + 1 new pieces.
C
      IF (MAXPTS.LT.(SINGUL+2)*NUM) THEN
         IFAIL=13
         RETURN
      END IF
C
C    Check EMAX.
C
      IF (EMAX.LT.1) THEN
         IFAIL=15
         RETURN
      END IF
C
C   Compute MAXSUB. This is a worst case computation. We assume that
C   the singular region will be cut CSING times, each time in
C   SINGUL + 1 new subregions. CSING is the maximum number of times
C   that this can happen, given MAXPTS and NUM. Then we assume that only
C   bisections (of non-singular regions) will take place. C1 is the
C   maximum number of times this can happen, given MAXPTS, NUM and CSING
C
      CSING = (MAXPTS-NUM)/((SINGUL+1)*NUM)
      C1 = (MAXPTS-NUM*(1+(SINGUL+1)*CSING))/(2*NUM)
      IF ( C1.EQ.0 ) THEN
         C1 = -1
      END IF
      MAXSUB = 2 + CSING*SINGUL + C1
C
C   Check size of double precision workspace.
C
      LIMIT= 3 + 17*NUMFUN + WRKSUB*(NDIM+NUMFUN+1)*2 + EMAX
     +     + NUMFUN*(3*WRKSUB+9+EMAX) + (EMAX+1)**2 + 3*NDIM
      IF ( NW.LT.LIMIT ) THEN
         IFAIL = 16
         RETURN
      END IF
C
C   Check if MAXSUB can be achieved.
C
      IF (MAXSUB.GT.WRKSUB) THEN
         IFAIL=17
         RETURN
      END IF
C
C***END DECHHR
C
      END
      SUBROUTINE DECUHR ( NDIM, NUMFUN, A, B, MINPTS, MAXPTS, FUNSUB,
     +     SINGUL, ALPHA, LOGF, EPSABS, EPSREL, KEY, WRKSUB, NW, RESTAR,
     +     EMAX, RESULT, ABSERR, NEVAL, IFAIL, WORK, IWORK)
C***BEGIN PROLOGUE DECUHR
C***DATE WRITTEN   930816   (YYMMDD)
C***REVISION DATE  940802   (YYMMDD)
C***CATEGORY NO. H2B1A1
C***AUTHOR
C            Terje O. Espelid, Department of Informatics,
C            University of Bergen,  Hoyteknologisenteret,
C            N-5020 Bergen, Norway
C            Email..  terje@ii.uib.no
C
C            Alan Genz,
C            Department of Mathematics
C            Washington State University
C            Pullman, WA 99164-3113, USA
C            Email..  genz@gauss.math.wsu.edu
C
C***KEYWORDS automatic multidimensional integrator,
C            singular integrands,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive with extrapolation.
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals
C
C      B(1) B(2)     B(NDIM)
C     I    I    ... I       (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
C      A(1) A(2)     A(NDIM)  1  2      NUMFUN
C
C       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN,
C              I   I  1  2      NDIM
C
C            hopefully satisfying for each component of I the following
C            claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C            The vector of definite integrals are all assumed to have
C            a singularity of dimension SINGUL at the point,
C            X(1) = A(1), X(2) = A(2), ... , X(SINGUL) = A(SINGUL).
C            The singularity is assumed to be caused by a homogeneous
C            function of degree ALPHA at the point. A logarithmic
C            singularity at the same point can also be handled.
C
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions with singularities.
C            DECUHR is a driver for the integration routine
C            DEADHR, which repeatedly subdivides the region
C            of integration. DEADHR estimates the integrals and the
C            errors over the subregions and subdivides the subregion
C            with greatest estimated errors until the error request
C            is met or MAXPTS function evaluations have been used.
C            The subdivision is done in such a way that we at any
C            stage have only one region containing the singularity.
C            When the singular region is picked it is subdivided in
C            SINGUL + 1 new regions (cutting each of the SINGUL
C            first directions once), creating one new singular region
C            and SINGUL non-singular regions. When a non-singular
C            region is picked we divide this in two halves.
C
C            For NDIM = 2 the default integration rule is of
C            degree 13 and uses 65 evaluation points.
C            For NDIM = 3 the default integration rule is of
C            degree 11 and uses 127 evaluation points.
C            For NDIM greater then 3 the default integration rule
C            is of degree 9 and uses NUM evaluation points where
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            The degree 9 rule may also be applied for NDIM = 2
C            and NDIM = 3.
C            A rule of degree 7 is available in all dimensions.
C            The number of evaluation
C            points used by the degree 7 rule is
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C
C            When DECUHR computes estimates to a vector of
C            integrals, all components of the vector are given
C            the same treatment. That is, I(F ) and I(F ) for
C                                            J         K
C            J not equal to K, are estimated with the same
C            subdivision of the region of integration.
C            For integrals with enough similarity, we may save
C            time by applying DECUHR to all integrands in one call.
C            For integrals that varies continuously as functions of
C            some parameter, the estimates produced by DECUHR will
C            also vary continuously when the same subdivision is
C            applied to all components. This will generally not be
C            the case when the different components are given
C            separate treatment. It is essential that all components
C            have the same type of singular behavior. Thus only one
C            value for ALPHA is allowed, which thus has to be the same
C            for all components.
C
C            On the other hand this feature should be used with
C            caution when the different components of the integrals
C            require clearly different subdivisions.
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <=  15.
C     NUMFUN Integer.
C            Number of components of the integral.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
C              NUM = 65
C            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
C              NUM = 127
C            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
C              NUM = 1 + 2*NDIM + 6*NDIM*NDIM +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            Elseif (KEY = 4) Then
C              NUM = 1 + 2*NDIM*(NDIM+2) + 2**NDIM
C            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
C            For 3 < NDIM < 13 the minimum values for MAXPTS are:
C             NDIM =    4   5   6    7    8    9    10   11    12
C            KEY = 3:  459 819 1359 2151 3315 5067 7815 12351 20235
C            KEY = 4:  195 309  483  765 1251 2133 3795  7005 13299
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C     SINGUL Integer.
C            Dimension of the singularity
C     ALPHA  Real.
C            Degree of homogeneous function.
C            If input ALPHA <= -SINGUL, then ALPHA is estimated by the
C            code and this estimate is used for extrapolation.
C            SINGUL is assumed to be given correct.
C     LOGF   Integer.
C            LOGF = 0, then no logarithmic singularity
C            LOGF = 1, then a logarithmic singularity of order 1.
C            It's value will be estimated by the code if
C            input ALPHA <= -SINGUL.
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     WRKSUB Integer.
C            The maximum allowed number of subregions.
C     NW     Integer.
C            Defines the length of the working array WORK.
C            NW should be >=
C              3 + 17*NUMFUN + WRKSUB*(NDIM+NUMFUN+1)*2 + EMAX
C              + NUMFUN*(3*WRKSUB+9+EMAX) + (EMAX+1)**2 + 3*NDIM
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            In this case the only parameters for DECUHR that may
C            be changed (with respect to the previous call of DECUHR)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C       EMAX Integer.
C            The maximum number of extrapolation steps.
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute errors.
C     NEVAL  Integer.
C            Number of function evaluations used by DECUHR.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C            IFAIL = 1 if MAXPTS was too small for DECUHR
C              to obtain the required accuracy. In this case DECUHR
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C            IFAIL =  2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL =  3 if NDIM is less than 2 or NDIM greater than 15.
C            IFAIL =  4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL =  5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL =  6 if NUMFUN is less than 1.
C            IFAIL =  7 if A(J) >= B(J), for some J.
C            IFAIL =  8 if MAXPTS is less than MINPTS.
C            IFAIL =  9 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 10 if RESTAR > 1 or RESTAR < 0.
C            IFAIL = 11 if SINGUL > NDIM or SINGUL < 1.
C            IFAIL = 12 if LOG < 0 or LOG > 1.
C            IFAIL = 13 if MAXPTS is less than (SINGUL+2)*NUM.
C            IFAIL = 14 if ALPHA <=  -SINGUL. (Integral does not exist).
C            IFAIL = 15 if EMAX < 1.
C            IFAIL = 16 if NW is too small.
C            IFAIL = 17 if WRKSUB is too small.
C     WORK   Real array of dimension NW.
C            Used as working storage.
C            WORK(NW) = NSUB, the number of subregions in the data
C            structure.
C            WORK(1),...,WORK(NUMFUN*(WRKSUB+1)) contain
C              the estimated components of the integrals over the
C              subregions.
C            WORK(NUMFUN*(WRKSUB+1)+1),...,WORK(2*NUMFUN*(WRKSUB+1))
C              contain the estimated errors over the subregions.
C            WORK(2*NUMFUN*(WRKSUB+1)+1),...,WORK(2*NUMFUN*(WRKSUB+1)+
C              NDIM*(WRKSUB+1)) contain the centers of the subregions.
C            WORK((2*NUMFUN*+NDIM)*(WRKSUB+1)+1),...,
C              WORK((2*NUMFUN+2*NDIM)*(WRKSUB+1)) contain the
C              subregion half widths.
C            WORK((2*NUMFUN+2*NDIM)*(WRKSUB+1)+1),...,
C              WORK((2*NUMFUN+2*NDIM+1)*(WRKSUB+1))
C              contain the greatest errors in each subregion.
C            WORK((2*NUMFUN+2*NDIM+1)*(WRKSUB+1)+1),...,
C              WORK(2*(NUMFUN+NDIM+1)*(WRKSUB+1)) contain the
C              direction of subdivision in each subregion.
C            WORK(2*(NDIM+NUMFUN+1)*(WRKSUB+1)+1),...,
C              WORK(2*(NDIM+NUMFUN+1)*(WRKSUB+1)+LENW) is used as
C              temporary storage in DEADHR.
C            WORK(2*(NDIM+NUMFUN+1)*(WRKSUB+1)+LENW+1),...,
C             WORK((2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW) contain the
C             tail estimates, Q, over the singular regions.
C            WORK((2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW+1),...,
C             WORK((2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW+NUMFUN*WRKSUB)
C             contain the U-sequence.
C            WORK((2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW+NUMFUN*WRKSUB+1),
C             ...,
C             WORK((2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW+2*NUMFUN*WRKSUB),
C             contain the error estimates to the elements in the U-seq.
C            Define
C            BASE=(2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW+2*NUMFUN*WRKSUB
C            then we have:
C            WORK(BASE+1),...,WORK(BASE+NUMFUN), contain old Q
C             estimates.
C            WORK(BASE+NUMFUN+1),...,WORK(BASE+2*NUMFUN), contain
C             new Q estimates.
C            WORK(BASE+2*NUMFUN+1), ...,WORK(BASE+3*NUMFUN), contain
C             new U-estimates.
C            WORK(BASE+3*NUMFUN+1),...,WORK(BASE+4*NUMFUN), contain
C             error estimates of the U-sequence.
C            WORK(BASE+4*NUMFUN+1),...,WORK(BASE+5*NUMFUN), contain
C             estimates of the extrapolation error.
C            WORK(BASE+5*NUMFUN+1),...,WORK(BASE+5*NUMFUN+1+EMAX),
C             contain extrapolation denominators.
C            WORK(BASE+5*NUMFUN+1+EMAX+1),...,
C             WORK(BASE+5*NUMFUN+(1+EMAX)**2), contain a work array
C             used to compute the global array.
C            WORK(BASE+5*NUMFUN+(1+EMAX)**2+1),...,
C             WORK(BASE+5*NUMFUN+(1+EMAX)**2+NUMFUN*(EMAX+1)), contain
C             extrapolation tableau.
C            WORK(BASE+5*NUMFUN+(1+EMAX)**2+NUMFUN*(EMAX+1)+1),...,
C             WORK(BASE+5*NUMFUN+(1+EMAX)**2+NUMFUN*(EMAX+1)+NDIM),
C             contain a working array used in the rule evaluation
C           routine.
C     IWORK  Integer array of dimension 2*WRKSUB + NDIM.
C            Used as working storage.
C
C***LONG DESCRIPTION
C
C   The information for each subregion is contained in the
C   data structure WORK.
C   When passed on to DEADHR, WORK is split into nineteen arrays:
C   VALUES, ERRORS, CENTRS, HWIDTS, GREATE, DIR,
C   WORK, Q,U,E,QOLD,QNEW,UNEW,UERR,EXTERR,NE,BETA,T and DIFF.
C   The first 6 are connected to each subregion:
C    VALUES contains the estimated values of the integrals.
C    ERRORS contains the estimated errors.
C    CENTRS contains the centers of the subregions.
C    HWIDTS contains the half widths of the subregions.
C    GREATE contains the greatest estimated error for each subregion.
C    DIR    contains the directions for further subdivision.
C   Number 7 is a work array for the integration rules:
C    WORK is used as work array in DEADHR.
C   The next 11 are used by the extrapolation procedure:
C    Q contains the estimates over the singular regions(Tail-estimates).
C    U contains the terms in the series.
C    E contains the estimated errors in each U-term.
C    QOLD  contains the estimates for the previous singular region.
C    QNEW  contains the estimates for the current singular region.
C    UNEW  contains the estimates for the new U-term (or update
C          corrections).
C    UERR  contains the estimates for the error in the new U-term.
C    EXERR contains the estimates for the extrapolation error.
C    NE    work array used in the extrapolation process.
C    BETA  work array used to compute the global error.
C    T     contains the last row in the extrapolation tableau.
C   And finally a work array for the singular region.
C    DIFF  work array used in the rule evaluation routine.
C
C   The integer work array is spit in three: LIST, UPOINT and ORDER.
C    LIST   contains the pointers used in maintaining the heap.
C    UPOINT contains pointers from each subregion to the U-term
C           where this region is one part. This is used in the
C           the updating process.
C    ORDER  Pointer array giving the subdivision sequence of the
C           singular region. Connected to DIFF.
C   The data structures for the subregions are in DEADHR organized
C   as a heap, and the size of GREATE(I) defines the position of
C   region I in the heap. The heap is maintained by the program
C   DETRHR. The singular region is not part of the heap due to the
C   fact that the error estimates associated with this region will
C   have to be updated even when this region has not been subdivided.
C
C   The subroutine for estimating the integral and the error over
C   each subregion, DERLHR, uses WORK2 as a work array.
C
C***REFERENCES
C
C   T.O.Espelid and A.Genz, DECUHR: An Algorithm for Automatic
C   Integration of Singular Functions over a Hyperrectangular Region.
C   Numerical Algorithms 8(1994), PP. 201-220.
C
C   T.O.Espelid, On integrating Vertex Singularities using
C   Extrapolation,  BIT, 34:62-79, 1994.
C
C   T.O.Espelid, On integrating Singularities using non-Uniform
C   Subdivision and Extrapolation,  In Numerical Integration IV,
C   H. Brass and G. Hammerlin (Eds), Birkhauser Verlag, Basel,
C   Vol. 112:77-89, 1993.
C
C   J.Berntsen, T.O.Espelid and A.Genz, An Adaptive Algorithm
C   for the Approximate Calculation of Multiple Integrals,
C   ACM Transactions on Mathematical Software,Vol.17,No.4,
C   December 1991,Pages 437-451.
C
C   J.Berntsen, T.O.Espelid and A.Genz, DCUHRE: An Adaptive
C   Multidimensional Integration Routine for a Vector of
C   Integrals, ACM Transactions on Mathematical Software, Vol.17,No.4,
C   December 1991,Pages 452-456.
C
C***ROUTINES CALLED DECHHR, DEADHR, DECALP
C***END PROLOGUE DECUHR
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN,MINPTS,MAXPTS,KEY,NW,RESTAR,WRKSUB
      INTEGER NEVAL,IFAIL,SINGUL,IWORK(2*WRKSUB+NDIM),EMAX,LOGF
      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL,ALPHA
      DOUBLE PRECISION RESULT(NUMFUN),ABSERR(NUMFUN),WORK(NW)
C
C   Local variables.
C
C   MAXDIM Integer.
C          The maximum allowed value of NDIM.
C   MAXWT  Integer. The maximum number of weights used by the
C          integration rule.
C   WTLENG Integer.
C          The number of generators used by the selected rule.
C   WORK2  Real work space. The length
C          depends on the parameters MAXDIM and MAXWT.
C   MAXSUB Integer.
C          The maximum allowed number of subdivisions
C          for the given values of KEY, NDIM and MAXPTS.
C   NUM    Integer. The number of integrand evaluations needed
C          over each subregion.
C
      INTEGER MAXWT, WTLENG, MAXDIM, LENW2, MAXSUB
      INTEGER NUM, NSUB, LENW, KEYF
      PARAMETER ( MAXDIM = 15 )
      PARAMETER ( MAXWT = 14 )
      PARAMETER ( LENW2 = 2*MAXDIM*(MAXWT+1) + 12*MAXWT + 2*MAXDIM )
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I9, I10,
     +     I11, I12, I13, I14, I15, I16, I17, I18, I19
      INTEGER K1, K2, K3, K4, K5, K6, K7, K8
      DOUBLE PRECISION WORK2(LENW2)
C
C***FIRST EXECUTABLE STATEMENT DECUHR
C
C   Compute NUM, WTLENG, MAXSUB,
C   and check the input parameters.
C
      CALL DECHHR ( MAXDIM, NDIM, NUMFUN, A, B, ALPHA, SINGUL, LOGF,
     +     MINPTS, MAXPTS, EPSABS, EPSREL, KEY, NW, RESTAR, EMAX,
     +     WRKSUB, NUM, MAXSUB, KEYF, IFAIL, WTLENG )
      IF ( IFAIL .NE. 0 ) RETURN
C
C   Check if we have to estimate ALPHA
C
      IF (ALPHA.LE.-SINGUL) THEN
	   CALL DECALP ( NDIM, A, B, NUMFUN, ALPHA,
     +     SINGUL, LOGF, FUNSUB, WORK(NW-NDIM-NUMFUN), WORK(NW-NUMFUN) )
C
C    Check the computed ALPHA-value.
C
          IF (ALPHA.LE.-SINGUL) THEN
            IFAIL=14
            RETURN
          END IF
      END IF
C
C     Compute the size of the temporary work space needed in DEADHR.
C
      LENW = 16*NUMFUN
C
C     Split up the work space.
C
      I1 = 1
      I2 = I1+(WRKSUB+1)*NUMFUN
      I3 = I2+(WRKSUB+1)*NUMFUN
      I4 = I3+(WRKSUB+1)*NDIM
      I5 = I4+(WRKSUB+1)*NDIM
      I6 = I5+WRKSUB+1
      I7 = I6+WRKSUB+1
      I8 = I7+LENW
      I9 = I8+(WRKSUB+1)*NUMFUN
      I10 = I9+WRKSUB*NUMFUN
      I11 = I10+WRKSUB*NUMFUN
      I12 = I11+NUMFUN
      I13 = I12+NUMFUN
      I14 = I13+NUMFUN
      I15 = I14+NUMFUN
      I16 = I15+NUMFUN
      I17 = I16+EMAX
      I18 = I17+(EMAX+1)**2
      I19 = I18+(EMAX+1)*NUMFUN
      K1 = 1
      K2 = K1+WTLENG*NDIM
      K3 = K2+WTLENG*5
      K4 = K3+WTLENG
      K5 = K4+NDIM
      K6 = K5+NDIM
      K7 = K6+NDIM
      K8 = K7+3*WTLENG
C
C   On restart runs the number of subregions from the
C   previous call is assigned to NSUB.
C
      IF ( RESTAR .EQ. 1 ) NSUB = WORK(NW)
      CALL DEADHR (NDIM,NUMFUN,A,B,MAXSUB,FUNSUB,SINGUL,ALPHA,LOGF,
     +     EPSABS,EPSREL,KEYF,RESTAR,NUM,LENW,WTLENG,EMAX,MINPTS,MAXPTS,
     +     NSUB,RESULT,ABSERR,NEVAL,IFAIL,WORK(I1),WORK(I2),WORK(I3),
     +     WORK(I4),WORK(I5),WORK(I6),WORK(I7),WORK(I8),WORK(I9),
     +     WORK(I10),WORK(I11),WORK(I12),WORK(I13),WORK(I14),WORK(I15),
     +     WORK(I16),WORK(I17),WORK(I18),WORK(I19),WORK2(K1),WORK2(K2),
     +     WORK2(K3),WORK2(K4),WORK2(K5),WORK2(K6),WORK2(K7),WORK2(K8),
     +     IWORK(1),IWORK(WRKSUB+1),IWORK(2*WRKSUB+1))
      WORK(NW) = NSUB
      RETURN
C
C***END DECUHR
C
      END
      SUBROUTINE DEFSHR (NDIM,CENTER,HWIDTH,X,G,NUMFUN,FUNSUB,FULSMS,
     +     FUNVLS)
C***BEGIN PROLOGUE DEFSHR
C***KEYWORDS fully symmetric sum
C***PURPOSE  To compute fully symmetric basic rule sums
C***AUTHOR
C            Alan Genz,
C            Department of Mathematics
C            Washington State University
C            Pullman, WA 99164-3113, USA
C            Email..  genz@gauss.math.wsu.edu
C***LAST MODIFICATION 88-04-08
C***DESCRIPTION DEFSHR computes a fully symmetric sum for a vector
C            of integrand values over a hyper-rectangular region.
C            The sum is fully symmetric with respect to the center of
C            the region and is taken over all sign changes and
C            permutations of the generators for the sum.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   CENTER Real array of dimension NDIM.
C          The coordinates for the center of the region.
C   HWIDTH Real Array of dimension NDIM.
C          HWIDTH(I) is half of the width of dimension I of the region.
C   X      Real Array of dimension NDIM.
C          A work array.
C   G      Real Array of dimension NDIM.
C          The generators for the fully symmetric sum. These MUST BE
C          non-negative and non-increasing.
C   NUMFUN Integer.
C          Number of components for the vector integrand.
C   FUNSUB Externally declared subroutine.
C          For computing the components of the integrand at a point X.
C          It must have parameters (NDIM, X, NUMFUN, FUNVLS).
C           Input Parameters:
C            X      Real array of dimension NDIM.
C                   Defines the evaluation point.
C            NDIM   Integer.
C                   Number of variables for the integrand.
C            NUMFUN Integer.
C                   Number of components for the vector integrand.
C           Output Parameters:
C            FUNVLS Real array of dimension NUMFUN.
C                   The components of the integrand at the point X.
C   ON RETURN
C
C   FULSMS Real array of dimension NUMFUN.
C          The values for the fully symmetric sums for each component
C          of the integrand.
C   FUNVLS Real array of dimension NUMFUN.
C          A work array.
C
C***ROUTINES CALLED: FUNSUB
C
C***END PROLOGUE DEFSHR
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN
      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM),X(NDIM),G(NDIM),
     +     FULSMS(NUMFUN),FUNVLS(NUMFUN)
C
C   Local variables.
C
      INTEGER IXCHNG,LXCHNG,I,J,L
      DOUBLE PRECISION GL,GI
C
C***FIRST EXECUTABLE STATEMENT DEFSHR
C
      DO 10 J=1,NUMFUN
         FULSMS(J)=0
 10   CONTINUE
C
C     Compute centrally symmetric sum for permutation of G
C
 20   DO 30 I=1,NDIM
         X(I)=CENTER(I)+G(I)*HWIDTH(I)
 30   CONTINUE
 40   CALL FUNSUB (NDIM,X,NUMFUN,FUNVLS)
      DO 50 J=1,NUMFUN
         FULSMS(J)=FULSMS(J)+FUNVLS(J)
 50   CONTINUE
      DO 60 I=1,NDIM
         G(I)=-G(I)
         X(I)=CENTER(I)+G(I)*HWIDTH(I)
         IF (G(I).LT.0) GO TO 40
 60   CONTINUE
C
C     Find next distinct permutation of G and loop back for next sum.
C     Permutations are generated in reverse lexicographic order.
C
      DO 80 I=2,NDIM
         IF (G(I-1).GT.G(I)) THEN
            GI=G(I)
            IXCHNG=I-1
            DO 70 L=1,(I-1)/2
               GL=G(L)
               G(L)=G(I-L)
               G(I-L)=GL
               IF (GL.LE.GI) IXCHNG=IXCHNG-1
               IF (G(L).GT.GI) LXCHNG=L
 70         CONTINUE
            IF (G(IXCHNG).LE.GI) IXCHNG=LXCHNG
            G(I)=G(IXCHNG)
            G(IXCHNG)=GI
            GO TO 20
         END IF
 80   CONTINUE
C
C     Restore original order to generators
C
      DO 90 I=1,NDIM/2
         GI=G(I)
         G(I)=G(NDIM-I+1)
         G(NDIM-I+1)=GI
 90   CONTINUE
C
C***END DEFSHR
C
      END
      SUBROUTINE DEINHR (NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C***BEGIN PROLOGUE DEINHR
C***PURPOSE DEINHR computes abscissas and weights of the integration
C            rule and the null rules to be used in error estimation.
C            These are computed as functions of NDIM and KEY.
C***LAST MODIFICATION 88-04-08
C***DESCRIPTION DEINHR will for given values of NDIM and KEY compute or
C            select the correct values of the abscissas and
C            corresponding weights for different
C            integration rules and null rules and assign them to
C            G and W.
C            The heuristic error coefficients ERRCOF
C            will be computed as a function of KEY.
C            Scaling factors SCALES and normalization factors NORMS
C            used in the error estimation are computed.
C
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables.
C     KEY    Integer.
C            Key to selected local integration rule.
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C            It is assumed that the error is computed using:
C             IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
C               THEN ERROR = ERRCOF(3)*N1
C               ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
C             ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
C            where N1-N3 are the null rules, EP is the error for
C            the parent
C            subregion and ES is the error for the sibling subregion.
C     RULPTS Real array of dimension WTLENG.
C            A work array containing the number of points produced by
C            each generator of the selected rule.
C     SCALES Real array of dimension (3,WTLENG).
C            Scaling factors used to construct new null rules,
C            N1, N2 and N3,
C            based on a linear combination of two successive null rules
C            in the sequence of null rules.
C     NORMS  Real array of dimension (3,WTLENG).
C            2**NDIM/(1-norm of the null rule constructed by each of
C            the scaling factors.)
C
C***ROUTINES CALLED  D132RE,D113RE,D07HRE,D09HRE
C***END PROLOGUE DEINHR
C
C   Global variables.
C
      INTEGER NDIM,KEY,WTLENG
      DOUBLE PRECISION G(NDIM,WTLENG),W(5,WTLENG),ERRCOF(6)
      DOUBLE PRECISION RULPTS(WTLENG),SCALES(3,WTLENG)
      DOUBLE PRECISION NORMS(3,WTLENG)
C
C   Local variables.
C
      INTEGER I,J,K
      DOUBLE PRECISION WE(14)
C
C***FIRST EXECUTABLE STATEMENT DEINHR
C
C   Compute W, G and ERRCOF.
C
      IF (KEY.EQ.1) THEN
         CALL D132RE (WTLENG,W,G,ERRCOF,RULPTS)
      ELSE IF (KEY.EQ.2) THEN
         CALL D113RE (WTLENG,W,G,ERRCOF,RULPTS)
      ELSE IF (KEY.EQ.3) THEN
         CALL D09HRE (NDIM,WTLENG,W,G,ERRCOF,RULPTS)
      ELSE IF (KEY.EQ.4) THEN
         CALL D07HRE (NDIM,WTLENG,W,G,ERRCOF,RULPTS)
      END IF
C
C   Compute SCALES and NORMS.
C
      DO 40 K=1,3
         DO 30 I=1,WTLENG
            IF (ABS(W(K+1,I)).GT.0) THEN
               SCALES(K,I)=-W(K+2,I)/W(K+1,I)
            ELSE
               SCALES(K,I)=100
            END IF
            DO 10 J=1,WTLENG
               WE(J)=W(K+2,J)+SCALES(K,I)*W(K+1,J)
 10         CONTINUE
            NORMS(K,I)=0
            DO 20 J=1,WTLENG
               NORMS(K,I)=NORMS(K,I)+RULPTS(J)*ABS(WE(J))
 20         CONTINUE
            NORMS(K,I)=2**NDIM/NORMS(K,I)
 30      CONTINUE
 40   CONTINUE
C
C***END DEINHR
C
      END
      SUBROUTINE DERLHR (NDIM,CENTER,HWIDTH,WTLENG,G,W,ERRCOF,NUMFUN,
     +FUNSUB,SCALES,NORMS,X,NULL,BASVAL,RGNERR,DIRECT,GREATE,DIFF,ORDER)
C***BEGIN PROLOGUE DERLHR
C***KEYWORDS basic numerical integration rule
C***PURPOSE  To compute basic integration rule values.
C***LAST MODIFICATION 92-09-16
C***DESCRIPTION DERLHR computes basic integration rule values for a
C            vector of integrands over a hyper-rectangular region.
C            These are estimates for the integrals. DERLHR also computes
C            estimates for the errors and determines the coordinate axis
C            where the fourth difference for the integrands is largest.
C            In case the region contains the singularity DERLHR creates
C            a pointer array ordered according to the size of the fourth
C            difference along the first -DIRECT coordinate axes.
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   CENTER Real array of dimension NDIM.
C          The coordinates for the center of the region.
C   HWIDTH Real Array of dimension NDIM.
C          HWIDTH(I) is half of the width of dimension I of the region.
C   WTLENG Integer.
C          The number of weights in the basic integration rule.
C   G      Real array of dimension (NDIM,WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1,J), ..., G(NDIM,J) are the are the generators for the
C          points associated with the Jth weights.
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   ERRCOF Real array of dimension 6.
C          The error coefficients for the rules.
C          It is assumed that the error is computed using:
C           IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
C             THEN ERROR = ERRCOF(3)*N1
C             ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
C           ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
C          where N1-N4 are the null rules, EP is the error
C          for the parent
C          subregion and ES is the error for the sibling subregion.
C   NUMFUN Integer.
C          Number of components for the vector integrand.
C   FUNSUB Externally declared subroutine.
C          For computing the components of the integrand at a point X.
C          It must have parameters (NDIM,X,NUMFUN,FUNVLS).
C           Input Parameters:
C            X      Real array of dimension NDIM.
C                   Defines the evaluation point.
C            NDIM   Integer.
C                   Number of variables for the integrand.
C            NUMFUN Integer.
C                   Number of components for the vector integrand.
C           Output Parameters:
C            FUNVLS Real array of dimension NUMFUN.
C                   The components of the integrand at the point X.
C   SCALES Real array of dimension (3,WTLENG).
C          Scaling factors used to construct new null rules based
C          on a linear combination of two successive null rules
C          in the sequence of null rules.
C   NORMS  Real array of dimension (3,WTLENG).
C          2**NDIM/(1-norm of the null rule constructed by each of the
C          scaling factors.)
C   X      Real Array of dimension NDIM.
C          A work array.
C   NULL   Real array of dimension (NUMFUN, 8)
C          A work array.
C   DIFF   Real array of dimension NDIM.
C          A work array.
C   ORDER  Integer array of dimension NDIM.
C          Assumed initially: ORDER(I) = I, I=1,2,...,NDIM. Later
C          we assume ORDER to be unchanged from the last call of
C          DERLHR.
C   ON RETURN
C
C   BASVAL Real array of dimension NUMFUN.
C          The values for the basic rule for each component
C          of the integrand.
C   RGNERR Real array of dimension NUMFUN.
C          The error estimates for each component of the integrand.
C   DIRECT Real.
C          The coordinate axis where the fourth difference of the
C          integrand values is largest.
C   GREATE Real
C          The greatest error estimate among the components.
C   DIFF   Real array of dimension NDIM.
C          Stores the fourth differences along the axes (only for the
C          singular region).
C   ORDER  Integer array of order NDIM.
C          Pointer array: DIFF(ORDER(I)) < DIFF(ORDER(I-1)), for
C          I=2,3,...,-DIRECT. Used only in connection with the
C          singular region.
C
C***REFERENCES
C
C   A.C.Genz and A.A.Malik, An adaptive algorithm for numerical
C   integration over an N-dimensional rectangular region,
C   J.Comp.Appl.Math., 6:295-302, 1980.
C
C   T.O.Espelid, Integration Rules, Null Rules and Error
C   Estimation, Reports in Informatics 33, Dept. of Informatics,
C   Univ. of Bergen, 1988.
C
C   T.O.Espelid and A.Genz, DECUHR: An Algorithm for Automatic
C   Integration of Singular Functions over a Hyperrectangular Region.
C   Numerical Algorithms 8(1994), PP. 201-220.
C
C***ROUTINES CALLED: DEFSHR, FUNSUB
C
C***END PROLOGUE DERLHR
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER WTLENG,NUMFUN,NDIM,ORDER(NDIM)
      DOUBLE PRECISION CENTER(NDIM),X(NDIM),HWIDTH(NDIM),BASVAL(NUMFUN),
     +RGNERR(NUMFUN),NULL(NUMFUN,8),W(5,WTLENG),G(NDIM,WTLENG),ERRCOF(6)
     +,DIRECT,SCALES(3,WTLENG),NORMS(3,WTLENG),GREATE,DIFF(NDIM)
C
C   Local variables.
C
      DOUBLE PRECISION RGNVOL,DIFSUM,DIFMAX,FRTHDF
      INTEGER I,J,K,DIVAXN,SWAP
      DOUBLE PRECISION SEARCH,RATIO,NEW
C
C***FIRST EXECUTABLE STATEMENT DERLHR
C
C
C       Compute volume of subregion, initialize DIVAXN and rule sums;
C       compute fourth differences and new DIVAXN (RGNERR is used
C       for a work array here). The integrand values used for the
C       fourth divided differences are accumulated in rule arrays.
C
      RGNVOL=1
      DIVAXN=1
      DO 10 I=1,NDIM
         RGNVOL=RGNVOL*HWIDTH(I)
         X(I)=CENTER(I)
         IF (HWIDTH(I).GT.HWIDTH(DIVAXN)) DIVAXN=I
 10   CONTINUE
      CALL FUNSUB (NDIM,X,NUMFUN,RGNERR)
      DO 30 J=1,NUMFUN
         BASVAL(J)=W(1,1)*RGNERR(J)
         DO 20 K=1,4
            NULL(J,K)=W(K+1,1)*RGNERR(J)
 20      CONTINUE
 30   CONTINUE
      DIFMAX=0
      RATIO=(G(1,3)/G(1,2))**2
      DO 60 I=1,NDIM
         X(I)=CENTER(I)-HWIDTH(I)*G(1,2)
         CALL FUNSUB (NDIM,X,NUMFUN,NULL(1,5))
         X(I)=CENTER(I)+HWIDTH(I)*G(1,2)
         CALL FUNSUB (NDIM,X,NUMFUN,NULL(1,6))
         X(I)=CENTER(I)-HWIDTH(I)*G(1,3)
         CALL FUNSUB (NDIM,X,NUMFUN,NULL(1,7))
         X(I)=CENTER(I)+HWIDTH(I)*G(1,3)
         CALL FUNSUB (NDIM,X,NUMFUN,NULL(1,8))
         X(I)=CENTER(I)
         DIFSUM=0
         DO 50 J=1,NUMFUN
            FRTHDF = ABS( 2*(1-RATIO)*RGNERR(J)-(NULL(J,7)+NULL(J,8))
     +           +RATIO*(NULL(J,5)+NULL(J,6)) )
C
C     Ignore differences below roundoff
C
            IF (ABS(RGNERR(J))+FRTHDF/4.GT.ABS(RGNERR(J)))
     +           DIFSUM = DIFSUM + FRTHDF
            DO 40 K=1,4
               NULL(J,K)=NULL(J,K)+W(K+1,2)*(NULL(J,5)+NULL(J,6))+W(K+1,
     +          3)*(NULL(J,7)+NULL(J,8))
 40         CONTINUE
            BASVAL(J)=BASVAL(J)+W(1,2)*(NULL(J,5)+NULL(J,6))+W(1,3)*
     +       (NULL(J,7)+NULL(J,8))
 50      CONTINUE
         IF (DIFSUM.GT.DIFMAX) THEN
            DIFMAX=DIFSUM
            DIVAXN=I
         END IF
         DIFF(I)=DIFSUM
 60   CONTINUE
      IF (DIRECT.GE.0) DIRECT=DIVAXN
C
C     Finish computing the rule values.
C
      DO 90 I=4,WTLENG
         CALL DEFSHR (NDIM,CENTER,HWIDTH,X,G(1,I),NUMFUN,FUNSUB,RGNERR,
     +    NULL(1,5))
         DO 80 J=1,NUMFUN
            BASVAL(J)=BASVAL(J)+W(1,I)*RGNERR(J)
            DO 70 K=1,4
               NULL(J,K)=NULL(J,K)+W(K+1,I)*RGNERR(J)
 70         CONTINUE
 80      CONTINUE
 90   CONTINUE
C
C     Compute errors and find the great.
C
      GREATE=0
      DO 120 J=1,NUMFUN
C
C     We search for the null rule, in the linear space spanned by two
C     successive null rules in our sequence, which gives the greatest
C     error estimate among all normalized (1-norm) null rules in this
C     space.
C
         DO 110 I=1,3
            SEARCH=0
            DO 100 K=1,WTLENG
               SEARCH=MAX(SEARCH,ABS(NULL(J,I+1)+SCALES(I,K)*NULL(J,I))*
     +          NORMS(I,K))
 100        CONTINUE
            NULL(J,I)=SEARCH
 110     CONTINUE
         IF (ERRCOF(1)*NULL(J,1).LE.NULL(J,2).AND.ERRCOF(2)*NULL(J,2).
     +    LE.NULL(J,3)) THEN
            RGNERR(J)=ERRCOF(3)*NULL(J,1)
         ELSE
            RGNERR(J)=ERRCOF(4)*MAX(NULL(J,1),NULL(J,2),NULL(J,3))
         END IF
         RGNERR(J)=RGNVOL*RGNERR(J)
         BASVAL(J)=RGNVOL*BASVAL(J)
         GREATE=MAX(GREATE,RGNERR(J))
 120  CONTINUE
C
C     If this region is the singular one it will be cut -DIRECT
C     times in directions: 1,2,...,-DIRECT. To decide in which
C     order this subdivision will take place we will permute
C     the integers such that DIFF(ORDER(I)),I=1,2,...,-DIRECT
C     give the sequence in decreasing order. ORDER is initially
C     1,2,...,-DIRECT, but will later be set to the values used by
C     the previous singular region. The sorting code used is
C     sifting from below.
C
      IF (DIRECT.LT.0) THEN
         IDIRECT = DIRECT
         DO 140 I=2,-IDIRECT
            NEW=DIFF(ORDER(I))
            DO 130 J=I-1,1,-1
               IF (NEW.GT.DIFF(ORDER(J))) THEN
C
C     Swap pointers
C
                  SWAP=ORDER(J+1)
                  ORDER(J+1)=ORDER(J)
                  ORDER(J)=SWAP
               ELSE
                  GO TO 140
               END IF
 130        CONTINUE
 140     CONTINUE
      END IF
C
C***END DERLHR
C
      END
      SUBROUTINE DESBHR (NDIM,CENTRS,HWIDTS,DIR,DIRECT,POINTR,VACANT,
     +     ORDER)
C***BEGIN PROLOGUE DESBHR
C***REFER TO DECUHR
C***PURPOSE DESBHR cuts the given region one or more times.
C***LAST MODIFICATION 92-09-16
C***DESCRIPTION The number of new subregions depends on the region.
C
C           Each cut creates two new halves of the current region.
C
C           The singular region will be cut -DIRECT times in
C           directions 1,2,3,...,-DIRECT. The order of these
C           -DIRECT cuts is given by ORDER.
C           This gives -DIRECT + 1 new subregions.
C
C           A non-singular region will be cut once in direction
C           DIRECT. This gives two new subregions.
C***ROUTINES CALLED-NONE
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables.
C     CENTRS Real array of dimension (NDIM,0:*).
C            Used to store the centers of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,0:*).
C            Used to store the half widths of the stored subregions.
C     DIR    Real array of dimension (0:*).
C            DIR is used to store the directions for
C            further subdivision.
C     DIRECT If positive: the direction of subdivision.
C            If negative: cut direction 1,2,...,-DIRECT.
C            This will be the case if this region is the singular one.
C     POINTR Pointer to the position in the data structure where
C            the new subregions are to be stored.
C     VACANT Pointer to the region which has to be subdivided.
C            One of the new reions will be stored at this vacant element
C     ORDER  Integer array of pointers.
C            In case this is the singular region it will be cut in
C            directions 1,2,...,-DIRECT, in the order given by this
C            pointer array.
C
C   ON RETURN
C
C     CENTRS Real array of dimension (NDIM,0:*).
C            Unchanged array except for position VACANT and
C            positions POINTR, POINTR+1,... (depending on the
C            number of new subregions).
C     HWIDTS Real array of dimension (NDIM,0:*).
C            Unchanged array except for position VACANT and
C            positions POINTR, POINTR+1,... (depending of the
C            number of new subregions).
C     DIR    Real array of dimension (0:*)
C            Unchanged array except for position VACANT and
C            positions POINTR, POINTR+1,... (depending of the
C            number of new subregions). Saves the last direction
C            these subregions are cut. Except for the singular region
C            where the value is left unchanged equal to DIRECT.
C
C***END PROLOGUE DESBHR
C
C   Global variables.
C
      INTEGER NDIM,POINTR,DIRECT,VACANT,ORDER(NDIM)
      DOUBLE PRECISION CENTRS(NDIM,0:*),HWIDTS(NDIM,0:*),DIR(0:*)
C
C   Local variables.
C
      INTEGER I,J,NEXT
C
C   FIRST EXECUTABLE STATEMENT DESBHR
C
      IF (DIRECT.GT.0) THEN
C
C     When DIRECT is positive, indicating a non-singular region,
C     we only cut the region in two halves.
C
         DO 10 J=1,NDIM
            HWIDTS(J,POINTR)=HWIDTS(J,VACANT)
            CENTRS(J,POINTR)=CENTRS(J,VACANT)
 10      CONTINUE
         HWIDTS(DIRECT,VACANT)=HWIDTS(DIRECT,VACANT)/2
         HWIDTS(DIRECT,POINTR)=HWIDTS(DIRECT,VACANT)
         CENTRS(DIRECT,POINTR)=CENTRS(DIRECT,VACANT)+HWIDTS(DIRECT,
     +    POINTR)
         CENTRS(DIRECT,VACANT)=CENTRS(DIRECT,VACANT)-HWIDTS(DIRECT,
     +    VACANT)
         DIR(POINTR)=DIRECT
      ELSE
C
C     DIRECT negative signals the singular region: the absolute value
C     tells us how many co-ordinates are involved (always assumed to
C     be 1,2,...,-DIRECT).
C     In this case we slice the singular region -DIRECT times
C     in directions 1,2,...,-DIRECT. These cuts are performed according
C     to the array ORDER. We get -DIRECT+1 new pieces
C     and the reduced singular region ends up in the
C     the same position as the previous one: thus VACANT = 0.
C
         DO 30 I=1,-DIRECT
            NEXT=ORDER(I)
            HWIDTS(NEXT,VACANT)=HWIDTS(NEXT,VACANT)/2
            DO 20 J=1,NDIM
               HWIDTS(J,POINTR+I-1)=HWIDTS(J,VACANT)
               CENTRS(J,POINTR+I-1)=CENTRS(J,VACANT)
 20         CONTINUE
            DIR(POINTR+I-1)=NEXT
            CENTRS(NEXT,POINTR+I-1)=CENTRS(NEXT,VACANT)
     +           + HWIDTS(NEXT,VACANT)
            CENTRS(NEXT,VACANT)=CENTRS(NEXT,VACANT)-HWIDTS(NEXT,VACANT)
 30      CONTINUE
      END IF
      DIR(VACANT)=DIRECT
      END
      SUBROUTINE DETRHR (DVFLAG,SBRGNS,GREATE,LIST,NEW)
C***BEGIN PROLOGUE DETRHR
C***REFER TO DECUHR
C***PURPOSE  DETRHR maintains a heap of subregions.
C***AUTHOR   Terje O. Espelid, Department of Informatics,
C            University of Bergen,  Hoyteknologisenteret,
C            N-5020 Bergen, Norway
C            Email..  terje@ii.uib.no
C***LAST MODIFICATION 92-09-16
C***DESCRIPTION DETRHR maintains a heap of subregions.
C            The subregions are stored in a partially sorted
C            binary tree, ordered according to the size of the
C            greatest error estimates of each subregion(GREATE).
C            The subregion with greatest error estimate is in the
C            first position of the heap.
C
C   PARAMETERS
C
C     DVFLAG Integer.
C            If DVFLAG = 1, we remove the subregion with
C            greatest error from the heap.
C            If DVFLAG = 2, we insert a new subregion in the heap.
C     SBRGNS Integer.
C            Number of subregions in the heap.
C     GREATE Real array of dimension SBRGNS.
C            Used to store the greatest estimated errors in
C            all subregions.
C     LIST   Integer array of dimension SBRGNS.
C            Used as a partially ordered list of pointers to the
C            different subregions. This list is a heap where the
C            element on top of the list is the subregion with the
C            greatest error estimate.
C     NEW    Integer.
C            Index to the new region to be inserted in the heap.
C***ROUTINES CALLED-NONE
C***END PROLOGUE DETRHR
C
C   Global variables.
C
      INTEGER DVFLAG,NEW,SBRGNS,LIST(*)
      DOUBLE PRECISION GREATE(0:*)
C
C   Local variables.
C
C   GREAT  is used as intermediate storage for the greatest error of a
C          subregion.
C   SUBRGN Position of child/parent subregion in the heap.
C   SUBTMP Position of parent/child subregion in the heap.
      INTEGER SUBRGN,SUBTMP
      DOUBLE PRECISION GREAT
C
C***FIRST EXECUTABLE STATEMENT DTRTRI
C
C     If DVFLAG = 1, we will reduce the partial ordered list by the
C     element with greatest estimated error. Thus the element in
C     in the heap with index LIST(1) is vacant and can be used later.
C     Reducing the heap by one element implies that the last element
C     should be re-positioned.
C
      IF (DVFLAG.EQ.1) THEN
         GREAT=GREATE(LIST(SBRGNS))
         SBRGNS=SBRGNS-1
         SUBRGN=1
 10      SUBTMP=2*SUBRGN
         IF (SUBTMP.LE.SBRGNS) THEN
            IF (SUBTMP.NE.SBRGNS) THEN
C
C     Find max. of left and right child.
C
               IF (GREATE(LIST(SUBTMP)).LT.GREATE(LIST(SUBTMP+1))) THEN
                  SUBTMP=SUBTMP+1
               END IF
            END IF
C
C     Compare max.child with parent.
C     If parent is max., then done.
C
            IF (GREAT.LT.GREATE(LIST(SUBTMP))) THEN
C
C     Move the pointer at position SUBTMP up the heap.
C
               LIST(SUBRGN)=LIST(SUBTMP)
               SUBRGN=SUBTMP
               GO TO 10
            END IF
         END IF
C
C     Update the pointer.
C
         IF (SBRGNS.GT.0) THEN
            LIST(SUBRGN)=LIST(SBRGNS+1)
         END IF
      ELSE IF (DVFLAG.EQ.2) THEN
C
C     If DVFLAG = 2, find the position for the NEW region in the heap.
C
         GREAT=GREATE(NEW)
         SUBRGN=SBRGNS
 20      SUBTMP=SUBRGN/2
         IF (SUBTMP.GE.1) THEN
C
C     Compare max.child with parent.
C     If parent is max, then done.
C
            IF (GREAT.GT.GREATE(LIST(SUBTMP))) THEN
C
C     Move the pointer at position SUBTMP down the heap.
C
               LIST(SUBRGN)=LIST(SUBTMP)
               SUBRGN=SUBTMP
               GO TO 20
            END IF
         END IF
C
C     Set the pointer to the new region in the heap.
C
         LIST(SUBRGN)=NEW
      END IF
C
C***END DETRHR
C
      END
      SUBROUTINE DEXTHR (NUMFUN,ALPHA,LOGF,SINGUL,EXSTEP,N,T,UPDATE,
     +UNEW,QNEW,QOLD,FIRST,EMAX,UERR,RESULT,ABSERR,EXTERR,NE,BETA)
C***BEGIN PROLOGUE DEXTHR
C***KEYWORDS Linear extrapolation, homogeneous functions,
C     logarithmic singularities, error estimation.
C***PURPOSE To compute better estimates to a vector of approximations
C     to multidimensional integrals and to provide new and
C     updated error estimates.
C***LAST MODIFICATION 92-09-16
C***DESCRIPTION
C            The routine uses linear extrapolation to compute better
C            approximations to each component in a vector of
C            multidimensional integrals. All components are assumed to
C            be singular due to a homogeneous function of degree ALPHA
C            in a vertex of dimension SINGUL. In addition we may have
C            a logarithmic singularity in the same vertex(dim. SINGUL).
C            A series, with tail correction, approach is used, assuming
C            that the terms are given with estimates of the error in
C            each term. New error estimates are computed too. The
C            routine have two options: either a new extrapolation term
C            is provided and we take a new extrapolation step,
C            or we update one previously computed term in the series and
C            therefore have to update the extrapolation tableau.
C
C   ON ENTRY
C
C     NUMFUN Integer.
C            Number of components of the integral.
C     ALPHA  Real.
C            Degree of homogeneous function: ALPHA > -SINGUL.
C            This singularity is assumed the same in all integrands.
C     LOGF    Integer
C            If LOGF = 1 then there is a logarithmic singularity in all
C            integrands, else there is no logarithmic singularity.
C            We assume the logarithm to appear only in power of 1.
C     SINGUL Integer
C            The dimension of the singularity.
C     EXSTEP Integer.
C            The exponent sequence's step size in the error expansion.
C     N      Integer
C            The number of U-terms in the series.
C     T      Real array of dimension (NUMFUN,0:EMAX)
C            Contains the last row in the extrapolation tableau for
C            each function in the vector.
C     UPDATE Integer
C            = 0 then this is a new extrapolation step.
C            > 0 then this is  a step where we have to correct the
C            existing tableau. The value of UPDATE gives the index to th
C            u-elment that has been modified.
C     UNEW   Real array of dimension NUMFUN.
C            If UPDATE = 0 then this gives the next terms in the series,
C            else it is the correction to the u-values with index UPDATE
C     QNEW   Real array of dimension NUMFUN.
C            The new tail correction for all functions in the vector.
C     QOLD   Real array of dimension NUMFUN.
C            The last tail corrections for all functions in the vector.
C     FIRST  Logical.
C            Value .TRUE. indicates that this is the first time
C            this routine is called.
C     EMAX   Integer
C            The maximum allowed number of extrapolation steps.
C     UERR   Real array of dimension (NUMFUN,N)
C            The estimated errors of all U-terms in the series.
C   ON RETURN
C
C     T      Real array of dimension (NUMFUN,0:EMAX)
C            Contains the last row in the extrapolation tableau for
C            each function in the vector after the extrapolation. In
C            case this is an updating steps, then each element may
C            have been changed.
C     RESULT Real array of dimension NUMFUN
C            Contains the new approximations for the components
C            of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Returns the global errors for all components.
C            This includes both the pure extrapolation
C            error and the effect of not knowing the U-terms exactly.
C     EXTERR Real array of dimension NUMFUN.
C            These errors are associated with the singular region and
C            they are the pure extrapolation errors.
C     NE     Real array of dimension (0:EMAX).
C            A table of denominators to be used in the extrapolation.
C            Dummy parameter
C     BETA   Real Array of dimension (EMAX +1)(EMAX +1)
C            A table of coefficients to be used in the error estimation.
C
C***REFERENCES
C
C   T.O.Espelid and A.Genz, DECUHR: An Algorithm for Automatic
C   Integration of Singular Functions over a Hyperrectangular Region.
C   Numerical Algorithms 8(1994), PP. 201-220.
C
C   T.O.Espelid, On integrating Vertex Singularities using
C   Extrapolation,  BIT 34,1(1994) 62-79.
C
C   T.O.Espelid, On integrating Singularities using non-Uniform
C   Subdivision and Extrapolation,  in Numerical Integration IV, eds.
C   Brass and Hammerlin, Birhauser, ISNM Vol. 112(1993) 77-89.
C
C***ROUTINES CALLED NONE
C***END PROLOGUE DEXTHR
C
C   Global variables.
C
      INTEGER N,SINGUL,EMAX,UPDATE,NUMFUN,LOGF,EXSTEP
      DOUBLE PRECISION ALPHA,T(NUMFUN,0:EMAX),QOLD(NUMFUN),QNEW(NUMFUN),
     +UNEW(NUMFUN),UERR(NUMFUN,N),NE(EMAX),BETA(0:EMAX,0:EMAX),
     +RESULT(NUMFUN),ABSERR(NUMFUN),EXTERR(NUMFUN)
      LOGICAL FIRST
C
C   Local variables.
C
      INTEGER I,J,STEPS
      DOUBLE PRECISION SAVE1,SAVE2,CONST,ES
      PARAMETER (CONST=10)
C
C  CONST heuristic constant used in the error estimation.
C  STEP  integer; keeping track of the number of extrapolation
C        steps actually used.
C
C***FIRST EXECUTABLE STATEMENT DEXTHR
C
C     Check if this is the first time the routine is called
C
      IF (FIRST) THEN
C
C     Initialize the extrapolation tableaus
C
         DO 10 J=1,NUMFUN
            T(J,0)=QOLD(J)
 10      CONTINUE
C
C     Compute all denominators and take into account if there
C     is a logarithmic term present.
C
         NE(1)=2**(SINGUL+ALPHA)-1
         IF (LOGF.EQ.1) THEN
            DO 30 I=3,EMAX,2
               ES=NE(I-2)
               DO 20 J=1,EXSTEP
                  ES=2*ES+1
 20            CONTINUE
               NE(I)=ES
 30         CONTINUE
            DO 40 I=2,EMAX,2
               NE(I)=NE(I-1)
 40         CONTINUE
         ELSE
C
C     assuming log = 0, (note log > 1 has not been implemented)
C
            DO 60 I=2,EMAX,1
               ES=NE(I-1)
               DO 50 J=1,EXSTEP
                  ES=2*ES+1
 50            CONTINUE
               NE(I)=ES
 60         CONTINUE
         END IF
C
C     Initialize the beta-factors to be used in the error estimation.
C
         DO 90 J=0,EMAX
            BETA(0,J)=1
            DO 70 I=1,J
               BETA(I,J)=BETA(I,J-1)+(BETA(I,J-1)-BETA(I-1,J-1))/NE(J)
 70         CONTINUE
            DO 80 I=J+1,EMAX
               BETA(I,J)=0
 80         CONTINUE
 90      CONTINUE
      END IF
C
C     The number of extrapolation steps
C
      STEPS=MIN(N,EMAX)
C
C     Check what kind of step this is: new extrapolation or modifying
C
      IF (UPDATE.EQ.0) THEN
C
C     A new extrapolation step.
C
         DO 110 J=1,NUMFUN
            SAVE1=T(J,0)+(UNEW(J)+(QNEW(J)-QOLD(J)))
            DO 100 I=1,STEPS
               SAVE2=SAVE1+(SAVE1-T(J,I-1))/NE(I)
               T(J,I-1)=SAVE1
               SAVE1=SAVE2
 100        CONTINUE
            T(J,STEPS)=SAVE1
 110     CONTINUE
C
C     A modification step.
C
      ELSE IF (UPDATE.LT.N-STEPS) THEN
C
C     Simply add the correction to all elements in the tableau.
C
         DO 130 J=1,NUMFUN
            DO 120 I=0,STEPS
               T(J,I)=T(J,I)+UNEW(J)
 120        CONTINUE
 130     CONTINUE
      ELSE
         DO 150 J=1,NUMFUN
            DO 140 I=0,STEPS
               T(J,I)=T(J,I)+UNEW(J)*(1-BETA(N-UPDATE+1,I))
 140        CONTINUE
 150     CONTINUE
      END IF
C
C     Then compute the error estimates.
C     First the extrapolation error and then the U-effect.
C     The error is accumulated in qnew
C
      DO 180 J=1,NUMFUN
         EXTERR(J)=CONST*ABS(T(J,STEPS)-T(J,STEPS-1))
         QNEW(J)=EXTERR(J)
C
C     Note: The last U-errors are effected by the extrapolation-process
C
         DO 160 I=1,STEPS
            QNEW(J)=QNEW(J)+ABS(1-BETA(I,STEPS))*UERR(J,N+1-I)
 160     CONTINUE
         DO 170 I=STEPS+1,N
            QNEW(J)=QNEW(J)+UERR(J,N+1-I)
 170     CONTINUE
 180  CONTINUE
C
C   Define the results and the new errors. We update only those
C   components which have an improved error estimate.
C
      DO 190 J=1,NUMFUN
	  IF ( QNEW(J).LE.ABSERR(J)) THEN
            RESULT(J) = T(J,STEPS)
       	    ABSERR(J) = QNEW(J)
          END IF
 190  CONTINUE
C
C***END DEXTHR
C
      END