      SUBROUTINE DDAJAC (NEQ, X, Y, YPRIME, DELTA, CJ, H,
     +   IER, WT, E, WM, IWM, RES, IRES, UROUND, JAC, RPAR,
     +   IPAR, NTEMP)
c-----------------------------------------------------------------------
c     This routine has been modified to accept an almost block diagonal
c     (ABD) Jacobian matrix compatible with the linear algebra package,
c     COLROW.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, August 13, 2001.
c
c-----------------------------------------------------------------------
C***BEGIN PROLOGUE  DDAJAC
C***SUBSIDIARY
C***PURPOSE  Compute the iteration matrix for DDASSL and form the
C            LU-decomposition.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDAJAC-S, DDAJAC-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------------
C     THIS ROUTINE COMPUTES THE ITERATION MATRIX
C     PD=DG/DY+CJ*DG/DYPRIME (WHERE G(X,Y,YPRIME)=0).
C     HERE PD IS COMPUTED BY THE USER-SUPPLIED
C     ROUTINE JAC IF IWM(MTYPE) IS 1 OR 4, AND
C     IT IS COMPUTED BY NUMERICAL FINITE DIFFERENCING
C     IF IWM(MTYPE)IS 2 OR 5
C     THE PARAMETERS HAVE THE FOLLOWING MEANINGS.
C     Y        = ARRAY CONTAINING PREDICTED VALUES
C     YPRIME   = ARRAY CONTAINING PREDICTED DERIVATIVES
C     DELTA    = RESIDUAL EVALUATED AT (X,Y,YPRIME)
C                (USED ONLY IF IWM(MTYPE)=2 OR 5)
C     CJ       = SCALAR PARAMETER DEFINING ITERATION MATRIX
C     H        = CURRENT STEPSIZE IN INTEGRATION
C     IER      = VARIABLE WHICH IS .NE. 0
C                IF ITERATION MATRIX IS SINGULAR,
C                AND 0 OTHERWISE.
C     WT       = VECTOR OF WEIGHTS FOR COMPUTING NORMS
C     E        = WORK SPACE (TEMPORARY) OF LENGTH NEQ
C     WM       = REAL WORK SPACE FOR MATRICES. ON
C                OUTPUT IT CONTAINS THE LU DECOMPOSITION
C                OF THE ITERATION MATRIX.
C     IWM      = INTEGER WORK SPACE CONTAINING
C                MATRIX INFORMATION
C     RES      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
C                TO EVALUATE THE RESIDUAL FUNCTION G(X,Y,YPRIME)
C     IRES     = FLAG WHICH IS EQUAL TO ZERO IF NO ILLEGAL VALUES
C                IN RES, AND LESS THAN ZERO OTHERWISE.  (IF IRES
C                IS LESS THAN ZERO, THE MATRIX WAS NOT COMPLETED)
C                IN THIS CASE (IF IRES .LT. 0), THEN IER = 0.
C     UROUND   = THE UNIT ROUNDOFF ERROR OF THE MACHINE BEING USED.
C     JAC      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
C                TO EVALUATE THE ITERATION MATRIX (THIS ROUTINE
C                IS ONLY USED IF IWM(MTYPE) IS 1 OR 4)
C-----------------------------------------------------------------------
C***ROUTINES CALLED  DGBFA, BDGEFA
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901010  Modified three MAX calls to be all on one line.  (FNF)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C   901101  Corrected PURPOSE.  (FNF)
C***END PROLOGUE  DDAJAC
C
      INTEGER  NEQ, IER, IWM(*), IRES, IPAR(*), NTEMP
      DOUBLE PRECISION
     *   X, Y(*), YPRIME(*), DELTA(*), CJ, H, WT(*), E(*), WM(*),
     *   UROUND, RPAR(*)
      EXTERNAL  RES, JAC
C
      EXTERNAL  DGBFA, BDGEFA
C
      INTEGER  I, I1, I2, II, IPSAVE, ISAVE, J, K, L, LENPD, LIPVT,
     *   LML, LMTYPE, LMU, MBA, MBAND, MEB1, MEBAND, MSAVE, MTYPE, N,
     *   NPD, NPDM1, NROW
      DOUBLE PRECISION  DEL, DELINV, SQUR, YPSAVE, YSAVE
c-----------------------------------------------------------------------
      integer nconti
      parameter (nconti = 2)
      double precision zero
      parameter (zero = 0.d0)

      integer lnpde, lkcol, lnint
      integer npde, kcol, nint
      integer neq1, neq2
      integer npdbk1, npdbt1
      integer npdtp2, npdbk2, npdbt2
      integer lipvt2

c-----------------------------------------------------------------------
C
      PARAMETER (NPD=1)
      PARAMETER (LML=1)
      PARAMETER (LMU=2)
      PARAMETER (LMTYPE=4)
c-----------------------------------------------------------------------
      parameter (lnpde = 17)
      parameter (lkcol = 18)
      parameter (lnint = 19)
      parameter (lipvt = 21)
c-----------------------------------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT  DDAJAC
      IER = 0
      NPDM1=NPD-1
      MTYPE=IWM(LMTYPE)
      GO TO (100,200,300,400,500),MTYPE
C
C
C     DENSE USER-SUPPLIED MATRIX
100   LENPD=NEQ*NEQ
      DO 110 I=1,LENPD
110      WM(NPDM1+I)=0.0D0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      GO TO 230
C
C
C     DENSE FINITE-DIFFERENCE-GENERATED MATRIX
200   IRES=0
      NROW=NPDM1
      SQUR = SQRT(UROUND)
      DO 210 I=1,NEQ
         DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)),ABS(WT(I)))
         DEL=SIGN(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         Y(I)=Y(I)+DEL
         YPRIME(I)=YPRIME(I)+CJ*DEL
         CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0D0/DEL
         DO 220 L=1,NEQ
220      WM(NROW+L)=(E(L)-DELTA(L))*DELINV
      NROW=NROW+NEQ
      Y(I)=YSAVE
      YPRIME(I)=YPSAVE
210   CONTINUE
C
C
C     DO DENSE-MATRIX LU DECOMPOSITION ON PD
230      CALL BDGEFA(WM(NPD),NEQ,NEQ,IWM(LIPVT),IER)
      RETURN
C
C
C     DUMMY SECTION FOR IWM(MTYPE)=3
c-----------------------------------------------------------------------
  300 continue

      npde   = iwm(lnpde)
      kcol   = iwm(lkcol)
      nint   = iwm(lnint)

      neq1   = npde * (kcol * nint + nconti)
      neq2   = neq1 + npde * nint

      npdbk1 = npd + npde * npde * nconti
      npdbt1 = npdbk1 + npde * npde * nint * kcol * (kcol + nconti)
      npdtp2 = npdbt1 + npde * npde * nconti
      npdbk2 = npdtp2 + npde * npde * nconti
      npdbt2 = npdbk2 + npde * npde * nint * (kcol + 1) * (kcol + 1 +
     &                  nconti)

      lipvt2 = lipvt + neq1

      call jac(x,y,yprime,wm(npd),cj,rpar,ipar)

      call bcrdcmp(neq1, wm(npd), npde, 2*npde, wm(npdbk1), kcol*npde,
     &            (kcol+nconti)*npde, nint, wm(npdbt1), npde,
     &            iwm(lipvt), ier)

      call bcrdcmp(neq2, wm(npdtp2), npde, 2*npde, wm(npdbk2),
     &            (kcol+1)*npde, (kcol+1+nconti)*npde, nint,
     &             wm(npdbt2), npde, iwm(lipvt2), ier)

      return
c-----------------------------------------------------------------------
C
C     BANDED USER-SUPPLIED MATRIX
400   LENPD=(2*IWM(LML)+IWM(LMU)+1)*NEQ
      DO 410 I=1,LENPD
410      WM(NPDM1+I)=0.0D0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 550
C
C
C     BANDED FINITE-DIFFERENCE-GENERATED MATRIX
500   MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN(MBAND,NEQ)
      MEBAND=MBAND+IWM(LML)
      MEB1=MEBAND-1
      MSAVE=(NEQ/MBAND)+1
      ISAVE=NTEMP-1
      IPSAVE=ISAVE+MSAVE
      IRES=0
      SQUR=SQRT(UROUND)
      DO 540 J=1,MBA
         DO 510 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          WM(ISAVE+K)=Y(N)
          WM(IPSAVE+K)=YPRIME(N)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          Y(N)=Y(N)+DEL
510       YPRIME(N)=YPRIME(N)+CJ*DEL
      CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) RETURN
      DO 530 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0D0/DEL
          I1=MAX(1,(N-IWM(LMU)))
          I2=MIN(NEQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)+NPDM1
          DO 520 I=I1,I2
520         WM(II+I)=(E(I)-DELTA(I))*DELINV
530      CONTINUE
540   CONTINUE
C
C
C     DO LU DECOMPOSITION OF BANDED PD
550   CALL DGBFA(WM(NPD),MEBAND,NEQ,
     *    IWM(LML),IWM(LMU),IWM(LIPVT),IER)
      RETURN
C------END OF SUBROUTINE DDAJAC------
      END
