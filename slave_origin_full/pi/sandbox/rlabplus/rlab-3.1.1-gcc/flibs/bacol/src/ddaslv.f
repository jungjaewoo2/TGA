      SUBROUTINE DDASLV (NEQ, DELTA, WM, IWM, CJ)
c-----------------------------------------------------------------------
c     This routine has been modified to accept an almost block diagonal
c     (ABD) Jacobian matrix compatible with the linear algebra package,
c     COLROW.
c
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, September 4, 2001.
c
c-----------------------------------------------------------------------
C***BEGIN PROLOGUE  DDASLV
C***SUBSIDIARY
C***PURPOSE  Linear system solver for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDASLV-S, DDASLV-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------------
C     THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR
C     SYSTEM ARISING IN THE NEWTON ITERATION.
C     MATRICES AND REAL TEMPORARY STORAGE AND
C     REAL INFORMATION ARE STORED IN THE ARRAY WM.
C     INTEGER MATRIX INFORMATION IS STORED IN
C     THE ARRAY IWM.
C     FOR A DENSE MATRIX, THE LINPACK ROUTINE
C     BDGESL IS CALLED.
C     FOR A BANDED MATRIX,THE LINPACK ROUTINE
C     DGBSL IS CALLED.
C-----------------------------------------------------------------------
C***ROUTINES CALLED  DGBSL, BDGESL
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDASLV
C
      INTEGER  NEQ, IWM(*)
      DOUBLE PRECISION  DELTA(*), WM(*)
C
      EXTERNAL  DGBSL, BDGESL
C
      INTEGER  LIPVT, LML, LMU, LMTYPE, MEBAND, MTYPE, NPD
c-----------------------------------------------------------------------
      double precision cj

      integer nconti
      parameter (nconti = 2)

      integer lnpde, lkcol, lnint
      integer npde, kcol, nint
      integer neq1
      integer npdbk1, npdbt1
      integer npdtp2, npdbk2, npdbt2
      integer lipvt2

c-----------------------------------------------------------------------
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
C***FIRST EXECUTABLE STATEMENT  DDASLV
      MTYPE=IWM(LMTYPE)
      GO TO(100,100,300,400,400),MTYPE
C
C     DENSE MATRIX
100   CALL BDGESL(WM(NPD),NEQ,NEQ,IWM(LIPVT),DELTA,0)
      RETURN
C
C     DUMMY SECTION FOR MTYPE=3
300   CONTINUE
c-----------------------------------------------------------------------
      npde   = iwm(lnpde)
      kcol   = iwm(lkcol)
      nint   = iwm(lnint)

      neq1   = npde * (kcol * nint + nconti)

      npdbk1 = npd + npde * npde * nconti
      npdbt1 = npdbk1 + npde * npde * nint * kcol * (kcol + nconti)
      npdtp2 = npdbt1 + npde * npde * nconti
      npdbk2 = npdtp2 + npde * npde * nconti
      npdbt2 = npdbk2 + npde * npde * nint * (kcol + 1) * (kcol + 1 +
     &                  nconti)

      lipvt2 = lipvt + neq1

      call dscal(npde, cj, delta, 1)
      call dscal(npde, cj, delta(neq1-npde+1), 1)

      call bcrslve(wm(npd), npde, 2*npde, wm(npdbk1), kcol*npde,
     &            (kcol+nconti)*npde, nint, wm(npdbt1), npde,
     &            iwm(lipvt), delta, 0)

      call dscal(npde, cj, delta(neq1+1), 1)
      call dscal(npde, cj, delta(neq-npde+1), 1)

      call bcrslve(wm(npdtp2), npde, 2*npde, wm(npdbk2),
     &            (kcol+1)*npde, (kcol+1+nconti)*npde, nint,
     &            wm(npdbt2), npde, iwm(lipvt2), delta(neq1+1), 0)

c-----------------------------------------------------------------------
      RETURN
C
C     BANDED MATRIX
400   MEBAND=2*IWM(LML)+IWM(LMU)+1
      CALL DGBSL(WM(NPD),MEBAND,NEQ,IWM(LML),
     *  IWM(LMU),IWM(LIPVT),DELTA,0)
      RETURN
C------END OF SUBROUTINE DDASLV------
      END
