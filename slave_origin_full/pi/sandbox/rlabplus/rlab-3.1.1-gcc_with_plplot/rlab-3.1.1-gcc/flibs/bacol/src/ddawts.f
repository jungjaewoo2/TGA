c-----------------------------------------------------------------------
      SUBROUTINE BDAWTS (IWT, RTOL, ATOL, Y, WT, RPAR, IPAR)
c-----------------------------------------------------------------------
C***BEGIN PROLOGUE  BDAWTS
C***SUBSIDIARY
C***PURPOSE  Set error weight vector for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDAWTS-S, BDAWTS-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------------
C     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR
C     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
C     I=1,-,N.
C     RTOL AND ATOL ARE SCALARS IF IWT = 0,
C     AND VECTORS IF IWT = 1.
C-----------------------------------------------------------------------
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  BDAWTS
c-----------------------------------------------------------------------
c   Last modified by Rong Wang, April 1, 2001
c-----------------------------------------------------------------------
C
      INTEGER       IWT, IPAR(*)
      DOUBLE PRECISION  RTOL(*), ATOL(*), Y(*), WT(*), RPAR(*)
c-----------------------------------------------------------------------
c     Constant:
      integer nconti
      parameter  (nconti =  2)
c
      integer npde, kcol, nint, npts, j, ij
      integer inpde, ikcol, inint

c-----------------------------------------------------------------------
C
      INTEGER  I
      DOUBLE PRECISION  ATOLI, RTOLI
C
c-----------------------------------------------------------------------
c     Direct IPAR indices:
      parameter  (inpde =  1)
      parameter  (ikcol =  2)
      parameter  (inint =  3)
c
c-----------------------------------------------------------------------
C***FIRST EXECUTABLE STATEMENT  BDAWTS
c-----------------------------------------------------------------------
      npde = ipar(inpde)
      nint = ipar(inint)
      kcol = ipar(ikcol)

      npts = nint * kcol + nconti + nint * (kcol + 1) + nconti

      if (iwt .eq. 0) then
         rtoli = rtol(1)
         atoli = atol(1)
         do 10 i = 1, npde * npts
            wt(i) = rtoli * abs(y(i)) + atoli
   10    continue
      else
         do 30 i = 1, npts
            do 20 j = 1, npde
               rtoli = rtol(j)
               atoli = atol(j)
               ij = (i - 1) * npde + j
               wt(ij) = rtoli * abs(y(ij)) + atoli
   20       continue
   30    continue
      endif

c-----------------------------------------------------------------------
      RETURN
C-----------END OF SUBROUTINE BDAWTS------------------------------------
      END
