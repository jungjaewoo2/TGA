      DOUBLE PRECISION FUNCTION DDANRM (NEQ, V, WT, RPAR, IPAR)
C***BEGIN PROLOGUE  DDANRM
C***SUBSIDIARY
C***PURPOSE  Compute vector norm for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDANRM-S, DDANRM-D)
C***AUTHOR  PETZOLD, LINDA R., (LLNL)
c-----------------------------------------------------------------------
c
c Last modified by Rong Wang, Feb 22, 2002.
c
c-----------------------------------------------------------------------
C***DESCRIPTION
C-----------------------------------------------------------------------
C     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED
C     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH
C     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS
C     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.
C        DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
C-----------------------------------------------------------------------
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDANRM
C
      INTEGER  NEQ, IPAR(*)
      DOUBLE PRECISION  V(NEQ), WT(NEQ), RPAR(*)
C
      INTEGER  I
      DOUBLE PRECISION  SUM, VMAX
C
c-----------------------------------------------------------------------  
        double precision        temp
        integer                 itemp
        integer                 iwkdnm
        parameter              (iwkdnm = 49)
c                               rpar(ipar(iwkdnm)) is the work storage
c                               for the modification version of the
c                               subroutine DDANRM.
c
c-----------------------------------------------------------------------  
C***FIRST EXECUTABLE STATEMENT  DDANRM
      DDANRM = 0.0D0
      VMAX = 0.0D0
      itemp = ipar(iwkdnm) - 1
      DO 10 I = 1,NEQ
c-----------------------------------------------------------------------  
        itemp = itemp + 1
        rpar(itemp) = abs(v(i)/wt(i))
        if (rpar(itemp) .gt. vmax) vmax = rpar(itemp)
c-----------------------------------------------------------------------  
10      CONTINUE
      IF(VMAX .LE. 0.0D0) GO TO 30
      SUM = 0.0D0
c-----------------------------------------------------------------------  
      itemp = ipar(iwkdnm) - 1
c-----------------------------------------------------------------------  
      DO 20 I = 1,NEQ
c-----------------------------------------------------------------------  
        itemp = itemp + 1
        temp = rpar(itemp)/vmax
        sum = sum + temp * temp
20    continue
c-----------------------------------------------------------------------  
      DDANRM = VMAX*SQRT(SUM/NEQ)
30    CONTINUE
      RETURN
C------END OF FUNCTION DDANRM------
      END
