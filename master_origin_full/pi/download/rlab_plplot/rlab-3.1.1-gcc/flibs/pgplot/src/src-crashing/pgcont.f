C*PGCONT -- contour map of a 2D data array (contour-following)
C%void cpgcont(const float *a, int idim, int jdim, int i1, int i2, int j1, int j2, \
C% const float *c, int nc, const float *tr, const int * cparam);
C+
      SUBROUTINE PGCONT (A, IDIM, JDIM, I1, I2, J1, J2, C, NC, TR, CPARAM)
      INTEGER IDIM, JDIM, I1, J1, I2, J2, NC, CPARAM(*)
      REAL A(IDIM,JDIM), C(*), TR(6)
C
C Draw a contour map of an array.  The map is truncated if
C necessary at the boundaries of the viewport.  Each contour line
C is drawn with the current line attributes (color index, style, and
C width); except that if argument NC is positive (see below), the line
C style is set by PGCONT to 1 (solid) for positive contours or 2
C (dashed) for negative contours.
C
C Arguments:
C  A      (input) : data array.
C  IDIM   (input) : first dimension of A.
C  JDIM   (input) : second dimension of A.
C  I1, I2 (input) : range of first index to be contoured (inclusive).
C  J1, J2 (input) : range of second index to be contoured (inclusive).
C  C      (input) : array of size NC+2 where 1:NC are contour levels, while
C                   C(NC+1) is ZMIN and C(NC+2) is ZMAX for coloring the
C                   levels
C  NC     (input) : +/- number of contour levels (less than or equal
C                   to dimension of C). If NC is positive, it is the
C                   number of contour levels, and the line-style is
C                   chosen automatically as described above. If NC is
C                   negative, it is minus the number of contour
C                   levels, and the current setting of line-style is
C                   used for all the contours.
C  TR     (input) : array defining a transformation between the I,J
C                   grid of the array and the world coordinates.
C                   The world coordinates of the array point A(I,J)
C                   are given by:
C                     X = TR(1) + TR(2)*I + TR(3)*J
C                     Y = TR(4) + TR(5)*I + TR(6)*J
C                   Usually TR(3) and TR(5) are zero - unless the
C                   coordinate transformation involves a rotation or
C                   shear.
C  CPARAM(6) (input) : Defines how multiple contour lines are plotted
C                     CPARAM(2) - line type
C                     CPARAM(3) - line width
C                   Then plot all the contour lines
C                   CPARAM(1) = 1: using single line with properties
C                     CPARAM(4) - line color
C                   CPARAM(1) = 2: use grey mapping between
C                     CPARAM(4) - first/background colur
C                     CPARAM(5) - last/foreground colur
C                   CPARAM(1) = 3: use palette mapping with indices in
C                     CPARAM(4) - red color function
C                     CPARAM(5) - green color function
C                     CPARAM(6) - blue color function
C--
C (7-Feb-1983)
C (24-Aug-1984) Revised to add the option of not automatically
C       setting the line-style. Sorry about the ugly way this is
C       done (negative NC); this is the least incompatible way of doing
C       it (TJP).
C (21-Sep-1989) Changed to call PGCONX instead of duplicating the code
C       [TJP].
C (24-Nov-2017) Added input array to better specify the contour lines
C       coloring. [MK, part of rlabplus]
C-----------------------------------------------------------------------
      INCLUDE  'pgplot.inc'
      INTEGER  I,ILT,ILW,ILC,MININD
      REAL     ZMIN,ZMAX
      REAL     CRM,CGM,CBM,CR0,CG0,CB0,CR1,CG1,CB1,R,G,B
      REAL     A0,A1
      LOGICAL  PGNOTO
      EXTERNAL PGCP
C
      IF (PGNOTO('PGCONT')) RETURN
C
C Save TRANS matrix.
C
      DO 10 I=1,6
          TRANS(I) = TR(I)
   10 CONTINUE
C     C SAVE OLD VALUES FOR COLOR,LINEWIDTH,LINETYPE
      CALL GRQCI(ILC)
      CALL GRQLW(ILW)
      CALL GRQLS(ILT)
C     C SET PROVIDED VALUES FOR LINEWIDTH,LINETYPE
      CALL GRSLS(CPARAM(2))
      CALL GRSLW(CPARAM(3))
      IF (CPARAM(1).EQ.1) THEN
C         C case 1
C         C   single color 
          CALL GRSCI(CPARAM(4))
          CALL PGCONX (A, IDIM, JDIM, I1, I2, J1, J2, C, NC, PGCP)
      ELSE
C         C multi-color contour lines
          MININD=PGMNCI(PGID)
          CALL GRQCR(MININD, CRM, CGM, CBM)
          ZMIN = C(ABS(NC)+1)
          ZMAX = C(ABS(NC)+2)
          IF (CPARAM(1).EQ.2) THEN
C             C case 2
C             C   grey mapping
!               WRITE (*,*) 'PGCONT.F: GRAY MAPPING'
              CALL GRQCR(CPARAM(4), CR0, CG0, CB0)
              CALL GRQCR(CPARAM(5), CR1, CG1, CB1)
              DO 11 I=1,ABS(NC)
                  A0 = (C(I) - ZMIN) / (ZMAX-ZMIN)
                  IF (A0.LT.0.0) A0 = 0.0
                  IF (A0.GT.1.0) A0 = 1.0
!                   WRITE (*,*) 'A0 =',A0
                  A1 = 1 - A0
                  R = A0*CR0+A1*CR1
                  G = A0*CG0+A1*CG1
                  B = A0*CB0+A1*CB1
                  CALL GRSCR (MININD, R, G, B)
                  CALL GRSCI (MININD)
                  CALL PGCONX(A,IDIM,JDIM,I1,I2,J1,J2,C(I),-1,PGCP)
C
!                   WRITE(*,*) 'NLAB,LABEL=',NLAB,LABELS
!                   IF (NLAB.GT.0) THEN
!                       PGCINT = 1
!                       PGCMIN = 1
!                       PGCLAB = LABELS
! C                     C
! C                     C Use PGCONX with external function PGCL.
! C                     C
!                       CALL PGCONX (A,IDIM,JDIM,I1,I2,J1,J2,
!      1                             C(I),1,PGCL)
!                   ENDIF

   11         CONTINUE
          ELSE IF (CPARAM(1).EQ.3) THEN
C             C case 3
C             C   palette mapping
              DO 12 I=1,ABS(NC)
                  A0 = (C(I) - ZMIN) / (ZMAX-ZMIN)
                  IF (A0.LT.0.0) A0 = 0.0
                  IF (A0.GT.1.0) A0 = 1.0
                  A1 = 1 - A0
C                 C red gnuplot function
                  IF (CPARAM(4).GT.0) THEN
                      R = RGB1CH(A0, CPARAM(4),.TRUE.)
                  ELSE
                      R = RGB1CH(A1,-CPARAM(4),.TRUE.)
                  ENDIF
C                 C green gnuplot function
                  IF (CPARAM(5).GT.0) THEN
                      G = RGB1CH(A0, CPARAM(5),.TRUE.)
                  ELSE
                      G = RGB1CH(A1,-CPARAM(5),.TRUE.)
                  ENDIF
C                 C blue gnuplot function
                  IF (CPARAM(6).GT.0) THEN
                      B = RGB1CH(A0, CPARAM(6),.TRUE.)
                  ELSE
                      B = RGB1CH(A1,-CPARAM(6),.TRUE.)
                  ENDIF
                  CALL GRSCR (MININD, R, G, B)
                  CALL GRSCI (MININD)
                  CALL PGCONX(A,IDIM,JDIM,I1,I2,J1,J2,C(I),-1,PGCP)
   12         CONTINUE
          ENDIF
          CALL GRSCR(MININD, CRM, CGM, CBM)
      ENDIF
C
C Use PGCONX with external function PGCP, which applies the TRANS
C scaling.
C
C
      CALL GRSCI(ILC)
      CALL GRSLW(ILW)
      CALL GRSLS(ILT)
      END
