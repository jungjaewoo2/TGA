C*PGCONL -- label contour map of a 2D data array 
C%void cpgconl(const float *a, int idim, int jdim, int i1, int i2, \
C% int j1, int j2, const float * c, int nc, const float *tr, const int * cparam, \
C% int idxc, const char *label, float ts, int intval, int minint);
C+
      SUBROUTINE PGCONL (A, IDIM, JDIM, I1, I2, J1, J2, C, NC, TR, CPARAM,
     1                   IDXC, LABEL, TS, INTVAL, MININT)
      INTEGER IDIM, JDIM, I1, J1, I2, J2, INTVAL, NC, MININT, IDXC
      INTEGER CPARAM (*)
      REAL A(IDIM,JDIM), C(*), TR(6), TS, RTS
      CHARACTER*(*) LABEL
C
C Label a contour map drawn with routine PGCONT. Routine PGCONT should
C be called first to draw the contour lines, then this routine should be
C called to add the labels. Labels are written at intervals along the
C contour lines, centered on the contour lines with lettering aligned
C in the up-hill direction. Labels are opaque, so a part of the under-
C lying contour line is obscured by the label. Labels use the current
C attributes (character height, line width, color index, character
C font).
C
C The first 9 arguments are the same as those supplied to PGCONT, and
C should normally be identical to those used with PGCONT. Note that
C only one contour level can be specified; tolabel more contours, call
C PGCONL for each level.
C
C The Label is supplied as a character string in argument LABEL.
C
C The spacing of labels along the contour is specified by parameters
C INTVAL and MININT. The routine follows the contour through the
C array, counting the number of cells that the contour crosses. The
C first label will be written in the MININT'th cell, and additional
C labels will be written every INTVAL cells thereafter. A contour
C that crosses less than MININT cells will not be labelled. Some
C experimentation may be needed to get satisfactory results; a good
C place to start is INTVAL=20, MININT=10.
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
C  IDXC   (input) : index in the range 1:NC of the contour level wanted
C                   to be plotted C(IDXC)
C  LABEL  (input) : character string to be used to label the specified
C                   contour. Leading and trailing blank spaces are
C                   ignored.
C  TS     (input) : float size of text used for the labeling
C  INTVAL (input) : spacing along the contour between labels, in
C                   grid cells.
C  MININT (input) : contours that cross less than MININT cells
C                   will not be labelled.
C--
C  5-May-1994 - New routine; this routine is virtually identical to
C               PGCONT, but calls PGCONX with a different external
C               routine [TJP].
C  4-Feb-1997 - PGCONX requires an array argument, not scalar [TJP].
C-----------------------------------------------------------------------
      INCLUDE  'pgplot.inc'
      INTEGER  I,ILT,ILW,ILC,MININD
      LOGICAL  PGNOTO
      EXTERNAL PGCL
C
      IF (PGNOTO('PGCONL')) RETURN
C
C Save TRANS matrix and other parameters.
C
      DO 10 I=1,6
          TRANS(I) = TR(I)
   10 CONTINUE
C     C SAVE OLD VALUES FOR COLOR,LINEWIDTH,LINETYPE and CHARACTER HEIGHT
      CALL GRQCI(ILC)
      CALL GRQLW(ILW)
      CALL GRQLS(ILT)
      CALL PGQCH(RTS)
C     C SET PROVIDED VALUES FOR LINEWIDTH,LINETYPE
      CALL GRSLS(CPARAM(2))
      CALL GRSLW(CPARAM(3))
      CALL PGSCH(TS)
C     C specify label parameters:
      PGCINT = INTVAL
      PGCMIN = MININT
      PGCLAB = LABEL
C
C Use PGCONX with external function PGCL.
C
      IF (CPARAM(1).EQ.1) THEN
C         C case 1
C         C   single color 
          CALL GRSCI(CPARAM(4))
          CALL PGCONX(A, IDIM, JDIM, I1, I2, J1, J2, C(IDXC), 1, PGCL)
      ELSE
C         C multi-color contour lines
          MININD=PGMNCI(PGID)
          CALL GRQCR(MININD, CRM, CGM, CBM)
          ZMIN = C(ABS(NC)+1)
          ZMAX = C(ABS(NC)+2)
          IF (CPARAM(1).EQ.2) THEN
C             C case 2
C             C   gray map mapping
              CALL GRQCR(CPARAM(4), CR0, CG0, CB0)
              CALL GRQCR(CPARAM(5), CR1, CG1, CB1)
              A0 = (C(IDXC) - ZMIN) / (ZMAX-ZMIN)
              IF (A0.LT.0.0) A0 = 0.0
              IF (A0.GT.1.0) A0 = 1.0
              A1 = 1 - A0
              R = A0*CR0+A1*CR1
              G = A0*CG0+A1*CG1
              B = A0*CB0+A1*CB1
              CALL GRSCR (MININD, R, G, B)
              CALL GRSCI (MININD)
              CALL PGCONX(A,IDIM,JDIM,I1,I2,J1,J2,C(IDXC),1,PGCL)
          ELSE IF (CPARAM(1).EQ.3) THEN
C             C case 3
C             C   palette mapping
              A0 = (C(IDXC) - ZMIN) / (ZMAX-ZMIN)
              IF (A0.LT.0.0) A0 = 0.0
              IF (A0.GT.1.0) A0 = 1.0
              A1 = 1 - A0
C             C red gnuplot function
              IF (CPARAM(4).GT.0) THEN
                  R = RGB1CH(A0, CPARAM(4),.TRUE.)
              ELSE
                  R = RGB1CH(A1,-CPARAM(4),.TRUE.)
              ENDIF
C             C green gnuplot function
              IF (CPARAM(5).GT.0) THEN
                  G = RGB1CH(A0, CPARAM(5),.TRUE.)
              ELSE
                  G = RGB1CH(A1,-CPARAM(5),.TRUE.)
              ENDIF
C             C blue gnuplot function
              IF (CPARAM(6).GT.0) THEN
                  B = RGB1CH(A0, CPARAM(6),.TRUE.)
              ELSE
                  B = RGB1CH(A1,-CPARAM(6),.TRUE.)
              ENDIF
              CALL GRSCR (MININD, R, G, B)
              CALL GRSCI (MININD)
              CALL PGCONX(A,IDIM,JDIM,I1,I2,J1,J2,C(IDXC),-1,PGCL)
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
      CALL PGSCH(RTS)
      END
