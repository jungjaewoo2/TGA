C*PGCONT -- contour map of a 2D data array (contour-following)
C%void cpglgcont(const float x1, const float y1, const float x2, const float y2, \
C% const float zmin, const float zmax, int nz, const int * cparam);
C+
      SUBROUTINE PGLGCONT (X1,Y1,X2,Y2,NZ,CPARAM)
      REAL X1,Y1,X2,Y2,C(*)
      INTEGER NZ,CPARAM(*)
C
C Draw a contour map of an array.  The map is truncated if
C necessary at the boundaries of the viewport.  Each contour line
C is drawn with the current line attributes (color index, style, and
C width); except that if argument NZ is positive (see below), the line
C style is set by PGCONT to 1 (solid) for positive contours or 2
C (dashed) for negative contours.
C
C Arguments:
C  X1,Y1  (input) : lower left  corner of the area where to plot contour lines
C  X2,Y2  (input) : upper right corner of the area where to plot contour lines
C  ZMIN,ZMAX (input) : min and max levels of contours
C  NZ     (input) : number of contour levels.
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
C (24-Nov-2017) Added input array to better specify the contour lines
C       coloring. [MK, part of rlabplus]
C-----------------------------------------------------------------------
      INCLUDE  'pgplot.inc'
      INTEGER  I,ILT,ILW,ILC,MININD
      REAL     DX
      REAL     CRM,CGM,CBM,CR0,CG0,CB0,CR1,CG1,CB1,R,G,B
      REAL     A0,A1
      LOGICAL  PGNOTO
      EXTERNAL PGCP
C
      IF (PGNOTO('PGLGCONT')) RETURN
C
C Save what we are about to change.
C
      CALL GRQCI(ILC)
      CALL GRQLW(ILW)
      CALL GRQLS(ILT)
C
      CALL GRSLS(CPARAM(2))
      CALL GRSLW(CPARAM(3))
C
      IF (CPARAM(1).EQ.1) THEN
C         C case 1
C         C   single color 
          CALL GRSCI(CPARAM(4))
          CALL GRMOVA(X1,0.5*(Y1+Y2))
          CALL GRLINA(X2,0.5*(Y1+Y2))
      ELSE
C         C multi-color contour lines
          DX = (X2 - X1)/(REAL(ABS(NZ))-1.0)
C
          MININD=PGMNCI(PGID)
          CALL GRQCR(MININD, CRM, CGM, CBM)
C
          IF (CPARAM(1).EQ.2) THEN
C             C case 2
C             C   grey mapping
              CALL GRQCR(CPARAM(4), CR0, CG0, CB0)
              CALL GRQCR(CPARAM(5), CR1, CG1, CB1)
              DO 11 I=1,ABS(NZ)
                  A0 = (REAL(I) - 1.0) / REAL(ABS(NZ))
                  A1 = 1 - A0
                  R = A0*CR0+A1*CR1
                  G = A0*CG0+A1*CG1
                  B = A0*CB0+A1*CB1
                  CALL GRSCR(MININD, R, G, B)
                  CALL GRSCI (MININD)
                  CALL GRMOVA(X1 + DX * (REAL(I) - 1.0),Y1)
                  CALL GRLINA(X1 + DX * (REAL(I) - 1.0),Y2)
   11         CONTINUE
          ELSE IF (CPARAM(1).EQ.3) THEN
C             C case 3
C             C   palette mapping
              DO 12 I=1,ABS(NZ)
                  A0 = (REAL(I) - 1.0) / REAL(ABS(NZ))
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
                  CALL GRSCR(MININD, R, G, B)
                  CALL GRSCI (MININD)
                  CALL GRMOVA(X1 + DX * (REAL(I) - 1.0),Y1)
                  CALL GRLINA(X1 + DX * (REAL(I) - 1.0),Y2)
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
