C*PGRGBFOR-- palette-map of a 2D data array using GNUPLOT RGB Formulae
C%void cpgrgbfor(const float *a, int idim, int jdim, int i1, int i2, \
C% int j1, int j2, float fg, float bg, const float *tr, int * frgb);
C+
      SUBROUTINE PGRGBFOR (A, IDIM, JDIM, I1, I2, J1, J2,
     1                   FG, BG, TR, FRGB)
      INTEGER IDIM, JDIM, I1, I2, J1, J2, FRGB(*)
      REAL    A(IDIM,JDIM), FG, BG, TR(6)
C
C Draw gray-scale map of an array in current window. The subsection
C of the array A defined by indices (I1:I2, J1:J2) is mapped onto
C the view surface world-coordinate system by the transformation
C matrix TR. The resulting quadrilateral region is clipped at the edge
C of the window and shaded with the shade at each point determined
C by the corresponding array value.  The shade is a number in the
C range 0 to 1 obtained by linear interpolation between the background
C level (BG) and the foreground level (FG), i.e.,
C
C   shade = [A(i,j) - BG] / [FG - BG]
C
C The background level BG can be either less than or greater than the
C foreground level FG.  Points in the array that are outside the range
C BG to FG are assigned shade 0 or 1 as appropriate.
C
C PGGRAY uses two different algorithms, depending how many color
C indices are available in the color index range specified for images.
C (This range is set with routine PGSCIR, and the current or default
C range can be queried by calling routine PGQCIR).
C
C If 16 or more color indices are available, PGGRAY first assigns
C color representations to these color indices to give a linear ramp
C between the background color (color index 0) and the foreground color
C (color index 1), and then calls PGIMAG to draw the image using these
C color indices. In this mode, the shaded region is "opaque": every
C pixel is assigned a color.
C
C If less than 16 color indices are available, PGGRAY uses only
C color index 1, and uses  a "dithering" algorithm to fill in pixels,
C with the shade (computed as above) determining the faction of pixels
C that are filled. In this mode the shaded region is "transparent" and
C allows previously-drawn graphics to show through.
C
C The transformation matrix TR is used to calculate the world
C coordinates of the center of the "cell" that represents each
C array element. The world coordinates of the center of the cell
C corresponding to array element A(I,J) are given by:
C
C          X = TR(1) + TR(2)*I + TR(3)*J
C          Y = TR(4) + TR(5)*I + TR(6)*J
C
C Usually TR(3) and TR(5) are zero -- unless the coordinate
C transformation involves a rotation or shear.  The corners of the
C quadrilateral region that is shaded by PGGRAY are given by
C applying this transformation to (I1-0.5,J1-0.5), (I2+0.5, J2+0.5).
C
C Arguments:
C  A      (input)  : the array to be plotted.
C  IDIM   (input)  : the first dimension of array A.
C  JDIM   (input)  : the second dimension of array A.
C  I1, I2 (input)  : the inclusive range of the first index
C                    (I) to be plotted.
C  J1, J2 (input)  : the inclusive range of the second
C                    index (J) to be plotted.
C  FG     (input)  : the array value which is to appear with the
C                    foreground color (corresponding to color index 1).
C  BG     (input)  : the array value which is to appear with the
C                    background color (corresponding to color index 0).
C  TR     (input)  : transformation matrix between array grid and
C                    world coordinates.
C  FRGB   (input)  : transformation formula between gray shade x \in [0,1]
C                    and corresponding RGB color.
C--
C  2-Sep-1987: remove device-dependent code to routine GRGRAY (TJP).
C  7-Jun-1988: change documentation and argument names (TJP).
C 31-May-1989: allow 1-pixel wide arrays to be plotted (TJP).
C 17-Mar-1994: pass PG scaling info to lower routines (TJP).
C 15-Sep-1994: use PGITF attribute (TJP).
C  8-Feb-1995: use color ramp based on current foreground and background
C              colors (TJP).
C  6-May-1996: allow multiple devives (TJP).
C 18-Nov-2017: add Gnuplot-like RGB Formula support for palette maping (pm3d)
C              part of PGPLOT-rlabpatch revival project
C-----------------------------------------------------------------------
      INCLUDE 'pgplot.inc'
      INCLUDE 'grpckg1.inc'
      REAL A0,A1,CR,CG,CB,CR0,CG0,CB0,CR1,CG1,CB1
      REAL PA(6)
      LOGICAL PGNOTO
      INTEGER I,MININD,MAXIND
      INTRINSIC REAL
C
C Check inputs.
C
!       WRITE (*,*) 'BG=',BG
!       WRITE (*,*) 'FG=',FG
      IF (PGNOTO('PGRGBFOR')) RETURN
      IF (I1.LT.1 .OR. I2.GT.IDIM .OR. I1.GT.I2 .OR.
     1    J1.LT.1 .OR. J2.GT.JDIM .OR. J1.GT.J2) THEN
          CALL GRWARN('PGRGBFOR: invalid range I1:I2, J1:J2')
      ELSE IF (FG.EQ.BG) THEN
          CALL GRWARN('PGRGBFOR: foreground level = background level')
      ELSE
C
C Call lower-level routine to do the work.
C
          CALL PGBBUF
          CALL PGSAVE
          CALL PGSCI(1)
          PA(1) = TR(1)*PGXSCL(PGID) + PGXORG(PGID)
          PA(2) = TR(2)*PGXSCL(PGID)
          PA(3) = TR(3)*PGXSCL(PGID)
          PA(4) = TR(4)*PGYSCL(PGID) + PGYORG(PGID)
          PA(5) = TR(5)*PGYSCL(PGID)
          PA(6) = TR(6)*PGYSCL(PGID)
          MININD=PGMNCI(PGID)
          MAXIND=PGMXCI(PGID)
C
          DO 5 I=MININD,MAXIND
              A0 = REAL(I-MININD)/REAL(MAXIND-MININD)
              CALL RGBFORMULA(A0,FRGB,CR,CG,CB)
!               WRITE (*,*) 'A0=',A0
!               WRITE (*,*) 'CR,CG,CB=',CR,CG,CB
              CALL GRSCR(I, CR, CG, CB)
 5        CONTINUE
          CALL GRIMG0(A, IDIM, JDIM, I1, I2, J1, J2,
     :                FG,BG,PA,MININD,MAXIND,PGITF(PGID))
          CALL PGEBUF
          CALL PGUNSA
      END IF
C-----------------------------------------------------------------------
      END


      SUBROUTINE RGBFORMULA (A0, FRGB, R, G, B)
      INTEGER FRGB(*)
      REAL    A0, R, G, B, DA
C-----------------------------------------------------------------------
      REAL    A1, CR0, CG0, CB0, CR1, CG1, CB1
      INTRINSIC REAL
C
C
C
      IF (FRGB(1).EQ.0) THEN
C         C INPUT FRGB(3) = (0, BACKGROUND COLOR, FOREGROUND COLOR INDEX)
          A1 = 1.0 - A0
          CALL GRQCR(FRGB(2), CR0, CG0, CB0)
          CALL GRQCR(FRGB(3), CR1, CG1, CB1)
          R = A1*CR0+A0*CR1
          G = A1*CG0+A0*CG1
          B = A1*CB0+A0*CB1
          RETURN
      ELSE IF (FRGB(1).EQ.1) THEN
C         C INPUT FRGB(4) = (1, IDXFUNCR, IDXFUNCG, IDXFUNCB)
          A1 = 1.0 - A0
C         R:
          IF (FRGB(2).GT.0) THEN
              R = RGB1CH(A1, FRGB(2),.TRUE.)
          ELSE
              R = RGB1CH(A0,-FRGB(2),.TRUE.)
          ENDIF
C         G:
          IF (FRGB(3).GT.0) THEN
              G = RGB1CH(A1, FRGB(3),.TRUE.)
          ELSE
              G = RGB1CH(A0,-FRGB(3),.TRUE.)
          ENDIF
C         B:
          IF (FRGB(4).GT.0) THEN
              B = RGB1CH(A1, FRGB(4),.TRUE.)
          ELSE
              B = RGB1CH(A0,-FRGB(4),.TRUE.)
          ENDIF
          RETURN
      ELSE IF (FRGB(1).EQ.2) THEN
C         C INPUT FRGB(4) = (2, NUMCOLTABLE, FIRSTCOLOR)
          IF (A0.GE.1.0) THEN
            CALL GRQCR(FRGB(2)+FRGB(3)-1, R, G, B)
          ELSE IF (A0.LE.0.0) THEN
            CALL GRQCR(FRGB(3), R, G, B)
          ELSE
            DDX = A0 * REAL(FRGB(2))
            IDX = INT(DDX)
            IF (REAL(IDX).EQ.DDX) THEN
                CALL GRQCR(FRGB(3)+IDX, R, G, B)
            ELSE
                CALL GRQCR(FRGB(3)+IDX,   CR0, CG0, CB0)
                CALL GRQCR(FRGB(3)+IDX+1, CR1, CG1, CB1)
                A0 = (DDX - REAL(IDX))
                A1 = 1.0 - A0
                R = A0*CR0+A1*CR1
                G = A0*CG0+A1*CG1
                B = A0*CB0+A1*CB1
            ENDIF
          ENDIF
      ENDIF
      RETURN
      END


      REAL FUNCTION RGB1CH(X,N,NOTMOD)
C
C GNUPLOT palette function for individual channels
C as reported on 
C   http://gnuplot.cvs.sourceforge.net/viewvc/gnuplot/gnuplot/src/getcolor.c?revision=1.39&content-type=text%2Fplain
      INTEGER N
      LOGICAL NOTMOD
      REAL    X,PI
C-----------------------------------------------------------------------
      PARAMETER (PI=3.14159)
      IF (N.EQ.1) THEN
          RGB1CH = 0.5
      ELSE IF (N.EQ.2) THEN
          RGB1CH = 1.0
      ELSE IF (N.EQ.3) THEN
          RGB1CH = X
      ELSE IF (N.EQ.4) THEN
          RGB1CH = X * X
      ELSE IF (N.EQ.5) THEN
          RGB1CH = X * X * X
      ELSE IF (N.EQ.6) THEN
          RGB1CH = X * X * X * X
      ELSE IF (N.EQ.7) THEN
          RGB1CH = SQRT(X)
      ELSE IF (N.EQ.8) THEN
          RGB1CH = SQRT(SQRT(X))
      ELSE IF (N.EQ.9) THEN
          RGB1CH = SIN(0.5 * PI * X)
      ELSE IF (N.EQ.10) THEN
          RGB1CH = COS(0.5 * PI * X)
      ELSE IF (N.EQ.11) THEN
          RGB1CH = ABS(X - 0.5)
      ELSE IF (N.EQ.12) THEN
          RGB1CH = (2.0*X-1.0)*(2.0*X-1.0)
      ELSE IF (N.EQ.13) THEN
          RGB1CH = SIN(PI * X)
      ELSE IF (N.EQ.14) THEN
          RGB1CH = ABS(COS(PI*X))
      ELSE IF (N.EQ.15) THEN
          RGB1CH = SIN(2.0*PI*X)
      ELSE IF (N.EQ.16) THEN
          RGB1CH = COS(2.0*PI*X)
      ELSE IF (N.EQ.17) THEN
          RGB1CH = ABS(SIN(2.0*PI*X))
      ELSE IF (N.EQ.18) THEN
          RGB1CH = ABS(COS(2.0*PI*X))
      ELSE IF (N.EQ.19) THEN
          RGB1CH = ABS(SIN(4.0*PI*X))
      ELSE IF (N.EQ.20) THEN
          RGB1CH = ABS(COS(4.0*PI*X))
      ELSE IF (N.EQ.21) THEN
          RGB1CH = 3.0*X
      ELSE IF (N.EQ.22) THEN
          RGB1CH = 3.0*X-1.0
      ELSE IF (N.EQ.23) THEN
          RGB1CH = 3.0*X-2.0
      ELSE IF (N.EQ.24) THEN
          RGB1CH = ABS(3.0*X-1.0)
      ELSE IF (N.EQ.25) THEN
          RGB1CH = ABS(3.0*X-2.0)
      ELSE IF (N.EQ.26) THEN
          RGB1CH = ABS(1.5*X-0.5)
      ELSE IF (N.EQ.27) THEN
          RGB1CH = 1.5*X-1.0
      ELSE IF (N.EQ.28) THEN
          RGB1CH = ABS(1.5*X-0.5)
      ELSE IF (N.EQ.29) THEN
          RGB1CH = ABS(1.5*X-1.0)
      ELSE IF (N.EQ.30) THEN
          IF (X.LE.0.25) THEN
              RGB1CH = 0.0
          ELSE IF (X.GE.0.57) THEN
              RGB1CH = 1.0
          ELSE
              RGB1CH=X/0.32-0.78125
          ENDIF
      ELSE IF (N.EQ.31) THEN
          IF (X.LE.0.42) THEN
              RGB1CH = 0.0
          ELSE IF (X.GE.0.92) THEN
              RGB1CH = 1.0
          ELSE
              RGB1CH=2.0*X-0.84
          ENDIF
      ELSE IF (N.EQ.32) THEN
          IF (X.LE.0.42) THEN
              RGB1CH = 4.0*X
          ELSE IF (X.LE.0.92) THEN
              RGB1CH=1.84-2.0*X
          ELSE
              RGB1CH=X/0.08-11.5
          ENDIF
      ELSE IF (N.EQ.33) THEN
          RGB1CH = ABS(2.0*X-0.5)
      ELSE IF (N.EQ.34) THEN
          RGB1CH = 2.0*X
      ELSE IF (N.EQ.35) THEN
          RGB1CH = 2.0*X-0.5
      ELSE IF (N.EQ.36) THEN
          RGB1CH = 2.0*X-1.0
      ELSE
C         C I'm tempted to write some nasty message to the sorry
C         C luser who got this far:
C         C "why don't I reboot your computer while you consult the user manual?"
          RGB1CH = 0.0
      ENDIF
C
C     C i was tempted to put mod operation rather then clipping
C
      IF (NOTMOD) THEN
          IF (RGB1CH.GT.1.0) THEN
              RGB1CH=1.0
          ELSE IF (RGB1CH.LT.0.0) THEN
              RGB1CH=0.0
          ENDIF
      ELSE
101       CONTINUE
          IF (RGB1CH.GT.1.0) THEN
              RGB1CH = RGB1CH - 1.0
              GOTO 101
          ENDIF
102       CONTINUE
          IF (RGB1CH.LT.0.0) THEN
              RGB1CH = RGB1CH + 1.0
              GOTO 102
          ENDIF
      ENDIF
      RETURN
      END












