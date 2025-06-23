        SUBROUTINE BCRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,JOB)
C
C***************************************************************
C
C  C R S L V E  SOLVES THE LINEAR SYSTEM
C                       A*X = B
C  USING THE DECOMPOSITION ALREADY GENERATED IN  C R D C M P.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY  ...
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         THE PIVOT VECTOR FROM  C R D C M P
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C               JOB    - INTEGER, INDICATING:
C                      = 0: SOLVE A*X = B;
C                      NON-ZERO: SOLVE TRANSPOSE(A)*X = B.
C
C       *** ON RETURN  ...
C
C                    B - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR
C
C***************************************************************
C
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,B
        DOUBLE PRECISION DOTPRD,BJ,XINCRJ,BINCRJ,SWAP,BI
        INTEGER NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(1),
     *          JOB
        DIMENSION TOPBLK(NRWTOP,1),ARRAY(NRWBLK,NCLBLK,1),
     *          BOTBLK(NRWBOT,1),B(1)
        INTEGER NRWTP1,NRWBK1,NVRLP1,NRWBT1,NROWEL,NVRLP0,NBLKS1,
     *          NBKTOP,J,I,LOOP,INCR,INCRJ,INCRI,JPIVOT,JRWTOP,
     *          LL,L1,IPLUSN,INCRN,NRWTP0,NRWEL1,K,INCRTP,NRWBTL,
     *          IPVTN,NRWELL,IPVTI,L
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        NRWTP1 = NRWTOP+1
        NRWBK1 = NRWBLK+1
        NVRLP1 = NOVRLP+1
        NRWTP0 = NRWTOP-1
        NRWBT1 = NRWBOT+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
        NBLKS1 = NBLOKS+1
        NBKTOP = NRWBLK+NRWTOP
C
C       IF JOB IS NON-ZERO, TRANSFER TO THE SECTION DEALING WITH
C       TRANSPOSE(A)*X = B.
C
        IF ( JOB .NE. 0 ) GO TO 530
C
C***************************************************************
C
C               ****  FORWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST, IN TOPBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 130 J = 1,NRWTOP
           B(J) = B(J)/TOPBLK(J,J)
           IF(J.EQ.NRWTOP)GO TO 120
              BJ = -B(J)
              LOOP = J+1
              DO 110 I = LOOP,NRWTOP
                 B(I) = B(I)+TOPBLK(I,J)*BJ
110           CONTINUE
120        CONTINUE
130     CONTINUE
C
C       ********************************************************
C
C          ***  IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCR = 0
        DO 280 K = 1,NBLOKS
           INCRTP = INCR+NRWTOP
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 220 J = 1,NRWTOP
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 210 I = 1,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
210           CONTINUE
220        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 240 J = NRWTP1,NRWBLK
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 225
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
225           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 230 I = LOOP,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BINCRJ
230           CONTINUE
240        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 270 J = NRWBK1,NBKTOP
              INCRJ = INCR+J
              JRWTOP = J -NRWTOP
              B(INCRJ) = B(INCRJ)/ARRAY(JRWTOP,J,K)
              IF(J.EQ.NBKTOP)GO TO 260
                 XINCRJ = -B(INCRJ)
                 LOOP = J-NRWTP0
                 DO 250 I = LOOP,NRWBLK
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
250              CONTINUE
260           CONTINUE
270        CONTINUE
           INCR = INCR+NRWBLK
280     CONTINUE
C
C       ********************************************************
C
C          ***  FINALLY, IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCRTP = INCR+NRWTOP
        DO 320 J = 1,NRWTOP
           INCRJ = INCR+J
           XINCRJ = -B(INCRJ)
           DO 310 I = 1,NRWBOT
              INCRI = INCRTP+I
              B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
310        CONTINUE
320     CONTINUE
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 350
           DO 340 J = NRWTP1,NVRLP0
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 325
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
325           CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 330 I = LOOP,NRWBOT
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*BINCRJ
330           CONTINUE
340        CONTINUE
350     CONTINUE
C
C***************************************************************
C
C               ****  BACKWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 430 LL = 1,NRWBOT
           J = NVRLP1-LL
           INCRJ = INCR+J
           NRWBTL = NRWBT1-LL
           B(INCRJ) = B(INCRJ)/BOTBLK(NRWBTL,J)
           IF(LL.EQ.NRWBOT)GO TO 420
              XINCRJ = -B(INCRJ)
              LOOP = NRWBOT-LL
              DO 410 I = 1,LOOP
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
410           CONTINUE
420        CONTINUE
430     CONTINUE
C
C       ********************************************************
C
C          ***  THEN IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 490 L = 1,NBLOKS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           K = NBLKS1-L
           INCR = INCR-NRWBLK
           DO 450 L1 = NRWEL1,NRWBLK
              I = NRWBLK+NRWEL1-L1
              IPLUSN = I+NRWTOP
              LOOP = IPLUSN+1
              INCRN = INCR+IPLUSN
              DOTPRD = B(INCRN)
              DO 440 J = LOOP,NCLBLK
                 INCRJ = INCR+J
                 DOTPRD = DOTPRD-ARRAY(I,J,K)*B(INCRJ)
440           CONTINUE
              B(INCRN) = DOTPRD
              IPVTN = PIVOT(INCRN)
              IF(INCRN.EQ.IPVTN)GO TO 445
                 SWAP = B(INCRN)
                 B(INCRN) = B(IPVTN)
                 B(IPVTN) = SWAP
445           CONTINUE
450        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           INCRTP = INCR+NRWTOP
           DO 460 J = NRWBK1,NCLBLK
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 455 I = 1,NROWEL
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
455           CONTINUE
460        CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 480 LL = 1,NROWEL
              J = NRWBK1-LL
              INCRJ = INCR+J
              NRWELL = NRWEL1-LL
              B(INCRJ) = B(INCRJ)/ARRAY(NRWELL,J,K)
              IF(LL.EQ.NROWEL)GO TO 470
                 XINCRJ = -B(INCRJ)
                 LOOP = NROWEL-LL
                 DO 465 I = 1,LOOP
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
465              CONTINUE
470           CONTINUE
480        CONTINUE
490     CONTINUE
C
C       ********************************************************
C
C          ***  IN TOPBLK FINISH WITH....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 520 L = 1,NRWTOP
           I = NRWTP1-L
           LOOP = I+1
           DOTPRD = B(I)
           DO 510 J = LOOP,NOVRLP
              DOTPRD = DOTPRD-TOPBLK(I,J)*B(J)
510        CONTINUE
           B(I) = DOTPRD
           IPVTI = PIVOT(I)
           IF(I.EQ.IPVTI)GO TO 515
                 SWAP = B(I)
                 B(I) = B(IPVTI)
                 B(IPVTI) = SWAP
515        CONTINUE
520     CONTINUE
C
C       RETURN FROM THE SOLUTION OF A.X = B.
        RETURN
C
C       IF JOB IS NON-ZERO, SOLVE TRANSPOSE(A)*X = B:
C
  530   CONTINUE

C       FIRST, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U).

        DO 540 I = 1,NRWTOP
           IPVTI = PIVOT(I)
           IF ( I .NE. IPVTI ) THEN
              SWAP = B(I)
              B(I) = B(IPVTI)
              B(IPVTI) = SWAP
           ENDIF
           BI = -B(I)
           LOOP = I+1
           DO 535 J = LOOP,NOVRLP
              B(J) = B(J) + BI*TOPBLK(I,J)
  535      CONTINUE
  540   CONTINUE

C       IN EACH BLOCK, K = 1,..,NBLOKS:

        INCR = NRWTOP
        DO 590 K = 1,NBLOKS

C          FIRST, THE FORWARD SOLUTION.

           DO 550 J = 1,NROWEL
              INCRJ = INCR + J
              DO 545 I = 1,J-1
                 B(INCRJ) = B(INCRJ) - ARRAY(I,NRWTOP+J,K)*B(INCR+I)
  545         CONTINUE
              B(INCRJ) = B(INCRJ)/ARRAY(J,NRWTOP+J,K)
  550       CONTINUE

C           FORWARD MODIFICATION.

            DO 570 I = 1,NOVRLP
               INCRI = INCR + NROWEL + I
               LOOP = NRWBLK + I
               DO 560 J = 1,NROWEL
                  INCRJ = INCR + J
                  B(INCRI) = B(INCRI) - ARRAY(J,LOOP,K)*B(INCRJ)
  560          CONTINUE
  570       CONTINUE

C           NOW, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U). THIS
C           CORRESPONDS TO THE LOOP 540 ABOVE.

            INCR = INCR + NROWEL
            DO 580 I = 1,NRWTOP
               INCRI = INCR + I
               IPVTI = PIVOT(INCRI)
               IF ( INCRI .NE. IPVTI ) THEN
                  SWAP = B(INCRI)
                  B(INCRI) = B(IPVTI)
                  B(IPVTI) = SWAP
               ENDIF
               LOOP = NROWEL + I
               BI = -B(INCRI)
               DO 575 J = I+1,NOVRLP
                  INCRJ = INCR+J
                  L = NRWBLK + J
                  B(INCRJ) = B(INCRJ) + BI*ARRAY(LOOP,L,K)
  575          CONTINUE
  580       CONTINUE
            INCR = INCR + NRWTOP
  590   CONTINUE

C       FINALLY, FINISH WITH NRWBOT SOLUTIONS:

        DO 600 J = 1,NRWBOT
           INCRJ = INCR + J
           DO 595 I = 1,J-1
              B(INCRJ) = B(INCRJ) - BOTBLK(I,J+NRWTOP)*B(INCR+I)
  595      CONTINUE
           B(INCRJ) = B(INCRJ)/BOTBLK(J,J+NRWTOP)
  600   CONTINUE


C       NOW, THE BACKWARD PASS:


C       FIRST, BACKWARD SOLUTION IN BOTBLK:

        INCRJ = INCR + NRWBOT
        DO 610 J = 1,NRWBOT-1
           INCRJ = INCRJ - 1
           DO 605 I = NRWBOT-J+1,NRWBOT
              INCRI = INCR + I
              B(INCRJ) = B(INCRJ) - BOTBLK(I,NOVRLP-J)*B(INCRI)
  605      CONTINUE

           IF ( INCRJ .NE. PIVOT(INCRJ) ) THEN
              SWAP = B(INCRJ)
              B(INCRJ) = B(PIVOT(INCRJ))
              B(PIVOT(INCRJ)) = SWAP
           ENDIF
  610   CONTINUE

C       NOW DO THE DEFERRED OPERATIONS IN BOTBLOK:

        DO 620 J = 1,NRWTOP
           INCRJ = INCR - J + 1
           DO 615 I = 1,NRWBOT
              INCRI = INCR + I
              B(INCRJ) = B(INCRJ) - BOTBLK(I,NRWTP1-J)*B(INCRI)
  615      CONTINUE
  620   CONTINUE


C       NOW, IN EACH BLOCK, K = NBLOKS,..,1:
        DO 800 K = NBLOKS,1,-1

C          FIRST, THE BACKSUBSTITUIONS:

           DO 630 J = 1,NRWTOP
              INCRJ = INCR - J + 1
              LOOP = NBKTOP - J + 1
              DO 625 I = 1,J-1
                 INCRI = INCR - I + 1
                 B(INCRJ) = B(INCRJ) - ARRAY(NRWBLK-I+1,LOOP,K)*B(INCRI)
  625         CONTINUE
              B(INCRJ) = B(INCRJ)/ARRAY(NRWBLK-J+1,LOOP,K)
  630      CONTINUE

C          THEN THE BACKWARD SOLUTION IN THE KTH BLOCK:

           DO 650 J = 1,NROWEL
              INCRJ = INCR - NRWTOP -J + 1
              DO 645 I = 1,J+NRWTOP-1
                 INCRI = INCRJ + I
                 B(INCRJ) = B(INCRJ) -
     *           ARRAY(NRWBLK-NRWTOP-J+1+I,NRWBLK-J+1,K)*B(INCRI)
  645         CONTINUE
              IF ( INCRJ .NE. PIVOT(INCRJ) ) THEN
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(PIVOT(INCRJ))
                 B(PIVOT(INCRJ)) = SWAP
              ENDIF
  650      CONTINUE

C          NOW, THE DEFERRED OPERATIONS ON B:

           INCR = INCR - NRWBLK
           DO 660 J = 1,NRWTOP
              INCRJ = INCR + J - NRWTOP
              DO 655 I = 1,NRWBLK
                 INCRI = INCR + I
                 B(INCRJ) = B(INCRJ) - ARRAY(I,J,K)*B(INCRI)
  655        CONTINUE
  660      CONTINUE
  800   CONTINUE

C       FINALLY, THE LAST SET OF BACK-SUBSTITUTIONS IN TOPBLK:

        DO 900 J = 1,NRWTOP
           INCRJ = NRWTOP -J + 1
           DO 850 I = INCRJ+1,NRWTOP
              B(INCRJ) = B(INCRJ) - TOPBLK(I,INCRJ)*B(I)
  850      CONTINUE
           B(INCRJ) = B(INCRJ)/TOPBLK(INCRJ,INCRJ)
  900   CONTINUE
C
C       RETURN FROM THE SOLUTION OF A-TRANSPOSE.X = B

        RETURN
        END
