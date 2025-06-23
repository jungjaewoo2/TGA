*
*     lp.f: RLaB => LAPACK Interface
*
***********************************************************************

      SUBROUTINE RGEMV( TRANS, M, N, ALFA, A, LDA, X, INCX, BETA, Y,
     $                  INCY)
      INTEGER            TRANS
      INTEGER            M, N, LDA, INCX, INCY
      DOUBLE PRECISION   A( LDA, * ), X(*), Y(*)
      DOUBLE PRECISION   ALFA, BETA

      CALL DGEMV(TRANS, M, N, ALFA, A, LDA, X, INCX, BETA, Y, INCY)
      RETURN
      END


***********************************************************************

      SUBROUTINE RGETRF( M, N, A, LDA, IPIV, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )

      CALL DGETRF( M, N, A, LDA, IPIV, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGETRF( M, N, A, LDA, IPIV, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )

      CALL ZGETRF( M, N, A, LDA, IPIV, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      INTEGER            INFO, LDA, LWORK, N
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( LWORK )

      CALL DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      INTEGER            INFO, LDA, LWORK, N
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( LWORK )

      CALL ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      INTEGER            TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )

      CALL DGETRS( CHAR(TRANS), N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      INTEGER            TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )

      CALL ZGETRS( CHAR(TRANS), N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,
     $                   INFO )
      INTEGER            NORM
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )

      CALL DGECON( CHAR(NORM), N, A, LDA, ANORM, RCOND, WORK, IWORK,
     $                   INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK,
     $                   INFO )
      INTEGER            NORM
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), WORK( * )

      CALL ZGECON( CHAR(NORM), N, A, LDA, ANORM, RCOND, WORK, RWORK,
     $             INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
      INTEGER            JOB
      INTEGER            IHI, ILO, INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * ), SCALE( * )

      CALL DGEBAL( CHAR(JOB), N, A, LDA, ILO, IHI, SCALE, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
      INTEGER            JOB
      INTEGER            IHI, ILO, INFO, LDA, N
      DOUBLE PRECISION   SCALE( * )
      COMPLEX*16         A( LDA, * )

      CALL ZGEBAL( CHAR(JOB), N, A, LDA, ILO, IHI, SCALE, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            INFO, LDA, LWORK, M, N
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )

      CALL DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            INFO, LDA, LWORK, M, N
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( LWORK )

      CALL ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            INFO, K, LDA, LWORK, M, N
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )

      CALL DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            INFO, K, LDA, LWORK, M, N
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( LWORK )

      CALL ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      DOUBLE PRECISION FUNCTION RLANGE( NORM, M, N, A, LDA, WORK )
      INTEGER            NORM
      INTEGER            LDA, M, N
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
      DOUBLE PRECISION DLANGE

      RLANGE =  DLANGE( CHAR(NORM), M, N, A, LDA, WORK )
      RETURN
      END

***********************************************************************

      DOUBLE PRECISION FUNCTION XLANGE( NORM, M, N, A, LDA, WORK )
      INTEGER            NORM
      INTEGER            LDA, M, N
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( LDA, * )
      DOUBLE PRECISION ZLANGE

      XLANGE = ZLANGE( CHAR(NORM), M, N, A, LDA, WORK )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
     $                   WORK, LWORK, INFO )
      INTEGER            JOBU, JOBVT
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
      DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
     $                   VT( LDVT, * ), WORK( * )

      CALL DGESVD( CHAR(JOBU), CHAR(JOBVT), M, N, A, LDA, S, U, LDU,
     $             VT, LDVT, WORK, LWORK, INFO )

      RETURN
      END

***********************************************************************

      SUBROUTINE XGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
     $                   WORK, LWORK, RWORK, INFO )
      INTEGER            JOBU, JOBVT
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
      DOUBLE PRECISION   RWORK( * ), S( * )
      COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ),
     $                   WORK( * )

      CALL ZGESVD( CHAR(JOBU), CHAR(JOBVT), M, N, A, LDA, S, U, LDU,
     $             VT, LDVT, WORK, LWORK, RWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RPOTRF( UPLO, N, A, LDA, INFO )
      INTEGER            UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * )

      CALL DPOTRF( CHAR(UPLO), N, A, LDA, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XPOTRF( UPLO, N, A, LDA, INFO )
      INTEGER            UPLO
      INTEGER            INFO, LDA, N
      COMPLEX*16         A( LDA, * )

      CALL ZPOTRF( CHAR(UPLO), N, A, LDA, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
      INTEGER            JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WI( * ), WORK( * ), WR( * )

      CALL DGEEV( CHAR(JOBVL), CHAR(JOBVR), N, A, LDA, WR, WI, VL,
     $            LDVL, VR, LDVR, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
     $                  WORK, LWORK, RWORK, INFO )
      INTEGER            JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   W( * ), WORK( * )

      CALL ZGEEV( CHAR(JOBVL), CHAR(JOBVR), N, A, LDA, W, VL, LDVL,
     $            VR, LDVR, WORK, LWORK, RWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )

      CALL DGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( LWORK )

      CALL ZGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )

      CALL DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( LWORK )

      CALL ZUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
     $                  LWORK, INFO )
      INTEGER            JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * )

      CALL DSYGV( ITYPE, CHAR(JOBZ), CHAR(UPLO), N, A, LDA, B, LDB,
     $            W, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XHEGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
     $                  LWORK, RWORK, INFO )
      INTEGER            JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )

      CALL ZHEGV( ITYPE, CHAR(JOBZ), CHAR(UPLO), N, A, LDA, B, LDB,
     $            W, WORK, LWORK, RWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            JPVT( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )

      CALL DGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            JPVT( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )

      CALL ZGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      INTEGER            JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )

      CALL DSYEV( CHAR(JOBZ), CHAR(UPLO), N, A, LDA, W, WORK,
     $            LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
     $                  INFO )
      INTEGER            JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * )

      CALL ZHEEV( CHAR(JOBZ), CHAR(UPLO), N, A, LDA, W, WORK, LWORK,
     $            RWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
     $                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
      INTEGER            JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )

      CALL DGGEV( CHAR(JOBVL), CHAR(JOBVR), N, A, LDA, B, LDB, ALPHAR,
     $            ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
     $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
      INTEGER            JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( * )

      CALL ZGGEV( CHAR(JOBVL), CHAR(JOBVR), N, A, LDA, B, LDB, ALPHA,
     $            BETA, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
     $                   WORK, LWORK, INFO )
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
      DOUBLE PRECISION   RCOND
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )

      CALL DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
     $                   WORK, LWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
     $                   WORK, LWORK, RWORK, INFO )
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
      DOUBLE PRECISION   RCOND
      DOUBLE PRECISION   RWORK( * ), S( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )

      CALL ZGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
     $                   WORK, LWORK, RWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI,
     $                  VS, LDVS, WORK, LWORK, BWORK, INFO )
      INTEGER            JOBVS, SORT
      INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
      LOGICAL            BWORK( * )
      DOUBLE PRECISION   A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),
     $                   WR( * )
      LOGICAL            SELECT
      EXTERNAL           SELECT

      CALL DGEES( CHAR(JOBVS), CHAR(SORT), SELECT, N, A, LDA, SDIM,
     $            WR, WI, VS, LDVS, WORK, LWORK, BWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,
     $                  LDVS, WORK, LWORK, RWORK, BWORK, INFO )
      INTEGER            JOBVS, SORT
      INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
      LOGICAL            BWORK( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
      LOGICAL            SELECT
      EXTERNAL           SELECT

      CALL ZGEES( CHAR(JOBVS), CHAR(SORT), SELECT, N, A, LDA, SDIM,
     $            W, VS, LDVS, WORK, LWORK, RWORK, BWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,
     $                   INFO )
      INTEGER            COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )

      CALL DTREXC( CHAR(COMPQ), N, T, LDT, Q, LDQ, IFST, ILST, WORK,
     $                   INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
      INTEGER            COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
      COMPLEX*16         Q( LDQ, * ), T( LDT, * )

      CALL ZTREXC( CHAR(COMPQ), N, T, LDT, Q, LDQ, IFST, ILST, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
     $                   LDC, SCALE, INFO )
      INTEGER            TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      DOUBLE PRECISION   SCALE
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )

      CALL DTRSYL( CHAR(TRANA), CHAR(TRANB), ISGN, M, N, A, LDA,
     $             B, LDB, C, LDC, SCALE, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
     $                   LDC, SCALE, INFO )
      INTEGER            TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      DOUBLE PRECISION   SCALE
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )

      CALL ZTRSYL( CHAR(TRANA), CHAR(TRANB), ISGN, M, N, A, LDA,
     $             B, LDB, C, LDC, SCALE, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      INTEGER            UPLO
      INTEGER            INFO, LDA, LWORK, N
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( LWORK )

      CALL DSYTRF( CHAR(UPLO), N, A, LDA, IPIV, WORK, LWORK, INFO )
      RETURN
      END

************************************************************************

      SUBROUTINE XHETRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
      INTEGER            UPLO
      INTEGER            INFO, LDA, LWORK, N
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( LWORK )

      CALL ZHETRF( CHAR(UPLO), N, A, LDA, IPIV, WORK, LWORK, INFO )
      RETURN
      END

************************************************************************

      SUBROUTINE XHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
     $                  LWORK, INFO )
*     .. Scalar Arguments ..
      INTEGER            UPLO
      INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( LWORK )
*     ..


      CALL ZHESV( CHAR(UPLO), N, NRHS, A, LDA, IPIV, B, LDB,
     $            WORK, LWORK, INFO )

      RETURN
      END

***********************************************************************

      SUBROUTINE RSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      INTEGER            UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )

      CALL DSYTRS( CHAR(UPLO), N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      INTEGER            UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )

      CALL ZHETRS( CHAR(UPLO), N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE RSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,
     $                   IWORK, INFO )
      INTEGER            UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )

      CALL DSYCON( CHAR(UPLO), N, A, LDA, IPIV, ANORM, RCOND, WORK,
     $             IWORK, INFO )
      RETURN
      END

***********************************************************************

      SUBROUTINE XHECON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,
     $                   INFO )
      INTEGER            UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( * )

      CALL ZHECON( CHAR(UPLO), N, A, LDA, IPIV, ANORM, RCOND, WORK,
     $             INFO )
      RETURN
      END

***********************************************************************

      DOUBLE PRECISION FUNCTION RLAMCH( CMACH )
      INTEGER CMACH

      RLAMCH = DLAMCH( CHAR(CMACH) )
      RETURN
      END

***********************************************************************
