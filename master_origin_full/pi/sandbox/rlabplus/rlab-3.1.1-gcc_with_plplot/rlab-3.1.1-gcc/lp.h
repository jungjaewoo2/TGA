
/* lp.h */

/*
 * Defines, and declarations for a
 * C => Fortran, RLaB => LAPACK interface.
 */
// Sources:
//    http://unifycelerity.googlecode.com/svn/trunk/Include/mcrinclude/lapack.h
#ifndef RLAB_LP_H
#define RLAB_LP_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
#define RGEMV  rgemv_
#define RGETRF rgetrf_
#define XGETRF xgetrf_
#define RGETRI rgetri_
#define XGETRI xgetri_
#define RGETRS rgetrs_
#define XGETRS xgetrs_
#define RGECON rgecon_
#define XGECON xgecon_
#define RGEBAL rgebal_
#define XGEBAL xgebal_
#define RGEQRF rgeqrf_
#define XGEQRF xgeqrf_
#define RORGQR rorgqr_
#define XUNGQR xungqr_
#define RGESVD rgesvd_
#define XGESVD xgesvd_
#define RPOTRF rpotrf_
#define XPOTRF xpotrf_
#define RGEEV rgeev_
#define XGEEV xgeev_
#define RGEHRD rgehrd_
#define XGEHRD xgehrd_
#define RORGHR rorghr_
#define XUNGHR xunghr_
#define RSYGV rsygv_
#define XHEGV xhegv_
#define RGEQPF rgeqpf_
#define XGEQPF xgeqpf_
#define RSYEV rsyev_
#define XHEEV xheev_
#define RGEGV rgegv_
#define XGEGV xgegv_
#define RGELSS rgelss_
#define XGELSS xgelss_
#define RGETRF rgetrf_
#define XGETRF xgetrf_
#define RGEES rgees_
#define XGEES xgees_
#define RTREXC rtrexc_
#define XTREXC xtrexc_
#define RTRSYL rtrsyl_
#define XTRSYL xtrsyl_
#define RSYTRF rsytrf_
#define XSYTRF xsytrf_
#define RSYCON rsycon_
#define XSYCON xsycon_
#define RSYTRS rsytrs_
#define XSYTRS xsytrs_
#define XHETRF xhetrf_
#define XHETRS xhetrs_
#define XHECON xhecon_
#define XHESV  xhesv_
#define RLANGE rlange_
#define XLANGE xlange_
#define RLAMCH rlamch_

#endif

#ifdef HAVE_FORTRAN_UND_FRONT
#define RGEMV  _rgemv
#define RGETRF _rgetrf
#define XGETRF _xgetrf
#define RGETRI _rgetri
#define XGETRI _xgetri
#define RGETRS _rgetrs
#define XGETRS _xgetrs
#define RGECON _rgecon
#define XGECON _xgecon
#define RGEBAL _rgebal
#define XGEBAL _xgebal
#define RGEQRF _rgeqrf
#define XGEQRF _xgeqrf
#define RORGQR _rorgqr
#define XUNGQR _xungqr
#define RGESVD _rgesvd
#define XGESVD _xgesvd
#define RPOTRF _rpotrf
#define XPOTRF _xpotrf
#define RGEEV _rgeev
#define XGEEV _xgeev
#define RGEHRD _rgehrd
#define XGEHRD _xgehrd
#define RORGHR _rorghr
#define XUNGHR _xunghr
#define RSYGV _rsygv
#define XHEGV _xhegv
#define RGEQPF _rgeqpf
#define XGEQPF _xgeqpf
#define RSYEV _rsyev
#define XHEEV _xheev
#define RGEGV _rgegv
#define XGEGV _xgegv
#define RGELSS _rgelss
#define XGELSS _xgelss
#define RGETRF _rgetrf
#define XGETRF _xgetrf
#define RGEES _rgees
#define XGEES _xgees
#define RTREXC _rtrexc
#define XTREXC _xtrexc
#define RTRSYL _rtrsyl
#define XTRSYL _xtrsyl
#define RSYTRF _rsytrf
#define XSYTRF _xsytrf
#define RSYCON _rsycon
#define XSYCON _xsycon
#define RSYTRS _rsytrs
#define XSYTRS _xsytrs
#define XHETRF _xhetrf
#define XHETRS _xhetrs
#define XHECON _xhecon
#define XHESV  _xhesv
#define RLANGE _rlange
#define XLANGE _xlange
#define RLAMCH _rlamch
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
/* Do nothing, the existing code is OK */
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
#define RGEMV  rgemv
#define RGETRF rgetrf
#define XGETRF xgetrf
#define RGETRI rgetri
#define XGETRI xgetri
#define RGETRS rgetrs
#define XGETRS xgetrs
#define RGECON rgecon
#define XGECON xgecon
#define RGEBAL rgebal
#define XGEBAL xgebal
#define RGEQRF rgeqrf
#define XGEQRF xgeqrf
#define RORGQR rorgqr
#define XUNGQR xungqr
#define RGESVD rgesvd
#define XGESVD xgesvd
#define RPOTRF rpotrf
#define XPOTRF xpotrf
#define RGEEV rgeev
#define XGEEV xgeev
#define RGEHRD rgehrd
#define XGEHRD xgehrd
#define RORGHR rorghr
#define XUNGHR xunghr
#define RSYGV rsygv
#define XHEGV xhegv
#define RGEQPF rgeqpf
#define XGEQPF xgeqpf
#define RSYEV rsyev
#define XHEEV xheev
#define RGEGV rgegv
#define XGEGV xgegv
#define RGELSS rgelss
#define XGELSS xgelss
#define RGETRF rgetrf
#define XGETRF xgetrf
#define RGEES rgees
#define XGEES xgees
#define RTREXC rtrexc
#define XTREXC xtrexc
#define RTRSYL rtrsyl
#define XTRSYL xtrsyl
#define RSYTRF rsytrf
#define XSYTRF xsytrf
#define RSYCON rsycon
#define XSYCON xsycon
#define RSYTRS rsytrs
#define XSYTRS xsytrs
#define XHETRF xhetrf
#define XHETRS xhetrs
#define XHECON xhecon
#define XHESV  xhesv
#define RLANGE rlange
#define XLANGE xlange
#define RLAMCH rlamch

#endif

extern int RGEMV  (int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);
extern int RGETRF (int *, int *, double *, int *, int *, int *);
extern int XGETRF (int *m, int *n, Complex *a, int *lda, int *ipiv, int *info);
extern int RGETRI (int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
extern int XGETRI (int *n, Complex *a, int *lda, int *ipiv, Complex *work, int *lwork, int *info);
extern int RGETRS (int *, int *, int *, double *, int *, int *, double *, int *, int *);
extern int XGETRS (int    *trans,
                   int    *n,
                   int    *nrhs,
                   Complex *a,
                   int    *lda,
                   int    *ipiv,
                   Complex *b,
                   int    *ldb,
                   int    *info);
extern int RGECON (int *, int *, double *, int *, double *, double *, double *, int *, int *);
extern int XGECON (int    *norm,
                   int    *n,
                   Complex *a,
                   int    *lda,
                   double *anorm,
                   double *rcond,
                   Complex *work,
                   double *rwork,
                   int    *info);
extern int RGEBAL (int    *job,
                   int    *n,
                   double *a,
                   int    *lda,
                   int    *ilo,
                   int    *ihi,
                   double *scale,
                   int    *info);
extern int XGEBAL (int    *job,
                   int    *n,
                   Complex *a,
                   int    *lda,
                   int    *ilo,
                   int    *ihi,
                   double *scale,
                   int    *info);
extern int RGEQRF (int    *m,
                   int    *n,
                   double *a,
                   int    *lda,
                   double *tau,
                   double *work,
                   int    *lwork,
                   int    *info);
extern int XGEQRF (int    *m,
                   int    *n,
                   Complex *a,
                   int    *lda,
                   Complex *tau,
                   Complex *work,
                   int    *lwork,
                   int    *info);
extern int RORGQR (int    *m,
                   int    *n,
                   int    *k,
                   double *a,
                   int    *lda,
                   double *tau,
                   double *work,
                   int    *lwork,
                   int    *info);
extern int XUNGQR (int    *m,
                   int    *n,
                   int    *k,
                   Complex *a,
                   int    *lda,
                   Complex *tau,
                   Complex *work,
                   int    *lwork,
                   int    *info);
extern int RGESVD (int *, int *, int *, int *, double *, int *, double *,
                   double *, int *, double *, int *, double *, int *, int *);
extern int XGESVD (int *jobu, int *jobvt, int *m, int *n, Complex *a,
                   int *lda, double *s, Complex *u, int *ldu, Complex *vt,
                   int *ldvt, Complex *work, int *lwork, double *rwork, int *info);
extern int RPOTRF (int  *uplo, int *n, double *a, int *lda, int *info);
extern int XPOTRF (int  *uplo, int *n, Complex *a, int *lda, int *info);
extern int RGEEV (int *, int *, int *, double *, int *, double *, double *,
                  double *, int *, double *, int *, double *,
                  int *, int *);
extern int XGEEV (int *jobvl, int *jobvr, int *n, Complex *a, int *lda, Complex *w, Complex *vl,
                  int *ldvl, Complex *vr, int *ldvr, Complex *work, int *lwork, double *rwork, int *info);
extern int RGEHRD (int    *n,
                   int    *ilo,
                   int    *ihi,
                   double *a,
                   int    *lda,
                   double *tau,
                   double *work,
                   int    *lwork,
                   int    *info);
extern int XGEHRD (int    *n,
                   int    *ilo,
                   int    *ihi,
                   Complex *a,
                   int    *lda,
                   Complex *tau,
                   Complex *work,
                   int    *lwork,
                   int    *info);
extern int RORGHR (int    *n,
                   int    *ilo,
                   int    *ihi,
                   double *a,
                   int    *lda,
                   double *tau,
                   double *work,
                   int    *lwork,
                   int    *info);
extern int XUNGHR (int    *n,
                   int    *ilo,
                   int    *ihi,
                   Complex *a,
                   int    *lda,
                   Complex *tau,
                   Complex *work,
                   int    *lwork,
                   int    *info);
extern int RSYGV (int *, int *, int *, int *, double *, int *, double *, int *,
                  double *, double *, int *, int *);
extern int XHEGV (int    *itype,
                  int    *jobz,
                  int    *uplo,
                  int    *n,
                  Complex *a,
                  int    *lda,
                  Complex *b,
                  int    *ldb,
                  double *w,
                  Complex *work,
                  int    *lwork,
                  double *rwork,
                  int    *info);
extern int RGEQPF (int    *m,
                   int    *n,
                   double *a,
                   int    *lda,
                   int    *jpvt,
                   double *tau,
                   double *work,
                   int    *info);
extern int XGEQPF (int    *m,
                   int    *n,
                   Complex *a,
                   int    *lda,
                   int    *jpvt,
                   Complex *tau,
                   Complex *work,
                   double *rwork,
                   int    *info);
extern int RSYEV (int *, int *, int *, double *, int *, double *, double *, int *, int *);
extern int XHEEV (int *jobz, int *uplo, int *n, Complex *a, int *lda, double *w,
                  Complex *work, int *lwork, double *rwork, int *info);
extern int RGEGV (int *, int *, int *, double *, int *, double *, int *,
                  double *, double *, double *,
                  double *, int *, double *, int *, double *, int *, int *);
extern int XGEGV (int   *jobvl,
                  int   *jobvr,
                  int    *n,
                  Complex *a,
                  int    *lda,
                  Complex *b,
                  int    *ldb,
                  Complex *alpha,
                  Complex *beta,
                  Complex *vl,
                  int    *ldvl,
                  Complex *vr,
                  int    *ldvr,
                  Complex *work,
                  int    *lwork,
                  double *rwork,
                  int    *info);
extern int RGELSS (int *, int *, int *, double *, int *, double *, int *, double *,
                   double *, int *, double *, int *, int *);
extern int XGELSS (int    *m,
                   int    *n,
                   int    *nrhs,
                   Complex *a,
                   int    *lda,
                   Complex *b,
                   int    *ldb,
                   double *s,
                   double *rcond,
                   int    *rank,
                   Complex *work,
                   int    *lwork,
                   double *rwork,
                   int    *info);
extern int RGETRF (int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern int XGETRF (int *m, int *n, Complex *a, int *lda, int *ipiv, int *info);
extern int RGEES (int   *jobvs,
                  int   *sort,
                  void  (*)(),
                  int    *n,
                  double *a,
                  int    *lda,
                  int    *sdim,
                  double *wr,
                  double *wi,
                  double *vs,
                  int    *ldvs,
                  double *work,
                  int    *lwork,
                  int    *bwork,
                  int    *info);
extern int XGEES (int   *jobvs,
                  int   *sort,
                  void    (*)(),
                  int    *n,
                  Complex *a,
                  int    *lda,
                  int    *sdim,
                  Complex *w,
                  Complex *vs,
                  int    *ldvs,
                  Complex *work,
                  int    *lwork,
                  double *rwork,
                  int    *bwork,
                  int    *info);
extern int RTREXC (int    *compq,
                   int    *n,
                   double *t,
                   int    *ldt,
                   double *q,
                   int    *ldq,
                   int    *ifst,
                   int    *ilst,
                   double *work,
                   int    *info);
extern int XTREXC (int   *compq,
                   int    *n,
                   Complex *t,
                   int    *ldt,
                   Complex *q,
                   int    *ldq,
                   int    *ifst,
                   int    *ilst,
                   int    *info);
extern int RTRSYL (int *, int *, int *, int *, int *, double *, int *,
                   double *, int *, double *, int *, double *, int *);
extern int XTRSYL (int   *trana,
                   int   *tranb,
                   int    *isgn,
                   int    *m,
                   int    *n,
                   Complex *a,
                   int    *lda,
                   Complex *b,
                   int    *ldb,
                   Complex *c,
                   int    *ldc,
                   double *scale,
                   int    *info);
extern int RSYTRF (int *, int *, double *, int *, int *, double *, int *, int *);
extern int XSYTRF (int   *uplo,
                   int    *n,
                   Complex *a,
                   int    *lda,
                   int    *ipiv,
                   Complex *work,
                   int    *lwork,
                   int    *info);
extern int RSYCON (int *, int *, double *, int *, int *, double *, double *, double *, int *, int *);
extern int XSYCON (int   *uplo,
                   int    *n,
                   Complex *a,
                   int    *lda,
                   int    *ipiv,
                   double *anorm,
                   double *rcond,
                   Complex *work,
                   int    *info);
extern int RSYTRS (int *, int *, int *, double *, int *, int *, double *, int *, int *);
extern int XSYTRS (int   *uplo,
                   int    *n,
                   int    *nrhs,
                   Complex *a,
                   int    *lda,
                   int    *ipiv,
                   Complex *b,
                   int    *ldb,
                   int    *info);
extern int XHETRF (int    *uplo,
                   int    *n,
                   Complex *a,
                   int    *lda,
                   int    *ipiv,
                   Complex *work,
                   int    *lwork,
                   int    *info);
extern int XHETRS (int   *uplo,
                   int    *n,
                   int    *nrhs,
                   Complex *a,
                   int    *lda,
                   int    *ipiv,
                   Complex *b,
                   int    *ldb,
                   int    *info);
extern int XHECON (int   *uplo,
                   int    *n,
                   Complex *a,
                   int    *lda,
                   int    *ipiv,
                   double *anorm,
                   double *rcond,
                   Complex *work,
                   int    *info);
extern int XHESV  (int    *uplo,
                   int    *n,
                   int    *nrhs,
                   Complex *a,
                   int    *lda,
                   int    *ipiv,
                   Complex *b,
                   int    *ldb,
                   Complex *work,
                   int    *lwork,
                   int    *info);
extern double RLANGE (int *, int *, int *, double *, int *, double *);
extern double XLANGE (int   *norm,
                      int    *m,
                      int    *n,
                      Complex *a,
                      int    *lda,
                      double *work);
extern double RLAMCH ();

#endif /* RLAB_LP_H */
