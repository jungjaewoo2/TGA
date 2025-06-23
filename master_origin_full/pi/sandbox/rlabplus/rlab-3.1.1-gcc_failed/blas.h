/* bl.h */

/*
 * Defines, and declarations for a
 * C => fortran, RLaB => BLAS interface.
 */

#ifndef RLAB_BL_H
#define RLAB_BL_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
# define DDOT   ddot_
# define ZDOTU  zdotu_
# define DWAP   dswap_
# define DAXPY  daxpy_
# define DGEMM  dgemm_
# define ZGEMM  zgemm_
# define ZGEEV  zgeev_
# define ZGGEVX zggevx_
# define ZGGEV  zggev
# define ZCOPY  zcopy_
# define ZAXPY  zaxpy_
# define ZGEMV  zgemv_
# define ZSCAL  zscal_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
# define DDOT   _ddot
# define ZDOTU  _zdotu
# define DSWAP  _dswap
# define DAXPY  _daxpy
# define DGEMM  _dgemm
# define ZGEMM  _zgemm
# define ZGEMM  _zgemm
# define ZGEEV  _zgeev
# define ZGGEVX _zggevx
# define ZGGEV  _zggev
# define ZCOPY  _zcopy
# define ZAXPY  _zaxpy
# define ZGEMV  _zgemv
# define ZSCAL  _zscal
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
/* Do nothing, the existing code is OK */
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
# define DDOT   ddot
# define ZDOTU  zdotu
# define DSWAP  dswap
# define DAXPY  daxpy
# define DGEMM  dgemm
# define ZGEMM  zgemm
# define ZGEMM  zgemm
# define ZGEEV  zgeev
# define ZGGEVX zggevx
# define ZGGEV  zggev
# define ZCOPY  zcopy
# define ZAXPY  zaxpy
# define ZGEMV  zgemv
# define ZSCAL  zscal
#endif

extern int DSWAP ();
extern double  DDOT  (int*, double*, int *, double*, int*);
extern Complex ZDOTU (int*, Complex*, int *, Complex*, int*);
extern int DAXPY (int *, double *, double *, int *, double *, int *);
extern int DGEMM (char *, char *, int *, int *, int *, double *, double *,
                  int *, double *, int *, double *, double *, int *);
extern int ZGEMM (char *, char *, int *, int *, int *, Complex *, Complex *,
                  int *, Complex *, int *, Complex *, Complex *, int *);
extern int ZGEEV(char*,char*, int*, Complex*,int*, Complex*, Complex*,int*,
                  Complex*,int*, Complex*,int*, double*, int*);
extern int ZGGEVX(char* balanc, char* jobvl, char* jobvr, char* sense,
                   int* n, Complex* a, int* lda, Complex* b, int* ldb,
                   Complex* alpha, Complex* beta,
                   Complex* vl, int*, Complex* vr, int*,
                   int *ilo, int *ihi, double *lscale, double *rscale,
                   double *abnrm, double *bbnrm,
                   double *rconde, double *rcondv,
                   Complex *wrk, int *lwork, double *rwork, int *iwork,
                   int *bwork, int *info);
extern int ZGGEV(char* jobvl, char* jobvr,
                  int* n, Complex* a, int* lda, Complex* b, int* ldb,
                  Complex* alpha, Complex* beta,
                  Complex* vl, int*, Complex* vr, int*,
                  Complex *wrk, int *lwork, double *rwork, int *info);
extern int ZCOPY(int*, const Complex*,int*, Complex*,int*);
extern int ZAXPY(int*, Complex*, Complex*,int*, Complex*,int*);
extern int ZGEMV(char*, int*,int*, Complex*, Complex*,int*, Complex*,int*,
                  Complex*, Complex*,int*);
extern int ZSCAL(int*, Complex*, Complex*,int*);

#endif  /* RLAB_BL_H */
