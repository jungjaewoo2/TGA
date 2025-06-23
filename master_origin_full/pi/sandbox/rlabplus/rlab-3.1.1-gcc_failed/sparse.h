// sparse.h

//
// Defines, and declarations for a
// C => fortran, RLaB => Sparse interface.
//

#ifndef RLAB_SPARSE_H
#define RLAB_SPARSE_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK

#define DGSADD dgsadd_
#define DGSSUB dgssub_
#define DGSMUL dgsmul_
#define DGSTRN dgstrn_

#define XGSADD xgsadd_
#define XGSMUL xgsmul_
#define XGSTRN xgstrn_

#define ZGSADD zgsadd_
#define ZGSSUB zgssub_
#define ZGSMUL zgsmul_
#define ZGSTRN zgstrn_

#endif

#ifdef HAVE_FORTRAN_UND_FRONT

#define DGSADD _dgsadd
#define DGSSUB _dgssub
#define DGSMUL _dgsmul
#define DGSTRN _dgstrn

#define XGSADD _xgsadd
#define XGSMUL _xgsmul
#define XGSTRN _xgstrn

#define ZGSADD _zgsadd
#define ZGSSUB _zgssub
#define ZGSMUL _zgsmul
#define ZGSTRN _zgstrn

#endif

#ifdef HAVE_FORTRAN_UPPERCASE
/* Do nothing, the existing code is OK */
#endif

#ifdef HAVE_FORTRAN_LOWERCASE

#define DGSADD dgsadd
#define DGSSUB dgssub
#define DGSMUL dgsmul
#define DGSTRN dgstrn

#define XGSADD xgsadd
#define XGSMUL xgsmul
#define XGSTRN xgstrn

#define ZGSADD zgsadd
#define ZGSSUB zgssub
#define ZGSMUL zgsmul
#define ZGSTRN zgstrn

#endif

extern int DGSADD (int *, int *, double *, int *, int *, double *, int *, int *,
                     int *, int *, double *, double *);
extern int DGSSUB (int *, int *, double *, int *, int *, double *, int *, int *,
                     int *, int *, double *, double *);
extern int DGSMUL (int *, int *, double *, int *, int *, double *, int *, int *,
                     int *, int *, double *, double * );
extern int DGSTRN (int *, int *, double *, int *, int *, int *, int *, double *);

extern int XGSADD (int *, int *, int *, int *, int *, int *, int *, int *, int * );
extern int XGSMUL (int *, int *, int *, int *, int *, int *, int *, int *, int *, int * );
extern int XGSTRN (int *, int *, int *, int *, int *, int *);

extern int ZGSADD (int *, int *, Complex *, int *, int *, Complex *, int *, int *,
                   int *, int *, Complex *, Complex *);
extern int ZGSSUB (int *, int *, Complex *, int *, int *, Complex *, int *, int *,
                   int *, int *, Complex *, Complex *);
extern int ZGSMUL (int *, int *, Complex *, int *, int *, Complex *, int *, int *,
                   int *, int *, Complex *, Complex *);
extern int ZGSTRN (int *, int *, Complex *, int *, int *, int *, int *, Complex *);

#endif  /* RLAB_SPARSE_H */
