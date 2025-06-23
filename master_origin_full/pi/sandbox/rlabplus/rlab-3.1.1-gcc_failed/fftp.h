/* fftp.h */

#ifndef RLAB_FFTP_H
#define RLAB_FFTP_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
# define CFFTI dcffti_
# define CFFTF dcfftf_
# define CFFTB dcfftb_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
# define CFFTI _dcffti
# define CFFTF _dcfftf
# define CFFTB _dcfftb
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
# define CFFTI DCFFTI
# define CFFTF DCFFTF
# define CFFTB DCFFTB
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
# define CFFTI dcffti
# define CFFTF dcfftf
# define CFFTB dcfftb
#endif

extern int CFFTI (int *, double *);
extern int CFFTF (int *, double *, double *);
extern int CFFTB (int *, double *, double *);

#endif  /* RLAB_FFTP_H */
