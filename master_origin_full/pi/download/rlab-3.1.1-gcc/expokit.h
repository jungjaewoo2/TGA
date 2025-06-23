/* expokit.h: All Functions ... */

/*  This file is a part of RLaB ("Our"-LaB) + rlabplus
   Copyright (C) 2007  Marijan Kostrun

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   See the file ./COPYING
   ********************************************************************** */

#ifndef RLAB_MEXP_H
#define RLAB_MEXP_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define DMEXPV dmexpv_
  #define DGPADM dgpadm_
  #define DSPADM dspadm_
  #define ZGPADM zgpadm_
  #define ZHPADM zhpadm_
  #define DGCHBV dgchbv_
  #define DSCHBV dschbv_
  #define ZGCHBV zgchbv_
  #define DNCHBV dnchbv_
  #define ZNCHBV znchbv_
  #define DGEXPV dgexpv_
  #define DSEXPV dsexpv_
  #define ZGEXPV zgexpv_
  #define ZHEXPV zhexpv_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define DMEXPV _dmexpv
  #define DGPADM _dgpadm
  #define DSPADM _dspadm
  #define ZGPADM _zgpadm
  #define ZHPADM _zhpadm
  #define DGCHBV _dgchbv
  #define DSCHBV _dschbv
  #define ZGCHBV _zgchbv
  #define DNCHBV _dnchbv
  #define ZNCHBV _znchbv
  #define DGEXPV _dgexpv
  #define DSEXPV _dsexpv
  #define ZGEXPV _zgexpv
  #define ZHEXPV _zhexpv
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define DMEXPV dmexpv
  #define DGPADM dgpadm
  #define DSPADM dspadm
  #define ZGPADM zgpadm
  #define ZHPADM zhpadm
  #define DGCHBV dgchbv
  #define DSCHBV dschbv
  #define ZGCHBV zgchbv
  #define DNCHBV dnchbv
  #define ZNCHBV znchbv
  #define DGEXPV dgexpv
  #define DSEXPV dsexpv
  #define ZGEXPV zgexpv
  #define ZHEXPV zhexpv
#endif

//
// FORTRAN77 solvers
//
// extern int DMEXPV ();
extern int DGPADM (int *, int *, double *, double *, int *, double *, int *, int *,
                   int *, int *, int *);
extern int DSPADM (int *, int *, double *, double *, int *, double *, int *, int *,
                   int *, int *, int *);
extern int ZGPADM (int *, int *, double *, Complex *, int *, Complex *, int *, int *,
                   int *, int *, int *);
extern int ZHPADM (int *, int *, double *, Complex *, int *, Complex *, int *, int *,
                   int *, int *, int *);
extern int DGCHBV (int *, double *, double *, int *, double *, double *, int *, int *);
extern int DSCHBV (int *, double *, double *, int *, double *, double *, int *, int *);
extern int ZGCHBV (int *, double *, Complex *, int *, Complex *, Complex *, int *, int *);
// extern int DNCHBV ();
// extern int ZNCHBV ();
extern int DGEXPV (int *, int *, double *, double *, double *, double *, double *,
                   double *, int *, int *, int *,
                   int (*)(double *, double *),
                   int *, int *);
extern int DSEXPV (int *, int *, double *, double *, double *, double *, double *,
                   double *, int *, int *, int *,
                   int (*)(double *, double *),
                   int *, int *);
extern int ZHEXPV (int *, int *, double *, Complex *, Complex *, double *, double *,
                   Complex *, int *, int *, int *,
                   int (*)(Complex *, Complex *),
                   int *, int *);
extern int ZGEXPV (int *, int *, double *, Complex *, Complex *, double *, double *,
                   Complex *, int *, int *, int *,
                   int (*)(Complex *, Complex *),
                   int *, int *);
#endif // RLAB_MEXP_H
