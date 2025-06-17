/* complex.h */
/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle
   Copyright (C) 2014  Marijan Kostrun

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

// In Memoriam:
//  struct _rlab_complex
//  {
//    double r;      /* Real part */
//    double i;      /* Imaginary part */
//  };
//  typedef struct _rlab_complex Complex;

#ifndef COMPLEX_H
#define COMPLEX_H

#include <complex.h>
typedef double _Complex Complex;

#define  RE(x)        (((double *) &(x))[0])
#define  IM(x)        (((double *) &(x))[1])
#define  SQR(x)       ((double) RE(x) * RE(x) + IM(x) * IM(x))
#define  EQCMPL(x,y)  ((RE(x) == RE(y)) && (IM(x) == IM(y)))
#define  CMPLZERO(x)  ((RE(x) == 0) && (IM(x) == 0))
#define  CMPLREAL(x)  (IM(x) == 0)
#define  CMPLIMAG(x)  (RE(x) == 0)

extern Complex complex_Mod     (Complex, Complex);
extern Complex complex_Pow_int (Complex, int);

#endif /* COMPLEX_H */
