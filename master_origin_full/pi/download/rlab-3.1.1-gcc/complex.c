/* complex.c */
/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle

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

/*
 * The functions with a capitolized name ( _Log ) take Complex as args.
 * They in turn call the lower level complex functions which take
 * doubles as args. At this point it appears that I will not use the
 * top-level functions. They will stay awhile, and if I don't use them
 * they will dissappear.
 */

#include "complex.h"
#include "util.h"

#include <math.h>


Complex complex_Pow_int(Complex c, int n)
{
  if (n < 0)
    return (1.0 / complex_Pow_int(c, -n));

  Complex result = 1;

  while (n > 1)
  {
    if (n % 2 == 1)
      result *= c;
    c *= c;
    n /= 2;
  }
  if (n > 0)
    result *= c;

  return result;
}

Complex complex_Mod (Complex c1, Complex c2)
{
  Complex q, r;

  q = c1 / c2;

  RE(q) = (int) RE(q);

  IM(q) = (int) IM(q);

  r = q * c2 - c1;

  return r;
}

