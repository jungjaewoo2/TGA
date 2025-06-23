/* math_macros.h */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 2001-2018 Marijan Kostrun

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
#ifndef MATH_MACROS_L
#define MATH_MACROS_L

#define SWAP(a, b, type) { type c; c = b; b = a; a = c; }

#define SWAP2(a1, b1, type1, a2, b2, type2) { SWAP(a1, b1, type1); SWAP(a2, b2, type2);}

#define ABS(x)    ((x) >= 0  ? (x) : -(x))

#define MIN(x,y)  ((x) < (y) ? (x) :  (y))

#define MAX(x,y)  ((x) > (y) ? (x) :  (y))

#define dprintf( ... ) { printf(THIS_FILE ": " THIS_SOLVER ": " __VA_ARGS__ ); }

static inline double sign( double x )
{
  if (x > 0.0) return  1.0;
  if (x < 0.0) return -1.0;
  return x;
}

#endif
