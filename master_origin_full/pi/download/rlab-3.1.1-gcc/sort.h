/* sort.h: Matrix Dense Real Functions ... */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle
   Copyright (C) 2014 Marijan Kostrun

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

#ifndef RLAB_SORT_H
#define RLAB_SORT_H
extern int  rlab_sort_get_method(void);
extern int  rlab_sort_set_method(int method);
extern int  rlab_sort_get_nans(void);
extern int  rlab_sort_set_nans(int method);
extern void r_sort  (double *v, int left, int right, double *ind);
extern void i_sort  (int    *v, int left, int right, double *ind);
extern void csort   (char *v[], int left, int right, double *ind);
extern void i_qsort (int *v, int left, int right);
#endif  /* RLAB_SORT_H */
