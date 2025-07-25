/*
 * mdrf3.h
 * Matrix-Dense-Real with and without FORTRAN77
 */

/*  This file is a part of rlabplus ("Our"-LaB)
    RLaB Copyright (C) 1995  Ian R. Searle,
    rlabplus (C) 2005-2018 M. Kostrun

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


int rotate3d(double *x, double *y, double *z, double kx, double ky, double kz, double theta);
int mdr_bwlabel (MDR *x, int conn, int * pixel_count, MDR **npix, MDR **bbox);


