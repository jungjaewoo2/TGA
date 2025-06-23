/* mdrf2.h: Matrix Dense Real Functions ... */

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

#ifndef RLAB_MDRF3_H
#define RLAB_MDRF3_H

#include "rlab.h"
#include "ent.h"
#include "btree.h"
#include "mdr.h"
#include "mdc.h"
#include "msr.h"
#include "msc.h"
#include "complex.h"


#include <stdio.h>
#include <math.h>

//
// eigs(a)
//
extern Btree * mdr_ar_EigsS (MDR * a, int nev, char *which, double sigma);
extern Btree * mdr_ar_EigsG (MDR * a, MDR *b, int nev, char *which, double sigma);
extern Btree * mdc_ar_EigsS (MDC * a, int nev, char *which, Complex sigma);
extern Btree * mdc_ar_EigsG (MDC * a, MDC *b, int nev, char *which, Complex sigma);
extern Btree * msr_ar_EigsS (MSR * a, int nev, char *which, double sigma);
extern Btree * msr_ar_EigsG (MSR * a, MSR * b, int nev, char *which, double sigma);
extern Btree * msc_ar_EigsS (MSC * a, int nev, char *which, Complex sigma);
extern Btree * msc_ar_EigsG (MSC * a, MSC * b, int nev, char *which, Complex sigma);

#endif /* RLAB_MDRF2_H */
