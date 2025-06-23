// msrf3.c Matrix Sparse Real Functions

//  This file is a part of RLaB ("Our"-LaB)
//
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
//   See the file ./COPYING
//   ********************************************************************** */

// sparskit2 support by Marijan Kostrun, VI/2005

#include "rlab.h"
#include "mem.h"
#include "msr.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdc.h"
#include "mds.h"
#include "btree.h"
#include "util.h"
#include "symbol.h"
#include "sort.h"
#include "mathl.h"

#include "class.h"

#include "fi.h"
#include "sparskit.h"


#include <stdio.h>
#include <math.h>

//
// write a matrix in sparse Harwell-Boing format
//
void *
msr_sparskit_prtmt (MSR * m, char *ime)
{
  int i, unit = 66, job = 2, ifmt = 11;
  char title[72] = "RLaB2 Rel.2 SPARKIT library";
  char key[8] = "rlab", mtype[3] = "SP", guesol[3] = "NN";
  int nrow = m->nr, ncol = m->nc;
  i = strlen (ime);
  OPENFILE (ime, &i, &unit);
  PRTMT (&nrow, &ncol, m->d, m->ja, m->ia, NULL, guesol, title, key, mtype,
	 &ifmt, &job, &unit);
  CLOSFILE (&unit);
  return 0;
}

//
// write a sparse matrix in coordinate format
//
void *
msr_sparskit_smms (MSR * m, char *ime)
{
  int unit = 66, i;
  int n = m->nr, first = 1, last = m->nr, mode = 0;

  i = strlen (ime);
  OPENFILE (ime, &i, &unit);
  SMMS (&n, &first, &last, &mode, m->d, m->ja, m->ia, &unit);
  CLOSFILE (&unit);
  return 0;
}

//
// write a sparse matrix in compressed sparse row format
//
void *
msr_sparskit_dump (MSR * m, char *ime)
{
  int unit = 66, i;
  int first = 1, last = m->nr, values = 1;

  i = strlen (ime);
  OPENFILE (ime, &i, &unit);
  DUMP (&first, &last, &values, m->d, m->ja, m->ia, &unit);
  CLOSFILE (&unit);
  return 0;
}

//
// plot a sparse matrix as a postscript file
//
void *
msr_sparskit_pspltm (MSR * m, char *ime)
{
  int unit = 66, i;
  int nrow = m->nr, ncol = m->nc, mode = 0, ptitle = 0, nlines = 0;
  char munt[3] = "in";
  char pltitle[17] = "RLaB2: spwrite()";
  double dsize = 8.5;

  i = strlen (ime);
  OPENFILE (ime, &i, &unit);
  PSPLTM (&nrow, &ncol, &mode, m->ja, m->ia, pltitle, &ptitle, &dsize,
	  munt, &nlines, NULL, &unit);
  CLOSFILE (&unit);
  return 0;
}
